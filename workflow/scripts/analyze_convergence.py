#!/usr/bin/env python3
"""
Cross-family convergence analysis over the 100-bin normalised coordinate.

This replaces the standalone analysis_v4 script that produced the manuscript's
headline numbers but was never part of the workflow, so `snakemake` reproduced
none of them. Three substantive changes were made while porting it:

1. Site selection follows one stated rule. The previous version hardcoded nine
   site numbers for three datasets (Andes L, Ebola GP, Ebola L). Those literals
   matched no computable criterion: applying the documented "significant in >= 2
   of three site-level methods" rule yields zero sites for Andes L and Ebola L,
   and one of the hardcoded sites (Andes L 245) has FEL alpha=2.68, beta=0.0,
   i.e. maximal purifying selection.

2. Bin co-occurrence is tested against a null model. It was previously reported
   as a raw count of overlaps, so "three groups share bin 12" carried no p-value.

3. Domain enrichment conditions on sites the test could actually reach.
   Contrast-FEL returns a result only where both clades vary; counting untested
   positions as "tested and not selected" inflates the odds ratio, because
   testable sites are themselves concentrated in the N-terminus.

Both raw and FDR-corrected results are written so that the reported threshold is
an explicit choice rather than an artefact of which column was read.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

PROTEIN_ROLE = {
    "GnGc": "Entry",
    "G_protein": "Entry",
    "F_protein": "Entry_Helper",
    "GP": "Entry",
    "L_protein": "Replication",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--all-sites", nargs="+", required=True, dest="all_sites",
                   help="Per virus_group/protein *_all_sites.tsv files")
    p.add_argument("--config", default="config/config.yaml")
    p.add_argument("--outdir", required=True)
    p.add_argument("--pvalue", type=float, default=0.05)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--permutations", type=int, default=20000)
    p.add_argument("--gard-summaries", nargs="*", default=[], dest="gard",
                   help="Per dataset *_gard_summary.txt; sites closer than "
                        "--breakpoint-margin codons to a detected breakpoint are "
                        "flagged, and datasets whose GARD run failed are marked "
                        "unknown rather than clean")
    p.add_argument("--breakpoint-margin", type=int, default=10,
                   dest="bp_margin")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def group_metadata(cfg):
    viruses = cfg.get("viruses", {})
    meta = {}
    for group, spec in cfg.get("virus_groups", {}).items():
        cats = [viruses[v]["category"] for v in spec.get("viruses", [])
                if v in viruses]
        non_res = [c for c in cats if c != "Reservoir"]
        meta[group] = {
            "family": spec.get("family", "Unknown"),
            "category": non_res[0] if non_res else "Reservoir",
            "has_contrast": "Reservoir" in cats and bool(set(cats) - {"Reservoir"}),
        }
    return meta


def load_gard(paths):
    """
    Map (group, protein) -> {"status": ..., "breakpoints": [codon, ...]}.

    A failed GARD run yields status FAILED, never an empty breakpoint list, so a
    crash cannot be read downstream as "no recombination detected". This is the
    distinction the previous pipeline collapsed, which is how the manuscript came
    to state that screening had confirmed the sites were clean.
    """
    out = {}
    for path in paths:
        p = Path(path)
        group = p.parent.name
        protein = p.name.replace("_gard_summary.txt", "")
        try:
            rows = list(pd.read_csv(path, sep="\t").itertuples())
        except Exception:
            out[(group, protein)] = {"status": "UNREADABLE", "breakpoints": []}
            continue
        if not rows:
            out[(group, protein)] = {"status": "UNREADABLE", "breakpoints": []}
            continue
        r = rows[0]
        status = getattr(r, "status", "OK")
        raw = str(getattr(r, "breakpoint_positions", "None"))
        bps = []
        if raw not in ("None", "nan", "NA", ""):
            for tok in raw.split(","):
                try:
                    bps.append(int(int(tok) / 3))   # nucleotide -> codon
                except ValueError:
                    pass
        out[(group, protein)] = {"status": status, "breakpoints": bps}
    return out


def annotate_recombination(df, gard, margin):
    """Flag sites near a breakpoint, and mark datasets whose GARD run failed."""
    status, near = [], []
    for group, protein, site in zip(df.virus_group, df.protein, df.site):
        info = gard.get((group, protein))
        if info is None:
            status.append("NOT_RUN"); near.append(False); continue
        status.append(info["status"])
        near.append(any(abs(site - bp) <= margin for bp in info["breakpoints"]))
    df["gard_status"] = status
    df["near_breakpoint"] = near
    return df


def infer_group_protein(path):
    """<...>/04_selection/<virus_group>/<protein>_all_sites.tsv"""
    p = Path(path)
    return p.parent.name, p.name.replace("_all_sites.tsv", "")


def load_sites(paths, meta):
    """
    Read every site table and mark differential sites under one rule.

    Groups that have a reservoir arm use the branch contrast (Contrast-FEL).
    Groups without one cannot run it at all, so they fall back to agreement
    between site-level models. Which branch applies is decided by the group's
    composition in the config, never by the group's name.
    """
    frames = []
    for path in paths:
        group, protein = infer_group_protein(path)
        if group not in meta:
            print(f"WARNING: {group} is absent from the config; skipping {path}",
                  file=sys.stderr)
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            print(f"WARNING: {path} has no rows; skipping", file=sys.stderr)
            continue
        df["virus_group"] = group
        df["protein"] = protein
        df["role"] = PROTEIN_ROLE.get(protein, "Other")
        df["family"] = meta[group]["family"]
        df["category"] = meta[group]["category"]
        df["mode"] = "contrast" if meta[group]["has_contrast"] else "site_level"
        frames.append(df)

    if not frames:
        sys.exit("ERROR: no usable site tables were read.")
    return pd.concat(frames, ignore_index=True)


def mark_differential(df, pval, fdr):
    for col in ("cfel_p", "cfel_beta_h2h", "cfel_beta_ref", "cfel_q",
                "meme_q", "fel_q", "meme_p", "fel_p", "fubar_pp"):
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    consensus = df.get("consensus_pos", pd.Series(False, index=df.index))
    consensus = consensus.astype(str).str.lower().eq("true")

    # ONE rule for every group. Previously groups with a reservoir comparator
    # were scored by Contrast-FEL at raw p < 0.05 while Ebola and Sudan, which
    # have no comparator, were scored by consensus_pos. That is a 6.8-fold
    # difference in stringency (129 sites against 19 on the same six groups),
    # and it falls exactly along the H2H/comparator split, so any comparison of
    # differential-site burden between groups measured which rule the group
    # received rather than its biology.
    #
    # consensus_pos is the rule that every group can be scored by, so it is the
    # one applied. The branch contrast is retained as an annotation on the sites
    # it can speak to, never as an alternative entry route.
    contrast_support = (
        df["mode"].eq("contrast")
        & df["cfel_p"].lt(pval)
        & df["cfel_beta_h2h"].gt(df["cfel_beta_ref"])
    ).fillna(False)
    df["contrast_support"] = contrast_support
    df["differential"] = consensus.fillna(False)

    # "Testable" means the site models were fitted at this codon, which is true
    # wherever a p-value was produced. The previous definition counted every
    # Ebola and Sudan site as testable through a mode short-circuit while
    # discarding 8335 contrast sites whose Contrast-FEL p-value was exactly 1 --
    # p = 1 is a result, not a missing one.
    testable = pd.Series(False, index=df.index)
    for col in ("meme_p", "fel_p", "fubar_pp"):
        if col in df:
            testable = testable | pd.to_numeric(df[col], errors="coerce").notna()
    df["testable"] = testable

    # FDR on the same quantity the rule uses, so the correction is reachable for
    # every group. Keying it on cfel_q made it structurally impossible for the
    # two groups without a contrast, whose cfel_q is 1.0 by construction.
    q_cols = [c for c in ("meme_q", "fel_q") if c in df]
    if q_cols:
        q_min = pd.concat([pd.to_numeric(df[c], errors="coerce") for c in q_cols],
                          axis=1).min(axis=1)
        df["differential_fdr"] = df["differential"] & q_min.lt(fdr)
    else:
        df["differential_fdr"] = False
    return df


def add_bins(df):
    """
    Normalise each site onto a 100-bin coordinate using the true protein length.

    Length is taken as the maximum site index observed for that protein, which
    equals the alignment width in codons; this avoids the earlier bug where a
    partial first record set the length to 91 for a 2150-residue protein.
    """
    lengths = df.groupby(["virus_group", "protein"])["site"].max().rename("protein_length")
    df = df.merge(lengths, on=["virus_group", "protein"], how="left")
    df["bin_100"] = ((df["site"] - 1) / df["protein_length"] * 100).astype(int).clip(0, 99)
    return df


def positional_bias(df_diff, df_all, role, n_perm, rng):
    """
    Are this role's differential sites displaced toward the protein N-terminus?

    The bin co-occurrence test asks whether some single bin is shared by many
    groups. With one to six differential sites per group spread over 100 bins,
    two groups landing in the same bin is rare by construction, so that test has
    almost no power and a null result from it is not evidence of absence.

    This asks the same biological question -- do adaptive changes concentrate at
    one end of the protein -- using every site rather than only those that
    happen to coincide. The null permutes each group's sites within the
    positions that group actually has, so it preserves both how many sites a
    group contributes and where that group could place them.

    Entry and Replication are both reported. Replication is the control: a
    displacement seen in both roles is an artefact of the site-selection rule or
    of alignment geometry, not a property of entry proteins.
    """
    diff = df_diff[df_diff.role == role]
    pool = df_all[df_all.role == role]
    if len(diff) < 5 or pool.empty:
        return None

    rel_all = ((pool["site"] - 1)
               / (pool["protein_length"] - 1).clip(lower=1))
    rel_diff = ((diff["site"] - 1)
                / (diff["protein_length"] - 1).clip(lower=1))
    observed = float(rel_diff.mean())

    by_group = {g: sub.values for g, sub in rel_all.groupby(pool["virus_group"])}
    counts = diff.groupby("virus_group").size().to_dict()
    counts = {g: n for g, n in counts.items() if g in by_group}
    if not counts:
        return None

    hits = 0
    for _ in range(n_perm):
        drawn = []
        for g, n in counts.items():
            avail = by_group[g]
            drawn.extend(rng.choice(avail, size=min(n, len(avail)), replace=False))
        if float(np.mean(drawn)) <= observed:
            hits += 1

    return {
        "role": role,
        "n_sites": int(len(diff)),
        "n_groups": int(len(counts)),
        "n_families": int(diff["family"].nunique()) if "family" in diff else -1,
        "mean_relative_position": round(observed, 4),
        "background_mean_position": round(float(rel_all.mean()), 4),
        "permutations": n_perm,
        "p_value_n_terminal": (hits + 1) / (n_perm + 1),
    }


def hotspot_permutation(df_diff, df_all, role, n_perm, rng):
    """
    Probability that some bin in this role is shared by as many groups as the
    best observed bin, when each group's sites are placed at random over the
    bins the contrast could actually reach.

    Restricting the null to testable bins matters: testable positions are not
    uniform along the protein, so a uniform null would overstate significance.
    """
    diff = df_diff[df_diff.role == role]
    if diff.empty:
        return None

    observed = (diff.groupby(["bin_100", "virus_group"]).size()
                .reset_index(name="n")
                .groupby("bin_100")["virus_group"].nunique())
    if observed.empty:
        return None
    best = int(observed.max())

    pools, counts = {}, {}
    for group, sub in diff.groupby("virus_group"):
        testable = df_all[(df_all.virus_group == group)
                          & (df_all.role == role)
                          & (df_all.testable)]
        pool = testable["bin_100"].unique()
        if len(pool) == 0:
            pool = np.arange(100)
        pools[group] = pool
        counts[group] = sub["bin_100"].nunique()

    hits = 0
    for _ in range(n_perm):
        seen = {}
        for group, pool in pools.items():
            k = min(counts[group], len(pool))
            for b in rng.choice(pool, size=k, replace=False):
                seen[b] = seen.get(b, 0) + 1
        if seen and max(seen.values()) >= best:
            hits += 1
    return {"role": role, "max_groups_observed": best,
            "permutations": n_perm, "p_value": (hits + 1) / (n_perm + 1)}


def domain_enrichment(df_all, df_diff, domains, fdr):
    """
    Fisher's exact test per protein and domain, with the denominator limited to
    sites the contrast could reach.
    """
    records = []
    for (protein,), sub in df_all.groupby(["protein"]):
        key = protein.lower()
        if key not in domains:
            continue
        tested = sub[sub.testable]
        sel = df_diff[df_diff.protein == protein]
        for dom_name, (start, end) in domains[key].items():
            in_dom = tested[(tested.site >= start) & (tested.site <= end)]
            out_dom = tested[(tested.site < start) | (tested.site > end)]
            sel_in = len(sel[(sel.site >= start) & (sel.site <= end)])
            sel_out = len(sel) - sel_in
            table = [[sel_in, len(in_dom) - sel_in],
                     [sel_out, len(out_dom) - sel_out]]
            if min(len(in_dom), len(out_dom)) == 0:
                continue
            odds, p = fisher_exact(table)
            records.append({"protein": protein, "domain": dom_name,
                            "domain_start": start, "domain_end": end,
                            "selected_in_domain": sel_in,
                            "tested_in_domain": len(in_dom),
                            "selected_outside": sel_out,
                            "tested_outside": len(out_dom),
                            "odds_ratio": odds, "p_value": p})
    if not records:
        return pd.DataFrame()
    out = pd.DataFrame(records).sort_values("p_value")
    reject, q, _, _ = multipletests(out["p_value"], alpha=fdr, method="fdr_bh")
    out["q_value"] = q
    out["significant_fdr"] = reject
    return out


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    cfg = load_config(args.config)
    meta = group_metadata(cfg)
    domains = cfg.get("protein_domains", {})
    if not domains:
        print("WARNING: no protein_domains in config; domain enrichment skipped.",
              file=sys.stderr)

    df = load_sites(args.all_sites, meta)
    df = mark_differential(df, args.pvalue, args.fdr)
    df = add_bins(df)

    gard = load_gard(args.gard)
    df = annotate_recombination(df, gard, args.bp_margin)

    df_diff = df[df.differential].copy()
    df_diff.to_csv(Path(args.outdir) / "differential_sites.tsv",
                   sep="\t", index=False)

    # Reindexed over every (group, category, protein) actually analysed, so a
    # group with no differential site is a row of zero rather than an absence.
    universe = (df[["virus_group", "category", "protein"]]
                .drop_duplicates().sort_values(["virus_group", "protein"]))
    counts = (df_diff.groupby(["virus_group", "category", "protein"])
              .size().reset_index(name="n_differential"))
    counts = (universe.merge(counts, on=["virus_group", "category", "protein"],
                             how="left")
              .fillna({"n_differential": 0}))
    counts["n_differential"] = counts["n_differential"].astype(int)
    counts.to_csv(Path(args.outdir) / "differential_site_counts.tsv",
                  sep="\t", index=False)

    # Bin co-occurrence, with the null model the previous version lacked.
    pivot = (df_diff.groupby(["role", "bin_100", "virus_group"]).size()
             .reset_index(name="n")
             .pivot_table(index=["role", "bin_100"], columns="virus_group",
                          values="n", aggfunc="sum").fillna(0))
    if not pivot.empty:
        pivot["n_groups"] = (pivot > 0).sum(axis=1)
        hotspots = pivot[pivot.n_groups >= 2].sort_values("n_groups", ascending=False)
    else:
        hotspots = pd.DataFrame(columns=["role", "bin_100", "n_groups"])
    # Written on both paths. This is a declared rule output, and skipping it when
    # no site is differential failed the very last job in the DAG after every
    # expensive step had already succeeded.
    hotspots.to_csv(Path(args.outdir) / "hotspot_bins.tsv", sep="\t")

    perm = [r for r in (hotspot_permutation(df_diff, df, role, args.permutations, rng)
                        for role in sorted(df_diff.role.unique())) if r]
    # One test per functional role means the family is the set of roles, not each
    # role alone. Reporting three uncorrected p-values invites quoting whichever
    # is smallest: Replication at 0.035 becomes q = 0.10 over three tests, which
    # does not clear 0.05.
    if perm:
        raw = [r["p_value"] for r in perm]
        adj = multipletests(raw, method="fdr_bh")[1] if len(raw) > 1 else raw
        for r, q in zip(perm, adj):
            r["q_value_across_roles"] = float(q)
            r["n_roles_tested"] = len(perm)
    pd.DataFrame(perm).to_csv(Path(args.outdir) / "hotspot_permutation_test.tsv",
                              sep="\t", index=False)

    # Same hypothesis, tested with every site instead of only coincidences, and
    # with Replication as the control. Reported alongside the co-occurrence test,
    # never instead of it, and corrected across the roles examined.
    pos = [r for r in (positional_bias(df_diff, df, role, args.permutations, rng)
                       for role in sorted(df_diff.role.unique())) if r]
    if pos:
        raw = [r["p_value_n_terminal"] for r in pos]
        adj = multipletests(raw, method="fdr_bh")[1] if len(raw) > 1 else raw
        for r, q in zip(pos, adj):
            r["q_value_across_roles"] = float(q)
    pd.DataFrame(pos).to_csv(Path(args.outdir) / "positional_bias_test.tsv",
                             sep="\t", index=False)

    # Groups with no differential site must be reported as zero rather than be
    # absent: a missing row cannot be told apart from a group that was never
    # analysed, and Sudan_ebolavirus disappeared from both outputs that way.
    all_groups = sorted(df["virus_group"].dropna().unique())

    dom = domain_enrichment(df, df_diff, domains, args.fdr)
    if dom.empty:
        dom = pd.DataFrame(columns=["role", "protein", "domain", "n_differential",
                                    "n_tested", "odds_ratio", "p_value", "q_value"])
    # Also written unconditionally: a stale copy from an earlier run would
    # otherwise survive a rerun that found no enriched domain, and contradict
    # convergence_summary.json.
    dom.to_csv(Path(args.outdir) / "domain_enrichment.tsv", sep="\t", index=False)

    summary = {
        "positional_bias": pos,
        "n_sites_total": int(len(df)),
        "n_sites_testable": int(df.testable.sum()),
        "n_differential_raw_p": int(len(df_diff)),
        "n_differential_after_fdr": int(df.differential_fdr.sum()),
        "threshold_raw_p": args.pvalue,
        "groups": {g: int((df_diff.virus_group == g).sum()) for g in all_groups},
        "hotspot_bins": int(len(hotspots)),
        "permutation_tests": perm,
        "domains_significant_after_fdr": (
            int(dom.significant_fdr.sum()) if not dom.empty else 0),
        "recombination": {
            "datasets_screened_ok": int(sum(
                1 for v in gard.values() if v["status"] == "OK")),
            "datasets_failed_or_missing": int(sum(
                1 for v in gard.values() if v["status"] != "OK")),
            "differential_sites_near_breakpoint": int(df_diff.near_breakpoint.sum())
                if "near_breakpoint" in df_diff else 0,
            "breakpoint_margin_codons": args.bp_margin,
        },
    }
    with open(Path(args.outdir) / "convergence_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))
    if summary["n_differential_after_fdr"] == 0 and summary["n_differential_raw_p"]:
        print("\nNOTE: no differential site survives FDR correction. The counts "
              "above are uncorrected and must be reported as exploratory.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
