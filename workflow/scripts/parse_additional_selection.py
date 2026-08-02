"""
parse_additional_selection.py
==============================
Merges ALL HyPhy site-level selection results into one master TSV per protein.

Handles: MEME, FEL, FUBAR, SLAC, Contrast-FEL, PRIME
Plus: reads GARD summary to flag if recombination was detected.

Output columns (one row per codon site):
  site, meme_p, fel_p, fubar_pp, slac_ds, slac_dn, slac_dN_minus_dS, slac_p,
  cfell_beta_h2h, cfel_beta_ref, cfel_subs_fg, cfel_perm_p, cfel_p, cfel_qval, cfel_sig,
  prime_hydrophobicity_p, prime_isoelectric_point_p, prime_volume_p,
  prime_polarity_p, prime_charge_p, prime_composition_p,
  prime_sig_properties,
  n_methods_pos, n_methods_neg, consensus_pos, consensus_neg,
  recombination_warning, virus_group, protein
"""

import argparse
import json
import os
import sys

# ---------------------------------------------------------------------------
# Shared helpers (same pattern as parse_hyphy.py)
# ---------------------------------------------------------------------------

STATUS_SENTINELS = {}


def benjamini_hochberg(pvals):
    """
    BH q-values for a list of (key, p) pairs; returns {key: q}.

    Every site model here is applied once per codon -- 546 to 2350 tests per
    alignment -- and the pipeline reported raw p < 0.05 throughout. Across the
    study that is 75 MEME calls of which 2 survive correction, so the correction
    is not a formality.
    """
    usable = [(k, float(v)) for k, v in pvals
              if isinstance(v, (int, float)) and 0.0 <= float(v) <= 1.0]
    n = len(usable)
    if n == 0:
        return {}
    usable.sort(key=lambda kv: kv[1])
    q = {}
    prev = 1.0
    for i in range(n - 1, -1, -1):
        k, pv = usable[i]
        prev = min(prev, pv * n / (i + 1), 1.0)
        q[k] = prev
    return q


def load_json(path, label=""):
    """
    Load a HyPhy JSON, recording sentinel files instead of silently treating
    them as an analysis that ran and found nothing.

    The rules write {"status": "FAILED"} when HyPhy crashes and
    {"status": "NOT_APPLICABLE"} when a group has no branch contrast. Both used
    to arrive here as a dict without an "MLE" key, which every parser turns into
    an empty result -- indistinguishable from a clean negative. A crashed MEME
    run silently removed sites from the consensus while the rule stayed green.
    """
    key = label or os.path.basename(path or "?")
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        STATUS_SENTINELS[key] = "MISSING"
        return None
    try:
        with open(path) as f:
            data = json.load(f)
    except Exception as e:
        print(f"Warning: cannot read {path}: {e}", file=sys.stderr)
        STATUS_SENTINELS[key] = "UNREADABLE"
        return None
    if not data:
        STATUS_SENTINELS[key] = "EMPTY"
        return None
    if isinstance(data, dict) and "MLE" not in data and "status" in data:
        STATUS_SENTINELS[key] = str(data["status"])
        return None
    return data


def extract_rows(content):
    """HyPhy content is dict {'0': [rows...]} or plain list."""
    if isinstance(content, dict):
        key = list(content.keys())[0]
        return content[key]
    return content if isinstance(content, list) else []


def find_col(headers, keywords, default):
    """Search headers list (each item may be [name, desc]) for keywords."""
    for i, h in enumerate(headers):
        text = " ".join(h).lower() if isinstance(h, list) else str(h).lower()
        if any(k in text for k in keywords):
            return i
    return default


# ---------------------------------------------------------------------------
# Individual parsers
# ---------------------------------------------------------------------------

def parse_meme(data, pval_thresh):
    """Returns {site: meme_p} — p-value at index 6."""
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    rows = extract_rows(mle.get("content", {}))
    idx = find_col(headers, ["p-value", "pval"], 6)
    out = {}
    for i, row in enumerate(rows):
        try:
            out[i + 1] = float(row[idx]) if len(row) > idx else 1.0
        except (TypeError, ValueError):
            out[i + 1] = 1.0
    return out


def parse_fel(data, pval_thresh):
    """
    Returns (   {site: fel_p},   {sites where beta > alpha}   ).

    FEL's p-value is two-sided (beta != alpha), so significance alone does not
    indicate positive selection — purifying sites are equally significant. The
    second return value carries the direction, and callers must require it before
    counting a site toward positive-selection method agreement.
    """
    if not data or "MLE" not in data:
        return {}, set()
    mle = data["MLE"]
    headers = mle.get("headers", [])
    rows = extract_rows(mle.get("content", {}))
    idx = find_col(headers, ["p-value", "pval"], 4)
    a_idx = find_col(headers, ["alpha"], 0)
    b_idx = find_col(headers, ["beta"], 1)
    out = {}
    positive = set()
    for i, row in enumerate(rows):
        site = i + 1
        try:
            out[site] = float(row[idx]) if len(row) > idx else 1.0
        except (TypeError, ValueError):
            out[site] = 1.0
        try:
            if len(row) > max(a_idx, b_idx) and float(row[b_idx]) > float(row[a_idx]):
                positive.add(site)
        except (TypeError, ValueError):
            pass
    return out, positive


def parse_fubar(data, pp_thresh):
    """Returns {site: fubar_pp} — Prob[alpha<beta] at index 4."""
    if not data:
        return {}
    mle = data.get("MLE", data.get("posterior", {}))
    if not mle:
        return {}
    headers = mle.get("headers", [])
    rows = extract_rows(mle.get("content", {}))
    idx = find_col(headers, ["prob[alpha<beta]", "alpha<beta", "positive selection"], 4)
    out = {}
    for i, row in enumerate(rows):
        try:
            out[i + 1] = float(row[idx]) if len(row) > idx else 0.0
        except (TypeError, ValueError):
            out[i + 1] = 0.0
    return out


def parse_slac(data, pval_thresh):
    """
    SLAC MLE headers (actual order in HyPhy JSON):
      0: ES  (expected synonymous)
      1: EN  (expected non-synonymous)
      2: S   (observed synonymous)
      3: N   (observed non-synonymous)
      4: P[S] (expected proportion synonymous)
      5: dS
      6: dN
      7: dN-dS
      8: P[dN/dS > 1]  <- positive selection p-value
      9: P[dN/dS < 1]  <- negative selection p-value
     10: Total branch length
    Returns {site: {slac_ds, slac_dn, slac_dN_minus_dS, slac_p, slac_np}}

    CRITICAL: SLAC JSON nests per-site data as:
      MLE.content['0']['by-site']['AVERAGED']  -> list of per-site rows
    Do NOT use extract_rows() here — the outer dict has 'by-branch'/'by-site'
    string keys, not numeric keys.
    """
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    content = mle.get("content", {})
    partition = content.get("0", {})

    # Navigate to per-site data
    if not isinstance(partition, dict):
        return {}
    by_site = partition.get("by-site", {})
    if not isinstance(by_site, dict):
        return {}
    rows = by_site.get("AVERAGED", [])
    if not rows:
        return {}

    # Fixed column indices (verified against SLAC JSON headers)
    ds_idx  = 5   # dS
    dn_idx  = 6   # dN
    dif_idx = 7   # dN-dS
    pp_idx  = 8   # P[dN/dS > 1] — positive selection
    np_idx  = 9   # P[dN/dS < 1] — negative selection

    def _safe(v, default="NA"):
        if v is None:
            return default
        try:
            return float(v)
        except (TypeError, ValueError):
            return default

    out = {}
    for i, row in enumerate(rows):
        try:
            out[i + 1] = {
                "slac_ds":   _safe(row[ds_idx])  if len(row) > ds_idx  else "NA",
                "slac_dn":   _safe(row[dn_idx])  if len(row) > dn_idx  else "NA",
                "slac_dN_minus_dS": _safe(row[dif_idx]) if len(row) > dif_idx else "NA",
                "slac_p":    _safe(row[pp_idx], default=1.0) if len(row) > pp_idx else 1.0,
                "slac_np":   _safe(row[np_idx], default=1.0) if len(row) > np_idx else 1.0,
            }
        except (IndexError, TypeError, ValueError):
            out[i + 1] = {"slac_ds": "NA", "slac_dn": "NA", "slac_dN_minus_dS": "NA",
                          "slac_p": 1.0, "slac_np": 1.0}
    return out



def _num(row, idx, default="NA"):
    """Read a numeric cell, returning `default` when absent or unparseable."""
    try:
        return float(row[idx]) if len(row) > idx else default
    except (TypeError, ValueError):
        return default


def _permutation_ok(perm, threshold):
    """
    True unless HyPhy's permutation test contradicts the likelihood-ratio call.

    A value of -1 means the site was never permuted, which is not evidence
    either way, so the q-value stands alone. Any real value at or above the
    threshold is a refutation.
    """
    try:
        value = float(perm)
    except (TypeError, ValueError):
        return True
    if value < 0:
        return True
    return value < threshold


def parse_contrast_fel(data, pval_thresh, fdr_thresh, require_perm=True):
    """
    Contrast-FEL MLE headers (actual HyPhy output):
      0: alpha               (dS / synonymous rate)
      1: beta (Foreground)   (dN in H2H / test branches)   <-- FOREGROUND
      2: beta (background)   (dN in reservoir / background) <-- BACKGROUND
      3: subs (Foreground)   (substitution count)
      4: P-value (overall)   <-- main p-value
      5: Q-value (overall)   <-- FDR-corrected
      6: Permutation p-value
      7: Total branch length

    NOTE: Foreground=index 1, Background=index 2.
    The original parser had these SWAPPED (ref_idx=1,test_idx=2 with wrong keywords).
    Also find_col was matching wrong columns for p-value.
    """
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    rows = extract_rows(mle.get("content", {}))

    # Fixed indices verified against actual Contrast-FEL JSON headers
    fg_idx = 1   # beta (Foreground) = H2H branches
    bg_idx = 2   # beta (background) = reservoir branches
    p_idx  = 4   # P-value (overall)
    q_idx  = 5   # Q-value (overall)
    subs_idx = 3  # subs (Foreground): substitutions supporting the contrast
    perm_idx = 6  # Permutation p-value; HyPhy writes -1 when not evaluated

    # Dynamic override if header labels differ
    for i, h in enumerate(headers):
        txt = " ".join(h).lower() if isinstance(h, list) else str(h).lower()
        if "foreground" in txt and "beta" in txt:
            fg_idx = i
        elif "background" in txt and "beta" in txt:
            bg_idx = i
        elif "p-value (overall)" in txt or ("p-value" in txt and "overall" in txt):
            p_idx = i
        elif "q-value (overall)" in txt or ("q-value" in txt and "overall" in txt):
            q_idx = i
        elif "subs" in txt and "foreground" in txt:
            subs_idx = i
        elif "permutation" in txt:
            perm_idx = i

    out = {}
    for i, row in enumerate(rows):
        try:
            p = float(row[p_idx]) if len(row) > p_idx else 1.0
            q = float(row[q_idx]) if len(row) > q_idx else 1.0
            beta_fg = round(float(row[fg_idx]), 6) if len(row) > fg_idx else "NA"
            beta_bg = round(float(row[bg_idx]), 6) if len(row) > bg_idx else "NA"
            # Substitution support and the permutation p-value were previously
            # discarded. Without them a "differential" site can rest on a single
            # substitution against a background rate estimated at zero, and the
            # -1 sentinel that marks an unevaluated permutation is invisible.
            subs = _num(row, subs_idx, default="NA")
            perm = _num(row, perm_idx, default="NA")
            out[i + 1] = {
                "cfel_beta_h2h":      beta_fg,
                "cfel_beta_ref":      beta_bg,
                "cfel_subs_fg":       subs,
                "cfel_perm_p":        perm,
                "cfel_p":             round(p, 6),
                "cfel_q":             round(q, 6),
                # HyPhy runs a permutation test over branch assignments for
                # sites whose LRT is significant, precisely because the LRT can
                # be driven by one substitution against a background rate
                # estimated at zero. Calling a site significant on the q-value
                # alone discards that check: of the 11 sites that reached
                # q < 0.20 in this study, none reached p < 0.05 under
                # permutation and six sat at exactly 1.0. When HyPhy did not
                # permute a site it writes -1, and that is not evidence either
                # way, so the q-value stands alone there.
                "cfel_sig":           (q < fdr_thresh
                                       and (not require_perm
                                            or _permutation_ok(perm, pval_thresh))),
                "cfel_h2h_stronger":  (
                    isinstance(beta_fg, float) and isinstance(beta_bg, float)
                    and beta_fg > beta_bg and q < fdr_thresh
                ),
            }
        except (TypeError, ValueError):
            out[i + 1] = {"cfel_beta_h2h": "NA", "cfel_beta_ref": "NA",
                          "cfel_subs_fg": "NA", "cfel_perm_p": "NA",
                          "cfel_p": 1.0, "cfel_q": 1.0,
                          "cfel_sig": False, "cfel_h2h_stronger": False}
    return out


def parse_prime(data, pval_thresh):
    """
    PRIME MLE headers (actual HyPhy output — 3 properties shown, more possible):
      0:  alpha (dS)
      1:  beta  (overall non-syn rate)
      2:  FEL alpha
      3:  FEL beta
      4:  Total branch length
      5:  # subs
      6:  # aa
      7:  PRIME LogL
      8:  FEL LogL
      9:  p-value  (omnibus — ANY property important)   <-- OVERALL
      10: q-value
      11: R (redundancy)
      12: lambda1  (effect size, property 1)
      13: p1       (p-value, property 1)                <-- PROPERTY p-values
      14: LogL1
      15: lambda2  (effect size, property 2)
      16: p2
      17: LogL2
      18: lambda3
      19: p3
      20: LogL3
      ... (up to 5 properties: pattern is idx 12+3k for lambda, 13+3k for p)

    Property names come from header labels (e.g. Hydrophobicity_KyteDoolittle,
    Isoelectric_Point_pI, Volume_Angstrom3, ...). We read them dynamically.

    CRITICAL: overall p-value is index 9 (not 2 as previously assumed).
    Property p-values follow pattern: p_k = index 13 + 3*(k-1) for k=1,2,3,...
    """
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    rows = extract_rows(mle.get("content", {}))

    # ── Find overall p-value index ────────────────────────────────────────────
    overall_p_idx = 9   # default from verified JSON
    for i, h in enumerate(headers):
        txt = " ".join(h).lower() if isinstance(h, list) else str(h).lower()
        if ("p-value" in txt or "p value" in txt) and ("omni" in txt or "overall" in txt or "any" in txt):
            overall_p_idx = i
            break

    # ── Discover property p-value columns dynamically ─────────────────────────
    # Pattern: headers contain p-values for each property as 'p1', 'p2', ... or
    # labeled with the property name. We scan for headers where label matches 'pN'
    # or where description contains a recognisable property name.
    PROP_KEYWORDS = {
        "hydrophobicity": "hydrophobicity",
        "kytedoolittle":  "hydrophobicity",
        "isoelectric":    "isoelectric_point",
        "volume":         "volume",
        "polarity":       "polarity",
        "charge":         "charge",
        "composition":    "composition",
        "size":           "size",
        "flexibility":    "flexibility",
    }
    prop_cols = {}   # {friendly_name: col_idx}
    for i, h in enumerate(headers):
        label = h[0].lower() if isinstance(h, list) else str(h).lower()
        desc  = h[1].lower() if isinstance(h, list) and len(h) > 1 else ""
        # Must be a p-value column (label starts with 'p' followed by digit, or desc says p-value)
        is_pval_col = (
            (len(label) <= 3 and label.startswith("p") and label[1:].isdigit())
            or ("p-value" in desc and "non-zero" in desc)
        )
        if not is_pval_col:
            continue
        # Identify which property
        for kw, fname in PROP_KEYWORDS.items():
            if kw in desc:
                if fname not in prop_cols:   # take first match per property
                    prop_cols[fname] = i
                break

    # Standardise output column names to 5 canonical names
    CANONICAL = [
        "volume", "polarity", "charge", "hydrophobicity",
        "isoelectric_point", "composition", "size", "flexibility"
    ]

    out = {}
    for i, row in enumerate(rows):
        try:
            overall_p = float(row[overall_p_idx]) if len(row) > overall_p_idx else 1.0
            entry = {"prime_overall_p": round(overall_p, 6)}
            sig_props = []
            # PRIME's own omnibus test decides whether the per-property tests
            # should be looked at at all. Reading the properties regardless
            # flagged 25 sites across the study whose omnibus q was 0.29 or
            # worse, 23 of them at exactly 1.0.
            omnibus_ok = overall_p < pval_thresh
            for prop, col in prop_cols.items():
                p = float(row[col]) if len(row) > col else 1.0
                entry[f"prime_{prop}_p"] = round(p, 6)
                if p < pval_thresh:
                    sig_props.append(prop)
            # Fill canonical columns that weren't found
            for prop in CANONICAL:
                if f"prime_{prop}_p" not in entry:
                    entry[f"prime_{prop}_p"] = "NA"
            entry["prime_sig_properties"] = (
                "|".join(sig_props) if (sig_props and omnibus_ok) else "None")
            out[i + 1] = entry
        except (TypeError, ValueError):
            entry = {"prime_overall_p": 1.0, "prime_sig_properties": "None"}
            for prop in CANONICAL:
                entry[f"prime_{prop}_p"] = "NA"
            out[i + 1] = entry
    return out


def read_gard_warning(gard_summary_path):
    """
    Tri-state recombination flag: "Yes", "No" or "Unknown".

    A boolean cannot carry the difference between an alignment GARD screened and
    found clean and one it never screened. Puumala L is on the skip list and its
    summary says status=SKIPPED with recombination_detected=NA; returning False
    for that reported it as free of recombination, which is exactly what the
    config comment forbids.
    """
    if not gard_summary_path or not os.path.exists(gard_summary_path):
        return "Unknown"
    try:
        with open(gard_summary_path) as f:
            header = f.readline().rstrip("\n").split("\t")
            line = f.readline().rstrip("\n")
            if not line:
                return "Unknown"
            # Look the column up by name. This used to index position 2, which
            # stopped being recombination_detected the moment a status column was
            # added in front of it, so the test compared "OK" against "true" and
            # every alignment reported no recombination.
            row = dict(zip(header, line.split("\t")))
            if row.get("status", "").strip().upper() != "OK":
                return "Unknown"
            detected = row.get("recombination_detected", "").strip().lower()
            if detected == "true":
                return "Yes"
            if detected == "false":
                return "No"
            return "Unknown"
    except Exception:
        return "Unknown"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Merge all HyPhy site-level selection results into master TSV"
    )
    p.add_argument("--meme",         required=True)
    p.add_argument("--fel",          required=True)
    p.add_argument("--fubar",        required=True)
    p.add_argument("--slac",         required=True)
    p.add_argument("--contrast-fel", required=True, dest="contrast_fel")
    p.add_argument("--prime",        required=True)
    p.add_argument("--gard",         required=True,
                   help="GARD summary TSV (from parse_gard.py)")
    p.add_argument("--pvalue",       type=float, default=0.05)
    p.add_argument("--fubar-pp",     type=float, default=0.90, dest="fubar_pp")
    p.add_argument("--contrast-fdr", type=float, default=0.20, dest="contrast_fdr")
    p.add_argument("--require-permutation", dest="require_permutation",
                   default="yes", choices=["yes", "no"],
                   type=lambda v: str(v).strip().lower(),
                   help="Whether a Contrast-FEL call must also survive HyPhy's "
                        "permutation test. The manuscript as submitted states "
                        "only the FDR criterion, so 'no' reproduces it; 'yes' "
                        "adds the check HyPhy computes for exactly this purpose.")
    p.add_argument("--min-methods", type=int, default=2, dest="min_methods",
                   help="Number of site-level methods that must agree before a "
                        "site is called consensus. The manuscript states 2.")
    p.add_argument("--virus-group",  required=True, dest="virus_group")
    p.add_argument("--protein",      required=True)
    p.add_argument("--out",          required=True)
    return p.parse_args()


def main():
    args = parse_args()
    args.require_permutation = (args.require_permutation == "yes")

    meme_res   = parse_meme(load_json(args.meme, "meme"),   args.pvalue)
    fel_res, fel_pos = parse_fel(load_json(args.fel, "fel"), args.pvalue)
    fubar_res  = parse_fubar(load_json(args.fubar, "fubar"), args.fubar_pp)
    slac_res   = parse_slac(load_json(args.slac, "slac"),   args.pvalue)
    cfel_res   = parse_contrast_fel(load_json(args.contrast_fel, "contrast_fel"),
                                 args.pvalue, args.contrast_fdr,
                                 args.require_permutation)
    prime_res  = parse_prime(load_json(args.prime, "prime"), args.pvalue)
    gard_warn  = read_gard_warning(args.gard)

    # Determine total sites
    all_dicts = [meme_res, fel_res, fubar_res, slac_res, cfel_res, prime_res]
    total = max((max(d.keys(), default=0) for d in all_dicts if d), default=0)

    if STATUS_SENTINELS:
        for name, status in sorted(STATUS_SENTINELS.items()):
            print(f"[{args.virus_group}/{args.protein}] {name}: {status}",
                  file=sys.stderr)
    # NOT_APPLICABLE is a legitimate state for the two branch-aware tests in a
    # group with no reservoir comparator. Anything else, on any of the four site
    # models, means an analysis did not run -- and a missing site model silently
    # removes sites from the consensus while leaving the rule green.
    fatal = {n: st for n, st in STATUS_SENTINELS.items()
             if not (n in ("contrast_fel",) and st == "NOT_APPLICABLE")}
    if fatal:
        print(f"ERROR: {args.virus_group}/{args.protein}: these analyses produced "
              f"no usable output: {fatal}. Refusing to write a table that would "
              f"read as a clean negative.", file=sys.stderr)
        sys.exit(1)

    if total == 0:
        print("ERROR: no site data in any input file; refusing to write an empty "
              "table that downstream steps would treat as a real result.",
              file=sys.stderr)
        sys.exit(1)

    # FDR across codons, per alignment and per method.
    meme_q = benjamini_hochberg(list(meme_res.items()))
    fel_q_map = benjamini_hochberg([(k, v) for k, v in fel_res.items()])
    slac_q = benjamini_hochberg([(k, v["slac_p"]) for k, v in slac_res.items()
                                 if isinstance(v.get("slac_p"), float)])
    prime_q = benjamini_hochberg([(k, v["prime_overall_p"]) for k, v in prime_res.items()
                                  if isinstance(v.get("prime_overall_p"), float)])

    header = [
        "site",
        "meme_p", "meme_q", "fel_p", "fel_q", "fubar_pp",
        "slac_ds", "slac_dn", "slac_dN_minus_dS", "slac_p", "slac_q", "slac_np",
        "cfel_beta_h2h", "cfel_beta_ref", "cfel_subs_fg", "cfel_perm_p",
        "cfel_p", "cfel_q", "cfel_sig",
        "prime_overall_p", "prime_overall_q",
        "prime_hydrophobicity_p", "prime_isoelectric_point_p", "prime_volume_p",
        "prime_polarity_p", "prime_charge_p", "prime_composition_p",
        "prime_sig_properties",
        "n_methods_pos",  # how many methods detect positive selection
        "n_methods_neg",  # how many methods detect negative selection
        "episodic_only",  # bool: MEME-significant with no pervasive support
        "consensus_pos",  # bool: >= 2 methods agree on positive selection
        "consensus_neg",  # bool: >= 2 methods agree on negative selection
        "cfel_h2h_stronger",  # bool: H2H stronger than Reservoir at this site
        "recombination_warning",
        "virus_group", "protein",
    ]

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    n_pos = 0
    n_cfell_sig = 0

    with open(args.out, "w") as f:
        f.write("\t".join(header) + "\n")

        for site in range(1, total + 1):
            mp  = meme_res.get(site, 1.0)
            fp  = fel_res.get(site, 1.0)
            fup = fubar_res.get(site, 0.0)
            sl  = slac_res.get(site, {"slac_ds": "NA", "slac_dn": "NA",
                                       "slac_dN_minus_dS": "NA", "slac_p": 1.0, "slac_np": 1.0})
            cf  = cfel_res.get(site, {"cfel_beta_h2h": "NA", "cfel_beta_ref": "NA",
                                       "cfel_subs_fg": "NA", "cfel_perm_p": "NA",
                                       "cfel_p": 1.0, "cfel_q": 1.0, "cfel_sig": False})
            pr  = prime_res.get(site, {
                "prime_overall_p": 1.0,
                "prime_hydrophobicity_p": "NA", "prime_isoelectric_point_p": "NA",
                "prime_volume_p": "NA", "prime_polarity_p": "NA",
                "prime_charge_p": "NA", "prime_composition_p": "NA",
                "prime_sig_properties": "None",
            })

            # Count methods detecting positive selection. FEL requires beta > alpha:
            # its p-value is two-sided, so without the direction check purifying
            # sites are counted as positively selected.
            fel_is_pos = site in fel_pos
            pos_methods = sum([
                1 if mp  < args.pvalue else 0,    # MEME
                1 if (fp < args.pvalue and fel_is_pos) else 0,    # FEL
                1 if fup >= args.fubar_pp else 0, # FUBAR
                1 if (isinstance(sl["slac_p"], float) and sl["slac_p"] < args.pvalue) else 0,  # SLAC
            ])
            # Episodic selection is a separate finding, not a weaker version of
            # pervasive selection. MEME tests whether beta > alpha on a SUBSET of
            # branches; FEL and FUBAR both test whether it holds across the whole
            # tree. A genuinely episodic site is therefore MEME-significant and
            # FEL/FUBAR-negative by design, and requiring method agreement removes
            # exactly that category. Flagging it separately keeps consensus_pos
            # meaning "two independent methods agree" without discarding the sites
            # that motivated lowering the threshold in the first place.
            pervasive_hit = ((fp < args.pvalue and fel_is_pos)
                             or fup >= args.fubar_pp)
            episodic_only = (mp < args.pvalue) and not pervasive_hit

            # Count methods detecting negative selection (purifying)
            neg_methods = sum([
                1 if (isinstance(sl["slac_np"], float) and sl["slac_np"] < args.pvalue) else 0,
                1 if (isinstance(fp, float) and fp < args.pvalue and not fel_is_pos) else 0,
            ])

            # H2H stronger than reservoir at this site?
            h2h_stronger = False
            if cf["cfel_sig"] and cf["cfel_beta_h2h"] != "NA" and cf["cfel_beta_ref"] != "NA":
                try:
                    h2h_stronger = float(cf["cfel_beta_h2h"]) > float(cf["cfel_beta_ref"])
                except (TypeError, ValueError):
                    pass

            if pos_methods >= args.min_methods:
                n_pos += 1
            if cf["cfel_sig"]:
                n_cfell_sig += 1

            # %.6g rather than round(x, 6): the strongest MEME hit in the study
            # is p = 5.1e-10 and rounding wrote it as 0.0, erasing the magnitude
            # of the one result that survives correction.
            def g(x):
                return f"{x:.6g}" if isinstance(x, float) else x

            # HyPhy writes -1 in the permutation column for sites it did not
            # permute. Left as a number, any `perm_p < 0.05` filter downstream
            # treats those 15536 sites as maximally significant.
            perm_raw = cf["cfel_perm_p"]
            try:
                perm_out = "NA" if float(perm_raw) < 0 else g(float(perm_raw))
            except (TypeError, ValueError):
                perm_out = "NA"

            row = [
                site,
                g(mp), g(meme_q.get(site, 1.0)),
                g(fp), g(fel_q_map.get(site, 1.0)),
                g(fup),
                sl["slac_ds"], sl["slac_dn"], sl["slac_dN_minus_dS"],
                g(sl["slac_p"]), g(slac_q.get(site, 1.0)), g(sl["slac_np"]),
                cf["cfel_beta_h2h"], cf["cfel_beta_ref"],
                cf["cfel_subs_fg"], perm_out,
                cf["cfel_p"], cf["cfel_q"], cf["cfel_sig"],
                pr["prime_overall_p"], g(prime_q.get(site, 1.0)),
                pr["prime_hydrophobicity_p"], pr["prime_isoelectric_point_p"],
                pr["prime_volume_p"], pr["prime_polarity_p"],
                pr["prime_charge_p"], pr["prime_composition_p"],
                pr["prime_sig_properties"],
                pos_methods, neg_methods,
                episodic_only,
                pos_methods >= args.min_methods,  # consensus_pos
                neg_methods >= args.min_methods,  # consensus_neg
                h2h_stronger,
                gard_warn,
                args.virus_group, args.protein,
            ]
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"\n=== All-Site Selection Table: {args.virus_group} / {args.protein} ===")
    print(f"  Total sites     : {total}")
    print(f"  Consensus pos   : {n_pos} sites (>=2 methods agree)")
    print(f"  Contrast-FEL sig: {n_cfell_sig} sites (H2H vs Reservoir, FDR < {args.contrast_fdr})")
    print(f"  Recombination warning: {gard_warn}")
    print(f"  Output: {args.out}")


if __name__ == "__main__":
    main()
