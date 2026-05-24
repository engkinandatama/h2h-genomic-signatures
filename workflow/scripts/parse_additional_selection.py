"""
parse_additional_selection.py
==============================
Merges ALL HyPhy site-level selection results into one master TSV per protein.

Handles: MEME, FEL, FUBAR, SLAC, Contrast-FEL, PRIME
Plus: reads GARD summary to flag if recombination was detected.

Output columns (one row per codon site):
  site, meme_p, fel_p, fubar_pp, slac_ds, slac_dn, slac_dndS, slac_p,
  cfell_beta_h2h, cfel_beta_ref, cfel_p, cfel_qval, cfel_sig,
  prime_volume_p, prime_polarity_p, prime_charge_p, prime_hydrophobicity_p, prime_composition_p,
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

def load_json(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        with open(path) as f:
            data = json.load(f)
        return data if data else None
    except Exception as e:
        print(f"Warning: cannot read {path}: {e}", file=sys.stderr)
        return None


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
    """Returns {site: fel_p} — p-value at index 4."""
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    rows = extract_rows(mle.get("content", {}))
    idx = find_col(headers, ["p-value", "pval"], 4)
    out = {}
    for i, row in enumerate(rows):
        try:
            out[i + 1] = float(row[idx]) if len(row) > idx else 1.0
        except (TypeError, ValueError):
            out[i + 1] = 1.0
    return out


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
    SLAC MLE headers (typical order):
      0: ES (expected synonymous)
      1: EN (expected non-synonymous)
      2: S  (observed synonymous)
      3: N  (observed non-synonymous)
      4: dS
      5: dN
      6: dN-dS
      7: p-value (dN > dS)  ← positive selection
      8: p-value (dS > dN)  ← negative selection
    Returns {site: {ds, dn, dndS, p_pos, p_neg}}
    """
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    rows = extract_rows(mle.get("content", {}))

    ds_idx  = find_col(headers, ["ds", "synonymous rate"], 4)
    dn_idx  = find_col(headers, ["dn", "non-synonymous rate"], 5)
    dif_idx = find_col(headers, ["dn-ds", "dndS", "dn - ds"], 6)
    pp_idx  = find_col(headers, ["p-value (dn>ds)", "p-value (dn > ds)", "positive"], 7)
    np_idx  = find_col(headers, ["p-value (ds>dn)", "p-value (ds > dn)", "negative"], 8)

    out = {}
    for i, row in enumerate(rows):
        try:
            out[i + 1] = {
                "slac_ds":   float(row[ds_idx])  if len(row) > ds_idx  else "NA",
                "slac_dn":   float(row[dn_idx])  if len(row) > dn_idx  else "NA",
                "slac_dndS": float(row[dif_idx]) if len(row) > dif_idx else "NA",
                "slac_p":    float(row[pp_idx])  if len(row) > pp_idx  else 1.0,
                "slac_np":   float(row[np_idx])  if len(row) > np_idx  else 1.0,
            }
        except (TypeError, ValueError):
            out[i + 1] = {"slac_ds": "NA", "slac_dn": "NA", "slac_dndS": "NA",
                          "slac_p": 1.0, "slac_np": 1.0}
    return out


def parse_contrast_fel(data, pval_thresh, fdr_thresh):
    """
    Contrast-FEL MLE headers (typical order):
      0: alpha (dS)
      1: beta_ref (dN in reference/background)
      2: beta_test (dN in test/foreground = H2H)
      3: LRT
      4: p-value
      5: q-value (Benjamini-Hochberg corrected)
    Returns {site: {cfel_beta_h2h, cfel_beta_ref, cfel_p, cfel_q, cfel_sig}}
    """
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    rows = extract_rows(mle.get("content", {}))

    # Try to find foreground (test) and background (reference) columns
    # In Contrast-FEL: beta columns labeled as "Test" or "Reference"
    ref_idx  = find_col(headers, ["beta reference", "beta-", "background", "reference"], 1)
    test_idx = find_col(headers, ["beta test", "beta+", "foreground", "test"], 2)
    p_idx    = find_col(headers, ["p-value", "pval"], 4)
    q_idx    = find_col(headers, ["q-value", "fdr", "benjamini"], 5)

    out = {}
    for i, row in enumerate(rows):
        try:
            p = float(row[p_idx])  if len(row) > p_idx  else 1.0
            q = float(row[q_idx])  if len(row) > q_idx  else 1.0
            out[i + 1] = {
                "cfel_beta_h2h": round(float(row[test_idx]), 6) if len(row) > test_idx else "NA",
                "cfel_beta_ref": round(float(row[ref_idx]),  6) if len(row) > ref_idx  else "NA",
                "cfel_p":        round(p, 6),
                "cfel_q":        round(q, 6),
                "cfel_sig":      q < fdr_thresh,
            }
        except (TypeError, ValueError):
            out[i + 1] = {"cfel_beta_h2h": "NA", "cfel_beta_ref": "NA",
                          "cfel_p": 1.0, "cfel_q": 1.0, "cfel_sig": False}
    return out


def parse_prime(data, pval_thresh):
    """
    PRIME MLE headers (typical):
      0: alpha (dS)
      1: LRT (overall)
      2: p-value (overall)
      Then pairs of [LRT, p-value] for each property:
        Volume, Polarity, Charge, Hydrophobicity, Composition
    Returns {site: {prime_overall_p, prime_volume_p, prime_polarity_p,
                    prime_charge_p, prime_hydrophobicity_p, prime_composition_p,
                    prime_sig_properties}}
    """
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    rows = extract_rows(mle.get("content", {}))

    # Map property names to column indices by scanning headers
    prop_cols = {}
    for prop in ["volume", "polarity", "charge", "hydrophobicity", "composition"]:
        for i, h in enumerate(headers):
            text = " ".join(h).lower() if isinstance(h, list) else str(h).lower()
            if prop in text and "p-value" in text:
                prop_cols[prop] = i
                break

    overall_p_idx = find_col(headers, ["p-value", "overall p", "total"], 2)

    out = {}
    for i, row in enumerate(rows):
        try:
            entry = {
                "prime_overall_p": round(float(row[overall_p_idx]), 6)
                                   if len(row) > overall_p_idx else 1.0,
            }
            sig_props = []
            for prop, col in prop_cols.items():
                p = float(row[col]) if len(row) > col else 1.0
                entry[f"prime_{prop}_p"] = round(p, 6)
                if p < pval_thresh:
                    sig_props.append(prop)
            entry["prime_sig_properties"] = "|".join(sig_props) if sig_props else "None"
            # Fill missing properties
            for prop in ["volume", "polarity", "charge", "hydrophobicity", "composition"]:
                if f"prime_{prop}_p" not in entry:
                    entry[f"prime_{prop}_p"] = "NA"
            out[i + 1] = entry
        except (TypeError, ValueError):
            out[i + 1] = {
                "prime_overall_p": 1.0,
                "prime_volume_p": "NA", "prime_polarity_p": "NA",
                "prime_charge_p": "NA", "prime_hydrophobicity_p": "NA",
                "prime_composition_p": "NA",
                "prime_sig_properties": "None",
            }
    return out


def read_gard_warning(gard_summary_path):
    """Read GARD summary TSV and return True if recombination was detected."""
    if not gard_summary_path or not os.path.exists(gard_summary_path):
        return False
    try:
        with open(gard_summary_path) as f:
            f.readline()  # skip header
            line = f.readline().strip()
            if not line:
                return False
            parts = line.split("\t")
            # recombination_detected is column index 2
            return parts[2].strip().lower() == "true" if len(parts) > 2 else False
    except Exception:
        return False


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
    p.add_argument("--virus-group",  required=True, dest="virus_group")
    p.add_argument("--protein",      required=True)
    p.add_argument("--out",          required=True)
    return p.parse_args()


def main():
    args = parse_args()

    meme_res   = parse_meme(load_json(args.meme),   args.pvalue)
    fel_res    = parse_fel(load_json(args.fel),     args.pvalue)
    fubar_res  = parse_fubar(load_json(args.fubar), args.fubar_pp)
    slac_res   = parse_slac(load_json(args.slac),   args.pvalue)
    cfel_res   = parse_contrast_fel(load_json(args.contrast_fel), args.pvalue, args.contrast_fdr)
    prime_res  = parse_prime(load_json(args.prime), args.pvalue)
    gard_warn  = read_gard_warning(args.gard)

    # Determine total sites
    all_dicts = [meme_res, fel_res, fubar_res, slac_res, cfel_res, prime_res]
    total = max((max(d.keys(), default=0) for d in all_dicts if d), default=0)

    if total == 0:
        print("Warning: No site data found in any of the input files.", file=sys.stderr)

    header = [
        "site",
        "meme_p", "fel_p", "fubar_pp",
        "slac_ds", "slac_dn", "slac_dndS", "slac_p", "slac_np",
        "cfel_beta_h2h", "cfel_beta_ref", "cfel_p", "cfel_q", "cfel_sig",
        "prime_overall_p",
        "prime_volume_p", "prime_polarity_p", "prime_charge_p",
        "prime_hydrophobicity_p", "prime_composition_p",
        "prime_sig_properties",
        "n_methods_pos",  # how many methods detect positive selection
        "n_methods_neg",  # how many methods detect negative selection
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
                                       "slac_dndS": "NA", "slac_p": 1.0, "slac_np": 1.0})
            cf  = cfel_res.get(site, {"cfel_beta_h2h": "NA", "cfel_beta_ref": "NA",
                                       "cfel_p": 1.0, "cfel_q": 1.0, "cfel_sig": False})
            pr  = prime_res.get(site, {
                "prime_overall_p": 1.0,
                "prime_volume_p": "NA", "prime_polarity_p": "NA",
                "prime_charge_p": "NA", "prime_hydrophobicity_p": "NA",
                "prime_composition_p": "NA",
                "prime_sig_properties": "None",
            })

            # Count methods detecting positive selection
            pos_methods = sum([
                1 if mp  < args.pvalue else 0,    # MEME
                1 if fp  < args.pvalue else 0,    # FEL
                1 if fup >= args.fubar_pp else 0, # FUBAR
                1 if (isinstance(sl["slac_p"], float) and sl["slac_p"] < args.pvalue) else 0,  # SLAC
            ])
            # Count methods detecting negative selection (purifying)
            neg_methods = sum([
                1 if (isinstance(sl["slac_np"], float) and sl["slac_np"] < args.pvalue) else 0,
                1 if (isinstance(fp, float) and fp < args.pvalue and
                      "NA" not in str(sl["slac_dndS"]) and
                      isinstance(sl["slac_dndS"], float) and sl["slac_dndS"] < 0) else 0,
            ])

            # H2H stronger than reservoir at this site?
            h2h_stronger = False
            if cf["cfel_sig"] and cf["cfel_beta_h2h"] != "NA" and cf["cfel_beta_ref"] != "NA":
                try:
                    h2h_stronger = float(cf["cfel_beta_h2h"]) > float(cf["cfel_beta_ref"])
                except (TypeError, ValueError):
                    pass

            if pos_methods >= 2:
                n_pos += 1
            if cf["cfel_sig"]:
                n_cfell_sig += 1

            row = [
                site,
                round(mp,  6), round(fp, 6), round(fup, 6),
                sl["slac_ds"], sl["slac_dn"], sl["slac_dndS"],
                sl["slac_p"], sl["slac_np"],
                cf["cfel_beta_h2h"], cf["cfel_beta_ref"],
                cf["cfel_p"], cf["cfel_q"], cf["cfel_sig"],
                pr["prime_overall_p"],
                pr["prime_volume_p"], pr["prime_polarity_p"], pr["prime_charge_p"],
                pr["prime_hydrophobicity_p"], pr["prime_composition_p"],
                pr["prime_sig_properties"],
                pos_methods, neg_methods,
                pos_methods >= 2,  # consensus_pos
                neg_methods >= 2,  # consensus_neg
                h2h_stronger,
                gard_warn,
                args.virus_group, args.protein,
            ]
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"\n=== All-Site Selection Table: {args.virus_group} / {args.protein} ===")
    print(f"  Total sites     : {total}")
    print(f"  Consensus pos   : {n_pos} sites (≥2 methods agree)")
    print(f"  Contrast-FEL sig: {n_cfell_sig} sites (H2H vs Reservoir, FDR < {args.contrast_fdr})")
    print(f"  Recombination warning: {gard_warn}")
    print(f"  Output: {args.out}")


if __name__ == "__main__":
    main()
