"""
parse_branch_selection.py
=========================
Parses HyPhy BUSTED, aBSREL, and RELAX JSON outputs to extract
branch-specific selection signatures relevant to host-range change.

Key outputs:
  - BUSTED : p-value for gene-wide positive selection on Foreground branches
  - aBSREL : list of branches with significant episodic diversification (q < 0.05)
  - RELAX  : relaxation/intensification parameter k and p-value
             k > 1 => intensification, k < 1 => relaxation

Usage:
  python parse_branch_selection.py \
      --busted  results/.../GnGc_busted.json \
      --absrel   results/.../GnGc_absrel.json \
      --relax    results/.../GnGc_relax.json \
      --out      results/.../GnGc_branch_selection.tsv \
      --virus-group Andes_virus \
      --protein  GnGc
"""

import argparse
import json
import os
import sys


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_json(path):
    """Load a JSON file; return None if missing, empty, or invalid."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        with open(path) as f:
            data = json.load(f)
        return data if data else None
    except Exception as e:
        print(f"Warning: cannot read {path}: {e}", file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# BUSTED parser
# ---------------------------------------------------------------------------

def parse_busted(data):
    """
    Extract gene-level positive selection p-value from BUSTED output.

    BUSTED JSON structure (relevant keys):
      {
        "test results": {
          "p-value": <float>,
          "LRT": <float>
        }
      }

    Returns dict:
      {
        "busted_pvalue": float,   # p-value for positive selection on foreground
        "busted_lrt": float,      # likelihood ratio test statistic
        "busted_sig": bool        # True if p < 0.05
      }
    """
    result = {"busted_pvalue": "NA", "busted_lrt": "NA", "busted_sig": False}
    if not data:
        return result

    test = data.get("test results", {})
    pval = test.get("p-value", None)
    lrt  = test.get("LRT", None)

    if pval is not None:
        result["busted_pvalue"] = round(float(pval), 6)
        result["busted_lrt"]    = round(float(lrt), 4) if lrt is not None else "NA"
        result["busted_sig"]    = float(pval) < 0.05

    return result


# ---------------------------------------------------------------------------
# aBSREL parser
# ---------------------------------------------------------------------------

def parse_absrel(data, qvalue_threshold=0.05):
    """
    Extract branch-level episodic diversification from aBSREL output.

    aBSREL JSON (relevant keys):
      {
        "branch attributes": {
          "0": {                      # partition index
            "<branch_name>": {
              "Corrected P-value": <float>,
              "Uncorrected P-value": <float>,
              "Rate classes": <int>,
              ...
            }
          }
        }
      }

    Returns dict:
      {
        "absrel_sig_branches": int,          # number of branches with q < threshold
        "absrel_sig_branch_names": str,      # comma-separated names
        "absrel_foreground_sig": bool        # any Foreground branch significant?
      }
    """
    result = {
        "absrel_sig_branches": 0,
        "absrel_sig_branch_names": "None",
        "absrel_foreground_sig": False,
    }
    if not data:
        return result

    branch_attrs = data.get("branch attributes", {})
    sig_branches = []

    for partition_key, branches in branch_attrs.items():
        for branch_name, attrs in branches.items():
            qval = attrs.get("Corrected P-value", 1.0)
            if qval is None:
                qval = attrs.get("Uncorrected P-value", 1.0)
            try:
                if float(qval) < qvalue_threshold:
                    sig_branches.append(branch_name)
            except (TypeError, ValueError):
                pass

    foreground_sig = any(
        "foreground" in b.lower() or "h2h" in b.lower() or "spillover" in b.lower()
        for b in sig_branches
    )

    result["absrel_sig_branches"]      = len(sig_branches)
    result["absrel_sig_branch_names"]  = ",".join(sig_branches) if sig_branches else "None"
    result["absrel_foreground_sig"]    = foreground_sig
    return result


# ---------------------------------------------------------------------------
# RELAX parser
# ---------------------------------------------------------------------------

def parse_relax(data):
    """
    Extract relaxation/intensification parameter k and p-value from RELAX output.

    RELAX JSON (relevant keys):
      {
        "test results": {
          "relaxation or intensification parameter": <float>,  # k
          "p-value": <float>,
          "LRT": <float>
        }
      }

    Interpretation:
      k > 1  =>  intensification of selection on Foreground branches
      k < 1  =>  relaxation of selection on Foreground branches
      k = 1  =>  no change (null hypothesis)

    Returns dict:
      {
        "relax_k": float,         # relaxation parameter
        "relax_pvalue": float,
        "relax_direction": str,   # "intensification" | "relaxation" | "neutral"
        "relax_sig": bool
      }
    """
    result = {
        "relax_k": "NA",
        "relax_pvalue": "NA",
        "relax_direction": "NA",
        "relax_sig": False,
    }
    if not data:
        return result

    test = data.get("test results", {})
    k    = test.get("relaxation or intensification parameter", None)
    pval = test.get("p-value", None)

    if k is not None and pval is not None:
        k    = float(k)
        pval = float(pval)
        result["relax_k"]         = round(k, 4)
        result["relax_pvalue"]    = round(pval, 6)
        result["relax_sig"]       = pval < 0.05
        if pval < 0.05:
            result["relax_direction"] = "intensification" if k > 1 else "relaxation"
        else:
            result["relax_direction"] = "neutral (ns)"

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Parse BUSTED / aBSREL / RELAX HyPhy JSON outputs"
    )
    p.add_argument("--busted",  required=True, help="BUSTED JSON")
    p.add_argument("--absrel",  required=True, help="aBSREL JSON")
    p.add_argument("--relax",   required=True, help="RELAX JSON")
    p.add_argument("--out",     required=True, help="Output TSV summary file")
    p.add_argument("--virus-group", required=True, dest="virus_group")
    p.add_argument("--protein", required=True)
    p.add_argument("--qvalue",  type=float, default=0.05,
                   help="FDR-corrected q-value threshold for aBSREL branches (default: 0.05)")
    return p.parse_args()


def main():
    args = parse_args()

    busted_data = load_json(args.busted)
    absrel_data = load_json(args.absrel)
    relax_data  = load_json(args.relax)

    busted = parse_busted(busted_data)
    absrel = parse_absrel(absrel_data, qvalue_threshold=args.qvalue)
    relax  = parse_relax(relax_data)

    # Merge results
    row = {
        "virus_group":             args.virus_group,
        "protein":                 args.protein,
        **busted,
        **absrel,
        **relax,
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    header = [
        "virus_group", "protein",
        "busted_pvalue", "busted_lrt", "busted_sig",
        "absrel_sig_branches", "absrel_sig_branch_names", "absrel_foreground_sig",
        "relax_k", "relax_pvalue", "relax_direction", "relax_sig",
    ]
    with open(args.out, "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(str(row[h]) for h in header) + "\n")

    # Human-readable summary
    print(f"\n=== Branch-Specific Selection: {args.virus_group} / {args.protein} ===")
    print(f"  BUSTED  p-value : {row['busted_pvalue']}  (significant: {row['busted_sig']})")
    print(f"  aBSREL  sig branches : {row['absrel_sig_branches']}  ({row['absrel_sig_branch_names']})")
    print(f"  RELAX   k={row['relax_k']}  p={row['relax_pvalue']}  direction: {row['relax_direction']}")


if __name__ == "__main__":
    main()
