"""
parse_gard.py
=============
Parses HyPhy GARD (Genetic Algorithm Recombination Detection) JSON output.

GARD JSON structure (key fields):
  {
    "improvements": [              # list of recombination breakpoints found
      {
        "breakpoint": <int>,       # alignment column position of breakpoint
        "c-AIC improvement": <float>
      }, ...
    ],
    "GARD model": {
      "AIC-c": <float>
    }
  }

If "improvements" is empty => no significant recombination detected.

Output TSV columns:
  virus_group, protein, recombination_detected, n_breakpoints,
  breakpoint_positions, max_improvement, gard_aic
"""

import argparse
import json
import os
import sys


def load_json(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except Exception as e:
        print(f"Warning: cannot read {path}: {e}", file=sys.stderr)
        return None


def parse_gard(data):
    """
    Extract recombination summary from GARD JSON.

    Returns dict with:
      status                  : 'OK' | 'FAILED'
      recombination_detected  : bool, or 'NA' when status is FAILED
      n_breakpoints           : int or 'NA'
      breakpoint_positions    : str  (comma-separated column positions)
      max_c_aic_improvement   : float or 'NA'
      gard_aic                : float or 'NA'

    A failed run and a genuine "no breakpoints" result must never look alike: an
    empty or sentinel JSON yields status=FAILED with recombination_detected='NA',
    so that downstream code cannot read a crash as a clean negative.
    """
    result = {
        "status": "OK",
        "recombination_detected": False,
        "n_breakpoints": 0,
        "breakpoint_positions": "None",
        "max_c_aic_improvement": "NA",
        "gard_aic": "NA",
    }
    failed = dict(result, status="FAILED", recombination_detected="NA",
                  n_breakpoints="NA")

    # Empty dict, missing file, or the {"status": "FAILED"} sentinel written by the
    # GARD rule when HyPhy exits non-zero.
    if not data or data.get("status") == "FAILED":
        return failed

    improvements = data.get("improvements", [])
    if not improvements:
        # GARD ran but reported no improvements section at all: treat as unusable
        # rather than as evidence of no recombination.
        return failed

    # Filter to actual breakpoints (c-AIC improvement > 0)
    breakpoints = [
        imp for imp in improvements
        if isinstance(imp, dict) and imp.get("c-AIC improvement", 0) > 0
    ]

    if not breakpoints:
        return result

    result["recombination_detected"] = True
    result["n_breakpoints"] = len(breakpoints)
    result["breakpoint_positions"] = ",".join(
        str(int(bp.get("breakpoint", 0))) for bp in breakpoints
    )
    improvements_vals = [bp.get("c-AIC improvement", 0) for bp in breakpoints]
    result["max_c_aic_improvement"] = round(max(improvements_vals), 4)

    # GARD AIC from model section
    gard_model = data.get("GARD model", {})
    aic = gard_model.get("AIC-c", None)
    if aic is not None:
        result["gard_aic"] = round(float(aic), 4)

    return result


def parse_args():
    p = argparse.ArgumentParser(description="Parse HyPhy GARD JSON output")
    p.add_argument("--json",        required=True, help="GARD JSON file")
    p.add_argument("--out",         required=True, help="Output summary text file")
    p.add_argument("--virus-group", required=True, dest="virus_group")
    p.add_argument("--protein",     required=True)
    return p.parse_args()


def main():
    args = parse_args()
    data = load_json(args.json)
    result = parse_gard(data)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    header = [
        "virus_group", "protein", "status",
        "recombination_detected", "n_breakpoints",
        "breakpoint_positions", "max_c_aic_improvement", "gard_aic",
    ]
    row = {
        "virus_group": args.virus_group,
        "protein": args.protein,
        **result,
    }

    with open(args.out, "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(str(row[h]) for h in header) + "\n")

    # Human-readable console output. ASCII only: a cp1252 console raises
    # UnicodeEncodeError on emoji and would fail the rule.
    detected = result["recombination_detected"]
    print(f"\n=== GARD: {args.virus_group} / {args.protein} ===")
    if result["status"] == "FAILED":
        print("  [FAILED] GARD produced no usable output.")
        print("  Recombination status is UNKNOWN for this alignment. Do not report this")
        print("  as evidence that the no-recombination assumption holds.")
    elif detected:
        print("  [WARNING] RECOMBINATION DETECTED")
        print(f"  Breakpoints : {result['n_breakpoints']} at positions [{result['breakpoint_positions']}]")
        print(f"  Max cAIC improvement : {result['max_c_aic_improvement']}")
        print("  Downstream selection analyses may have inflated false positive rates.")
        print("  Consider using only the largest non-recombinant segment.")
    else:
        print("  [OK] No recombination breakpoints found.")
        print("  Assumption of no recombination holds for this alignment.")


if __name__ == "__main__":
    main()
