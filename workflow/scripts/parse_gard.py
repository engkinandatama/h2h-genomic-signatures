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


def _improvement_stages(improvements):
    """
    GARD writes 'improvements' as a dict keyed by stage index, not a list.

    The keys are the strings "0", "1", ... and each value holds that stage's
    incremental deltaAICc and the breakpoint set found at it. Iterating the dict
    directly yields those key strings, so a loop written for a list of records
    never matches anything and silently reports no recombination -- which is what
    this parser did for every dataset, including alignments whose own GARD log
    named the breakpoints it had found.
    """
    if isinstance(improvements, dict):
        return [improvements[k] for k in sorted(improvements, key=lambda x: int(x))]
    if isinstance(improvements, list):
        return list(improvements)
    return []


def _flatten_breakpoints(bps):
    """HyPhy matrices serialise as nested lists: [[5138]] or [[5097], [5636]]."""
    out = []
    if not bps:
        return out
    for item in bps if isinstance(bps, list) else [bps]:
        if isinstance(item, list):
            out.extend(int(float(x)) for x in item)
        else:
            out.append(int(float(item)))
    return out


def parse_gard(data):
    """
    Extract recombination summary from GARD JSON.

    Returns dict with:
      status                  : 'OK' | 'FAILED' | 'SKIPPED'
      recombination_detected  : bool, or 'NA' when not OK
      n_breakpoints           : int or 'NA'
      breakpoint_positions    : str  (comma-separated alignment columns)
      delta_aicc_vs_baseline  : float or 'NA'
      baseline_aicc           : float or 'NA'
      best_model_aicc         : float or 'NA'

    The breakpoint count comes from 'breakpointData', which holds one entry per
    partition of the best model and is therefore one longer than the number of
    breakpoints. That is the model GARD actually settled on, so it is preferred
    over re-deriving the count from the improvement stages.

    A failed run and a genuine "no breakpoints" result must never look alike: an
    empty or sentinel JSON yields status=FAILED with recombination_detected='NA',
    so that downstream code cannot read a crash as a clean negative.
    """
    result = {
        "status": "OK",
        "recombination_detected": False,
        "n_breakpoints": 0,
        "breakpoint_positions": "None",
        "delta_aicc_vs_baseline": "NA",
        "baseline_aicc": "NA",
        "best_model_aicc": "NA",
    }
    failed = dict(result, status="FAILED", recombination_detected="NA",
                  n_breakpoints="NA")

    if data and data.get("status") == "SKIPPED":
        # Not screened because compute was the constraint. Distinct from FAILED so
        # the Methods can say how many alignments were screened and how many were
        # not, and distinct from a clean negative in either case.
        return dict(failed, status="SKIPPED")
    if not data or data.get("status") == "FAILED":
        return failed

    if "improvements" not in data:
        # GARD ran but wrote no improvements section at all: treat as unusable
        # rather than as evidence of no recombination.
        return failed

    for key, field in (("baselineScore", "baseline_aicc"),
                       ("bestModelAICc", "best_model_aicc")):
        if data.get(key) is not None:
            result[field] = round(float(data[key]), 4)
    if result["baseline_aicc"] != "NA" and result["best_model_aicc"] != "NA":
        result["delta_aicc_vs_baseline"] = round(
            result["baseline_aicc"] - result["best_model_aicc"], 4)

    partitions = data.get("breakpointData") or {}
    stages = _improvement_stages(data.get("improvements"))
    positions = []
    for stage in reversed(stages):
        positions = _flatten_breakpoints(
            stage.get("breakpoints") if isinstance(stage, dict) else None)
        if positions:
            break

    n_from_partitions = max(len(partitions) - 1, 0) if partitions else None
    n_breakpoints = n_from_partitions if n_from_partitions is not None else len(positions)

    if n_breakpoints == 0:
        return result

    result["recombination_detected"] = True
    result["n_breakpoints"] = n_breakpoints
    if positions:
        result["breakpoint_positions"] = ",".join(str(p) for p in sorted(positions))
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
        "breakpoint_positions", "delta_aicc_vs_baseline",
        "baseline_aicc", "best_model_aicc",
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
    if result["status"] == "SKIPPED":
        print("  [SKIPPED] Not screened (params.gard.skip_datasets).")
        print("  Recombination status is UNKNOWN for this alignment.")
    elif result["status"] == "FAILED":
        print("  [FAILED] GARD produced no usable output.")
        print("  Recombination status is UNKNOWN for this alignment. Do not report this")
        print("  as evidence that the no-recombination assumption holds.")
    elif detected:
        print("  [WARNING] RECOMBINATION DETECTED")
        print(f"  Breakpoints : {result['n_breakpoints']} at positions [{result['breakpoint_positions']}]")
        print(f"  cAIC improvement over the single-partition model : {result['delta_aicc_vs_baseline']}")
        print("  Downstream selection analyses may have inflated false positive rates.")
        print("  Consider using only the largest non-recombinant segment.")
    else:
        print("  [OK] No recombination breakpoints found.")
        print("  Assumption of no recombination holds for this alignment.")


if __name__ == "__main__":
    main()
