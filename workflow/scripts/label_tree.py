import argparse
import sys
import os
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Label Phylogenetic Tree Branches for HyPhy MEME/FEL/aBSREL"
    )
    parser.add_argument("--tree", required=True, help="Input Newick tree from IQ-TREE")
    parser.add_argument("--out",  required=True, help="Output labeled Newick tree")
    parser.add_argument("--category",   default="",
                        help="Virus category (H2H, Spillover, Reservoir)")
    parser.add_argument("--h2h_taxa",   nargs="*", default=[],
                        help="List of specific H2H taxa names to label")
    parser.add_argument("--bootstrap-threshold", type=float, default=0.0,
                        dest="bootstrap_threshold",
                        help="Collapse internal branches with bootstrap support "
                             "below this value (0-100). Default: 0 (no filtering). "
                             "Recommended: 70 for standard analyses.")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Bootstrap filtering
# ---------------------------------------------------------------------------

def collapse_low_bootstrap_nodes(newick_str, threshold):
    """
    Remove internal node labels (bootstrap values) below `threshold`.

    IQ-TREE Newick format example:
      (A:0.1, B:0.2)85:0.3   <- internal node bootstrap=85
      (A:0.1, B:0.2):0.3     <- internal node without bootstrap

    Strategy: find numeric internal node labels between ')' and ':',
    and remove them if below threshold. This prevents HyPhy from being
    misled by poorly-supported splits in branch-specific selection tests.

    Note: This is label-removal only (not true polytomy collapsing),
    which is sufficient for HyPhy's likelihood-based framework.
    """
    if threshold <= 0:
        return newick_str

    def replace_internal_label(match):
        label = match.group(1)
        try:
            val = float(label)
            if 0 <= val <= 100:
                return ")" if val < threshold else f"){label}"
        except ValueError:
            pass
        return f"){label}"

    # Match: content between ')' and ':' that looks like a bootstrap number
    filtered = re.sub(r"\)([A-Za-z0-9_.]*?)(?=:)", replace_internal_label, newick_str)
    return filtered


# ---------------------------------------------------------------------------
# Branch labeling
# ---------------------------------------------------------------------------

def label_leaves(newick_str, tag="{Foreground}", target_taxa=None):
    """
    Label leaf nodes (sequence IDs) that represent H2H/Spillover taxa.
    Applies HyPhy branch annotation format: TaxonName{Foreground}.
    """
    def replace_leaf(match):
        node_name = match.group(1)
        # Skip pure numeric tokens (bootstrap values or branch lengths)
        if re.match(r"^[0-9.]+$", node_name):
            return match.group(0)

        if any(c.isalpha() for c in node_name):
            if target_taxa:
                if any(t.lower() in node_name.lower() for t in target_taxa):
                    return f"{node_name}{tag}"
                return match.group(0)

            # Default: label if sequence ID contains h2h or spillover
            name_lower = node_name.lower()
            if "h2h" in name_lower or "spillover" in name_lower:
                return f"{node_name}{tag}"

        return match.group(0)

    # Match node names immediately before ':branch_length'
    labeled_str = re.sub(r"([a-zA-Z0-9_|.-]+)(?=:)", replace_leaf, newick_str)
    return labeled_str


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not os.path.exists(args.tree):
        print(f"Error: Tree file {args.tree} tidak ditemukan.")
        sys.exit(1)

    with open(args.tree, "r") as f:
        tree_str = f.read().strip()

    # Step 1: Remove low-bootstrap internal nodes (if threshold set)
    if args.bootstrap_threshold > 0:
        print(f"Filtering internal nodes with bootstrap < {args.bootstrap_threshold}...")
        tree_str = collapse_low_bootstrap_nodes(tree_str, args.bootstrap_threshold)
        print("Bootstrap filtering complete.")

    # Step 2: Label H2H/Spillover leaves as Foreground
    print(f"Labeling tree {args.tree} for Foreground (H2H/Spillover) branches...")
    labeled_tree = label_leaves(tree_str, tag="{Foreground}", target_taxa=args.h2h_taxa)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, "w") as f:
        f.write(labeled_tree + "\n")

    print(f"Labeled tree saved to {args.out}")


if __name__ == "__main__":
    main()
