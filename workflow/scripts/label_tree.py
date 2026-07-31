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
    parser.add_argument("--foreground-prefixes", nargs="*", default=[],
                        dest="foreground_prefixes",
                        help="Virus names whose records form the foreground. The "
                             "merge step prefixes every header with its virus "
                             "name, so these are matched as exact prefixes.")
    parser.add_argument("--reference-prefixes", nargs="*", default=[],
                        dest="reference_prefixes",
                        help="Virus names whose records form the reference set.")
    parser.add_argument("--require-contrast", action="store_true",
                        dest="require_contrast",
                        help="Fail if either branch set ends up empty. Set for "
                             "groups that feed Contrast-FEL and RELAX.")
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

def label_leaves(newick_str, tag="{Foreground}", target_taxa=None,
                 fg_prefixes=None, ref_prefixes=None, counts=None):
    """
    Label leaf nodes (sequence IDs) that represent H2H/Spillover taxa.
    Applies HyPhy branch annotation format: TaxonName{Foreground}.

    When fg_prefixes/ref_prefixes are given, tips are partitioned by the virus
    each record came from, which merge_sequences.py has already written into the
    header as a prefix. That is authoritative. The substring heuristic below is
    only a fallback: it guessed the partition from the words "h2h", "spillover"
    and "reservoir" appearing anywhere in a name, and when a group's human-derived
    dataset was empty it produced a tree with a Reference set and no Foreground
    set at all. HyPhy then failed with "'Foreground' is not a valid choice",
    which says nothing about the missing sequences that actually caused it.
    """
    # Longest prefix first, so one virus name that begins with another cannot
    # capture the other's tips.
    ordered = sorted(
        [(p, tag) for p in (fg_prefixes or [])]
        + [(p, "{Reference}") for p in (ref_prefixes or [])],
        key=lambda kv: -len(kv[0]))

    def replace_leaf(match):
        node_name = match.group(1)
        # Skip pure numeric tokens (bootstrap values or branch lengths)
        if re.match(r"^[0-9.]+$", node_name):
            return match.group(0)

        if any(c.isalpha() for c in node_name):
            if ordered:
                for prefix, label in ordered:
                    if node_name.startswith(prefix):
                        if counts is not None:
                            counts[label] = counts.get(label, 0) + 1
                        return f"{node_name}{label}"
                if counts is not None:
                    counts["unmatched"] = counts.get("unmatched", 0) + 1
                return match.group(0)

            if target_taxa:
                if any(t.lower() in node_name.lower() for t in target_taxa):
                    return f"{node_name}{tag}"
                return match.group(0)

            # Default: tips are partitioned by the host the isolate came from.
            # Human- and outbreak-derived isolates form the foreground; reservoir
            # isolates form the reference. Both sets are labelled explicitly:
            # RELAX is invoked with --reference Reference, and relying on HyPhy's
            # undocumented "everything unlabelled" fallback made the reference set
            # implicit and version-dependent.
            #
            # NOTE ON SCOPE: only tips carry a host, so only tips are labelled.
            # The contrast therefore measures selection associated with the host
            # an isolate was sampled from, not with a lineage's transmission mode.
            name_lower = node_name.lower()
            if "h2h" in name_lower or "spillover" in name_lower:
                return f"{node_name}{tag}"
            if "reservoir" in name_lower:
                return f"{node_name}{{Reference}}"

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
    counts = {}
    labeled_tree = label_leaves(tree_str, tag="{Foreground}",
                                target_taxa=args.h2h_taxa,
                                fg_prefixes=args.foreground_prefixes,
                                ref_prefixes=args.reference_prefixes,
                                counts=counts)

    n_fg = counts.get("{Foreground}", 0)
    n_ref = counts.get("{Reference}", 0)
    n_un = counts.get("unmatched", 0)
    print(f"Tips labelled: {n_fg} Foreground, {n_ref} Reference, {n_un} unmatched.")
    if args.foreground_prefixes:
        print(f"  foreground prefixes: {', '.join(args.foreground_prefixes)}")
    if args.reference_prefixes:
        print(f"  reference prefixes : {', '.join(args.reference_prefixes)}")

    if args.require_contrast and (n_fg == 0 or n_ref == 0):
        # Stop here rather than write a tree that cannot support the contrast.
        # Contrast-FEL and RELAX would otherwise fail deep inside HyPhy with an
        # error that names the missing label but not the missing data.
        print(f"ERROR: this group is configured for a branch contrast, but the "
              f"tree has {n_fg} foreground and {n_ref} reference tips. One side "
              f"is empty, so no contrast exists.", file=sys.stderr)
        print(f"       Check that both viruses in the group produced sequences; "
              f"an empty fetch upstream is the usual cause.", file=sys.stderr)
        sys.exit(1)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, "w") as f:
        f.write(labeled_tree + "\n")

    print(f"Labeled tree saved to {args.out}")


if __name__ == "__main__":
    main()
