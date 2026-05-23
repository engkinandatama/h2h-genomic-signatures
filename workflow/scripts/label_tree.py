import argparse
import sys
import os
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Label Phylogenetic Tree Branches for HyPhy MEME/FEL")
    parser.add_argument("--tree", required=True, help="Input Newick tree from IQ-TREE")
    parser.add_argument("--out", required=True, help="Output labeled Newick tree")
    parser.add_argument("--category", default="", help="Virus category (H2H, Spillover, Reservoir)")
    parser.add_argument("--h2h_taxa", nargs='*', default=[], help="List of specific H2H taxa names to label")
    return parser.parse_args()

def label_leaves(newick_str, tag="{Foreground}", target_taxa=None):
    def replace_leaf(match):
        node_name = match.group(1)
        # Jika node hanya angka/titik (bootstrap value atau branch length), jangan dilabeli
        if re.match(r'^[0-9.]+$', node_name):
            return match.group(0)
        
        # Jika nama node mengandung huruf (menandakan ID sekuens daun)
        if any(c.isalpha() for c in node_name):
            # Jika ada daftar target_taxa spesifik, pastikan cocok
            if target_taxa:
                if any(t.lower() in node_name.lower() for t in target_taxa):
                    return f"{node_name}{tag}"
                else:
                    return match.group(0)
            
            # Default Opsi A: Labeli jika mengandung h2h atau spillover
            name_lower = node_name.lower()
            if "h2h" in name_lower or "spillover" in name_lower:
                return f"{node_name}{tag}"
            
        return match.group(0)
        
    # Cari nama node sebelum titik dua ':' yang menandai panjang cabang
    labeled_str = re.sub(r'([a-zA-Z0-9_|.-]+)(?=:)', replace_leaf, newick_str)
    return labeled_str

def main():
    args = parse_args()
    
    if not os.path.exists(args.tree):
        print(f"Error: Tree file {args.tree} tidak ditemukan.")
        sys.exit(1)
        
    with open(args.tree, 'r') as f:
        tree_str = f.read().strip()
        
    # Selalu jalankan labeling untuk pohon gabungan (Opsi A)
    print(f"Melabeli pohon {args.tree} untuk cabang Foreground (H2H/Spillover)...")
    labeled_tree = label_leaves(tree_str, tag="{Foreground}", target_taxa=args.h2h_taxa)
        
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, 'w') as f:
        f.write(labeled_tree + "\n")
        
    print(f"Tree disimpan ke {args.out}")

if __name__ == "__main__":
    main()

