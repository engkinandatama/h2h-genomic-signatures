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
            # Jika ada daftar target_taxa spesifik, pastikan ia cocok
            if target_taxa:
                if any(t in node_name for t in target_taxa):
                    return f"{node_name}{tag}"
                else:
                    return match.group(0)
            # Jika target_taxa kosong, labeli semua daun (untuk pohon homogen H2H/Spillover)
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
        
    # Tentukan apakah kita harus melabeli pohon ini
    # Jika kategori adalah H2H atau Spillover, kita tandai seluruh daun sebagai Foreground
    # Jika kategori adalah Reservoir, kita biarkan default (tanpa label Foreground)
    is_foreground = args.category.upper() in ["H2H", "SPILLOVER"]
    
    if is_foreground or args.h2h_taxa:
        print(f"Melabeli pohon {args.tree} sebagai Foreground...")
        labeled_tree = label_leaves(tree_str, tag="{Foreground}", target_taxa=args.h2h_taxa)
    else:
        print(f"Pohon {args.tree} dikategorikan sebagai Background (tidak dilabeli).")
        labeled_tree = tree_str
        
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as f:
        f.write(labeled_tree + "\n")
        
    print(f"Tree disimpan ke {args.out}")

if __name__ == "__main__":
    main()
