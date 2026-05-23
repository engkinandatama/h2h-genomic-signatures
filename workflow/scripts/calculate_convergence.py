import argparse
import os
import sys
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate Convergent Genomic Signatures across viruses (Levels 1-4)")
    parser.add_argument("--inputs", nargs='+', required=True, help="List of HyPhy selection results txt files")
    parser.add_argument("--out_matrix", required=True, help="Output convergence matrix CSV file")
    parser.add_argument("--out_jaccard", required=True, help="Output Jaccard similarity heatmap PNG")
    parser.add_argument("--out_correlation", required=True, help="Output selection-category correlation plot PNG")
    return parser.parse_args()

def parse_virus_metadata(virus_name):
    v_lower = virus_name.lower()
    # 1. Famili Virus
    family = "Unknown"
    if "andes" in v_lower or "puumala" in v_lower or "sin_nombre" in v_lower:
        family = "Hantaviridae"
    elif "nipah" in v_lower:
        family = "Paramyxoviridae"
    elif "ebola" in v_lower:
        family = "Filoviridae"
        
    # 2. Kategori
    category = "Unknown"
    if "h2h" in v_lower:
        category = "H2H"
    elif "spillover" in v_lower:
        category = "Spillover"
    elif "reservoir" in v_lower:
        category = "Reservoir"
        
    return family, category

def get_protein_length(results_path, virus, protein):
    # Cari file FASTA protein terkait untuk mendapatkan panjang riil protein
    # Path standar: results/Zoonotic_Convergent_Signatures/01_raw_fasta/{virus}/{protein}_filtered.faa
    # Kita cari direktori secara relatif dari path file hasil seleksi
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(results_path)))
    faa_path = os.path.join(base_dir, "01_raw_fasta", virus, f"{protein}_filtered.faa")
    
    if os.path.exists(faa_path):
        try:
            with open(faa_path, 'r') as f:
                # Ambil sekuens pertama
                f.readline() # Skip header
                seq = []
                for line in f:
                    if line.startswith('>'):
                        break
                    seq.append(line.strip())
                return len("".join(seq))
        except Exception:
            pass
    # Fallback default length jika file tidak ditemukan
    return 1000

def main():
    args = parse_args()
    
    print(f"Membaca {len(args.inputs)} file hasil seleksi HyPhy...")
    
    data_list = []
    
    # 1. Parsing file input dan ekstrak metadata
    for file_path in args.inputs:
        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} tidak ditemukan. Dilewati.")
            continue
            
        # Dapatkan nama virus dan protein dari path
        # Format path: .../04_selection/{virus}/{protein}_selection_results.txt
        parts = file_path.replace("\\", "/").split("/")
        if len(parts) < 3:
            continue
            
        virus = parts[-2]
        protein_file = parts[-1]
        protein = protein_file.replace("_selection_results.txt", "")
        
        family, category = parse_virus_metadata(virus)
        
        # Baca situs terseleksi
        sites = []
        with open(file_path, 'r') as f:
            f.readline() # Skip header
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts_line = line.split('\t')
                try:
                    site = int(parts_line[0])
                    sites.append(site)
                except ValueError:
                    continue
                    
        prot_len = get_protein_length(file_path, virus, protein)
        
        data_list.append({
            'virus': virus,
            'protein': protein,
            'family': family,
            'category': category,
            'protein_length': prot_len,
            'selected_sites': sites,
            'num_selected': len(sites),
            'selection_rate': len(sites) / prot_len if prot_len > 0 else 0
        })
        
    if not data_list:
        print("Error: Tidak ada data hasil seleksi yang berhasil diproses.")
        sys.exit(1)
        
    df_results = pd.DataFrame(data_list)
    
    # Simpan matriks konvergensi ke CSV
    df_results_out = df_results.copy()
    # Ubah list ke string koma untuk penampilan CSV
    df_results_out['selected_sites'] = df_results_out['selected_sites'].apply(lambda x: ",".join(map(str, x)))
    os.makedirs(os.path.dirname(args.out_matrix), exist_ok=True)
    df_results_out.to_csv(args.out_matrix, index=False)
    print(f"Matriks konvergensi disimpan di: {args.out_matrix}")
    
    # 2. Hitung Jaccard Similarity menggunakan Vektor Binning Persen (100 bin)
    # Hal ini mengatasi perbedaan panjang protein lintas virus
    bin_matrix = []
    labels = []
    
    for idx, row in df_results.iterrows():
        bins = np.zeros(100)
        prot_len = row['protein_length']
        for site in row['selected_sites']:
            # Konversi posisi residu ke persentase 0-99
            if prot_len > 0:
                bin_idx = int((site - 1) / prot_len * 100)
                bin_idx = min(99, max(0, bin_idx))
                bins[bin_idx] = 1
        bin_matrix.append(bins)
        labels.append(f"{row['virus']} - {row['protein']}")
        
    bin_matrix = np.array(bin_matrix)
    n_samples = len(labels)
    jaccard_mat = np.zeros((n_samples, n_samples))
    
    for i in range(n_samples):
        for j in range(n_samples):
            intersection = np.sum(np.logical_and(bin_matrix[i], bin_matrix[j]))
            union = np.sum(np.logical_or(bin_matrix[i], bin_matrix[j]))
            jaccard_mat[i, j] = intersection / union if union > 0 else 0.0
            
    df_jaccard = pd.DataFrame(jaccard_mat, index=labels, columns=labels)
    
    # Plot Jaccard Heatmap
    plt.figure(figsize=(12, 10))
    sns.set_theme(style="darkgrid")
    sns.heatmap(
        df_jaccard, 
        annot=True, 
        fmt=".2f", 
        cmap="coolwarm", 
        cbar_kws={'label': 'Jaccard Similarity Index'},
        linewidths=0.5
    )
    plt.title("Jaccard Similarity Heatmap of Positive Selection Sites (Normalized Bins)", fontsize=14, pad=15)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(args.out_jaccard), exist_ok=True)
    plt.savefig(args.out_jaccard, dpi=300)
    plt.close()
    print(f"Heatmap Jaccard disimpan di: {args.out_jaccard}")
    
    # 3. Uji Korelasi & Boxplot Perbandingan Kategori (H2H vs Spillover vs Reservoir)
    plt.figure(figsize=(8, 6))
    
    # Map kategori ke warna
    palette = {'H2H': '#ef4444', 'Spillover': '#f59e0b', 'Reservoir': '#3b82f6'}
    
    sns.boxplot(
        data=df_results, 
        x='category', 
        y='selection_rate', 
        palette=palette,
        hue='category',
        legend=False
    )
    sns.stripplot(
        data=df_results, 
        x='category', 
        y='selection_rate', 
        color='black', 
        alpha=0.5, 
        size=6
    )
    
    # Uji statistik Mann-Whitney U (H2H vs Reservoir)
    h2h_rates = df_results[df_results['category'] == 'H2H']['selection_rate'].values
    res_rates = df_results[df_results['category'] == 'Reservoir']['selection_rate'].values
    
    p_val_text = ""
    if len(h2h_rates) > 0 and len(res_rates) > 0:
        try:
            stat, p_val = mannwhitneyu(h2h_rates, res_rates, alternative='two-sided')
            p_val_text = f"Mann-Whitney U Test (H2H vs Reservoir):\np-value = {p_val:.4f}"
            if p_val < 0.05:
                p_val_text += " (Significant)"
            else:
                p_val_text += " (Not Significant)"
        except Exception as e:
            p_val_text = f"Statistical test error: {e}"
            
    # Tampilkan teks uji statistik di plot
    if p_val_text:
        plt.text(
            0.05, 0.95, 
            p_val_text, 
            transform=plt.gca().transAxes, 
            fontsize=10, 
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
        )
        
    plt.title("Comparison of Positive Selection Rates across Transmission Categories", fontsize=12, pad=15)
    plt.xlabel("Transmission Category", fontsize=10)
    plt.ylabel("Positive Selection Rate (Selected Sites / Protein Length)", fontsize=10)
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(args.out_correlation), exist_ok=True)
    plt.savefig(args.out_correlation, dpi=300)
    plt.close()
    print(f"Plot korelasi kategori disimpan di: {args.out_correlation}")

if __name__ == "__main__":
    main()
