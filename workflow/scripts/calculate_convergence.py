import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests

# Domain boundaries for Level 2 Domain Enrichment Analysis
DOMAINS = {
    "gngc": {
        "Gn_domain": (1, 400),
        "Gc_domain": (401, 650),
        "Tail_domain": (651, 1000)
    },
    "l_protein": {
        "N_terminal": (1, 599),
        "RdRp_catalytic": (600, 1200),
        "C_terminal": (1201, 2500)
    },
    "g_protein": {
        "Stalk_domain": (1, 299),
        "RBD_domain": (300, 600)
    },
    "f_protein": {
        "Fusion_peptide": (100, 150),
        "HR1_domain": (151, 250),
        "HR2_domain": (450, 500)
    },
    "gp": {
        "RBD_domain": (50, 150),
        "Mucin_like_domain": (300, 500)
    }
}

# Entry/Replication role mapping for homologous comparison (Level 3)
ROLE_MAP = {
    "GnGc": "Entry",
    "G_protein": "Entry",
    "F_protein": "Entry_Helper",
    "GP": "Entry",
    "L_protein": "Replication"
}

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate Convergent Genomic Signatures across viruses (Levels 1-4)")
    parser.add_argument("--inputs", nargs='+', required=True, help="List of HyPhy selection results txt files")
    parser.add_argument("--out_matrix", required=True, help="Output convergence matrix CSV file")
    parser.add_argument("--out_jaccard", required=True, help="Output Jaccard similarity heatmap PNG")
    parser.add_argument("--out_correlation", required=True, help="Output selection-category correlation plot PNG")
    return parser.parse_args()

def parse_group_metadata(group_name):
    metadata = {
        "Andes_virus": ("Hantaviridae", "H2H"),
        "Sin_Nombre_virus": ("Hantaviridae", "Spillover"),
        "Puumala_virus": ("Hantaviridae", "Spillover"),
        "Nipah_Bangladesh": ("Paramyxoviridae", "H2H"),
        "Nipah_Malaysia": ("Paramyxoviridae", "Spillover"),
        "Ebola_Zaire": ("Filoviridae", "H2H")
    }
    return metadata.get(group_name, ("Unknown", "Unknown"))

def get_protein_length(results_path, virus_group, protein):
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(results_path)))
    # Try reading the merged protein FASTA file first
    faa_path = os.path.join(base_dir, "02_aligned", virus_group, f"{protein}_merged.faa")
    
    if os.path.exists(faa_path):
        try:
            with open(faa_path, 'r') as f:
                f.readline() # Skip header
                seq = []
                for line in f:
                    if line.startswith('>'):
                        break
                    seq.append(line.strip())
                return len("".join(seq))
        except Exception:
            pass
    return 1000 # Fallback

def main():
    args = parse_args()
    
    print(f"Membaca {len(args.inputs)} file hasil seleksi...")
    
    data_list = []
    
    # 1. Parsing file input dan ekstrak metadata
    for file_path in args.inputs:
        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} tidak ditemukan. Dilewati.")
            continue
            
        parts = file_path.replace("\\", "/").split("/")
        if len(parts) < 3:
            continue
            
        virus_group = parts[-2]
        protein_file = parts[-1]
        protein = protein_file.replace("_selection_results.txt", "")
        
        family, category = parse_group_metadata(virus_group)
        
        # Check if the group was skipped in the pipeline (empty or missing codon alignment)
        codon_aln_path = file_path.replace("04_selection", "02_aligned").replace("_selection_results.txt", "_codon_aligned.fasta")
        if not os.path.exists(codon_aln_path) or os.path.getsize(codon_aln_path) == 0:
            print(f"Information: {virus_group} - {protein} dilewati (insufficient sequences). Mengabaikan dari analisis statistik.", flush=True)
            continue
        
        # Baca situs terseleksi (baris di file hanya situs yang signifikan)
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
                    
        prot_len = get_protein_length(file_path, virus_group, protein)
        
        data_list.append({
            'virus_group': virus_group,
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
    df_results_out['selected_sites'] = df_results_out['selected_sites'].apply(lambda x: ",".join(map(str, x)))
    os.makedirs(os.path.dirname(args.out_matrix) or ".", exist_ok=True)
    df_results_out.to_csv(args.out_matrix, index=False)
    print(f"Matriks konvergensi disimpan di: {args.out_matrix}")
    
    # =========================================================================
    # LEVEL 1: Perbandingan Jumlah Positively Selected Sites (H2H vs Spillover)
    # =========================================================================
    print("\n--- LEVEL 1: Mann-Whitney U test (H2H vs Spillover) ---")
    h2h_data = df_results[df_results['category'] == 'H2H']
    spillover_data = df_results[df_results['category'] == 'Spillover']
    
    h2h_rates = h2h_data['selection_rate'].values
    spillover_rates = spillover_data['selection_rate'].values
    
    level1_p = 1.0
    if len(h2h_rates) > 0 and len(spillover_rates) > 0:
        try:
            stat, level1_p = mannwhitneyu(h2h_rates, spillover_rates, alternative='two-sided')
            # Uji per protein (Bonferroni corrected)
            print(f"Pooled Mann-Whitney U test: U={stat}, p-value={level1_p:.6f}")
        except Exception as e:
            print(f"Error Level 1: {e}")
            
    # =========================================================================
    # LEVEL 2: Domain Enrichment Analysis (Fisher's Exact Test)
    # =========================================================================
    print("\n--- LEVEL 2: Domain Enrichment Analysis ---")
    domain_tests = []
    
    # Kelompokkan data berdasarkan protein
    for prot_name, prot_df in df_results.groupby('protein'):
        prot_lower = prot_name.lower()
        if prot_lower not in DOMAINS:
            continue
            
        prot_domains = DOMAINS[prot_lower]
        for dom_name, (start, end) in prot_domains.items():
            # Hitung jumlah selected dan non-selected sites untuk H2H vs Spillover
            h2h_sel = 0
            h2h_non = 0
            spill_sel = 0
            spill_non = 0
            
            # Filter baris untuk H2H
            h2h_rows = prot_df[prot_df['category'] == 'H2H']
            for _, row in h2h_rows.iterrows():
                length = row['protein_length']
                # Batasi domain dengan panjang protein riil
                dom_start = max(1, start)
                dom_end = min(length, end)
                dom_len = max(0, dom_end - dom_start + 1)
                
                sel_in_dom = sum(1 for s in row['selected_sites'] if dom_start <= s <= dom_end)
                non_in_dom = max(0, dom_len - sel_in_dom)
                
                h2h_sel += sel_in_dom
                h2h_non += non_in_dom
                
            # Filter baris untuk Spillover
            spill_rows = prot_df[prot_df['category'] == 'Spillover']
            for _, row in spill_rows.iterrows():
                length = row['protein_length']
                dom_start = max(1, start)
                dom_end = min(length, end)
                dom_len = max(0, dom_end - dom_start + 1)
                
                sel_in_dom = sum(1 for s in row['selected_sites'] if dom_start <= s <= dom_end)
                non_in_dom = max(0, dom_len - sel_in_dom)
                
                spill_sel += sel_in_dom
                spill_non += non_in_dom
                
            # Fisher's Exact Test
            table = [[h2h_sel, h2h_non], [spill_sel, spill_non]]
            odds_ratio = 1.0
            p_val = 1.0
            if (h2h_sel + h2h_non) > 0 and (spill_sel + spill_non) > 0:
                try:
                    odds_ratio, p_val = fisher_exact(table, alternative='two-sided')
                except Exception:
                    pass
                    
            domain_tests.append({
                'protein': prot_name,
                'domain': dom_name,
                'h2h_selected': h2h_sel,
                'h2h_non_selected': h2h_non,
                'spillover_selected': spill_sel,
                'spillover_non_selected': spill_non,
                'odds_ratio': odds_ratio,
                'p_value': p_val
            })
            
    if domain_tests:
        df_domains = pd.DataFrame(domain_tests)
        # Lakukan koreksi FDR Benjamini-Hochberg
        pvals = df_domains['p_value'].values
        reject, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        df_domains['FDR_corrected_p'] = pvals_corrected
        df_domains['significant'] = reject
        
        print(df_domains[['protein', 'domain', 'h2h_selected', 'spillover_selected', 'odds_ratio', 'p_value', 'FDR_corrected_p', 'significant']].to_string(index=False))
    else:
        print("Tidak ada domain fungsional yang terdefinisi untuk protein dalam dataset.")
        
    # =========================================================================
    # LEVEL 3: Jaccard Similarity + Permutation Test (n=10,000)
    # =========================================================================
    print("\n--- LEVEL 3: Jaccard Similarity & Permutation Test (10K Permutations) ---")
    
    # 6 virus groups
    groups = ["Andes_virus", "Sin_Nombre_virus", "Puumala_virus", "Nipah_Bangladesh", "Nipah_Malaysia", "Ebola_Zaire"]
    n_groups = len(groups)
    
    # Kita bangun profil binning 200 bin untuk setiap group
    # 100 bin untuk Entry protein, 100 bin untuk Replication protein
    group_vectors = {}
    for g in groups:
        vector = np.zeros(200)
        g_df = df_results[df_results['virus_group'] == g]
        
        # Entry protein
        entry_row = g_df[g_df['protein'].apply(lambda p: ROLE_MAP.get(p) == "Entry")]
        if not entry_row.empty:
            row = entry_row.iloc[0]
            prot_len = row['protein_length']
            for site in row['selected_sites']:
                if prot_len > 0:
                    bin_idx = min(99, max(0, int((site - 1) / prot_len * 100)))
                    vector[bin_idx] = 1
                    
        # Replication protein (L_protein atau L)
        rep_row = g_df[g_df['protein'].apply(lambda p: ROLE_MAP.get(p) == "Replication")]
        if not rep_row.empty:
            row = rep_row.iloc[0]
            prot_len = row['protein_length']
            for site in row['selected_sites']:
                if prot_len > 0:
                    bin_idx = min(99, max(0, int((site - 1) / prot_len * 100)))
                    vector[100 + bin_idx] = 1
                    
        group_vectors[g] = vector
        
    # Hitung matriks Jaccard dan p-value empiris via Permutation Test (10K)
    jaccard_mat = np.zeros((n_groups, n_groups))
    pval_mat = np.zeros((n_groups, n_groups))
    
    np.random.seed(42)
    n_permutations = 10000
    
    for i in range(n_groups):
        for j in range(n_groups):
            g_i = groups[i]
            g_j = groups[j]
            v_i = group_vectors[g_i]
            v_j = group_vectors[g_j]
            
            # Hitung Jaccard observed
            intersection = np.sum(np.logical_and(v_i, v_j))
            union = np.sum(np.logical_or(v_i, v_j))
            j_obs = intersection / union if union > 0 else 0.0
            jaccard_mat[i, j] = j_obs
            
            if i == j:
                pval_mat[i, j] = 0.0
                continue
                
            # Hitung jumlah active bins (k_i dan k_j)
            k_i = int(np.sum(v_i))
            k_j = int(np.sum(v_j))
            
            # Permutation test
            # Kita acak k_i dan k_j bins dalam vektor 200 elemen sebanyak 10.000 kali
            if k_i > 0 and k_j > 0:
                count_geq = 0
                for _ in range(n_permutations):
                    # Acak posisi bin
                    rand_i = np.zeros(200)
                    rand_i[np.random.choice(200, k_i, replace=False)] = 1
                    
                    rand_j = np.zeros(200)
                    rand_j[np.random.choice(200, k_j, replace=False)] = 1
                    
                    # Hitung Jaccard acak
                    inter_rand = np.sum(np.logical_and(rand_i, rand_j))
                    union_rand = np.sum(np.logical_or(rand_i, rand_j))
                    j_rand = inter_rand / union_rand if union_rand > 0 else 0.0
                    
                    if j_rand >= j_obs:
                        count_geq += 1
                        
                pval_mat[i, j] = count_geq / n_permutations
            else:
                pval_mat[i, j] = 1.0
                
    df_jaccard = pd.DataFrame(jaccard_mat, index=groups, columns=groups)
    df_pvals = pd.DataFrame(pval_mat, index=groups, columns=groups)
    
    print("\nMatriks Jaccard Similarity:")
    print(df_jaccard)
    print("\nMatriks P-value Empiris:")
    print(df_pvals)
    
    # Plot Jaccard Heatmap
    plt.figure(figsize=(10, 8))
    # Buat label dengan asterisk untuk significance (p < 0.05)
    annot_labels = np.empty((n_groups, n_groups), dtype=object)
    for i in range(n_groups):
        for j in range(n_groups):
            val = jaccard_mat[i, j]
            p = pval_mat[i, j]
            if i != j and p < 0.05:
                annot_labels[i, j] = f"{val:.2f}*"
            else:
                annot_labels[i, j] = f"{val:.2f}"
                
    sns.heatmap(
        df_jaccard, 
        annot=annot_labels, 
        fmt="", 
        cmap="Oranges", 
        cbar_kws={'label': 'Jaccard Similarity Index'},
        linewidths=0.5
    )
    plt.title("Jaccard Similarity matrix (Entry + Replication Bins)\n* indicates significant overlap (permutation p < 0.05)", fontsize=12, pad=15)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(args.out_jaccard) or ".", exist_ok=True)
    plt.savefig(args.out_jaccard, dpi=300)
    plt.close()
    print(f"Heatmap Jaccard disimpan di: {args.out_jaccard}")
    
    # =========================================================================
    # LEVEL 4: Korelasi Intensitas Seleksi dengan Transmission Phenotype
    # =========================================================================
    print("\n--- LEVEL 4: Spearman Correlation & Leave-One-Out Sensitivity ---")
    
    # Hitung total selection rate per group (jumlah selected sites di Entry + Replication dibagi total length)
    group_rates = []
    for g in groups:
        g_df = df_results[df_results['virus_group'] == g]
        if g_df.empty:
            print(f"Warning: Tidak ada data hasil seleksi untuk {g}. Kelompok ini akan diabaikan dari korelasi Spearman.", flush=True)
            continue
            
        _, category = parse_group_metadata(g)
        score = 2 if category == "H2H" else 1 # H2H=2, Spillover=1
        
        total_selected = g_df['num_selected'].sum()
        total_length = g_df['protein_length'].sum()
        rate = total_selected / total_length if total_length > 0 else 0.0
        
        group_rates.append({
            'virus_group': g,
            'category': category,
            'score': score,
            'selection_rate': rate
        })
        
    if not group_rates:
        print("\nWarning: Tidak ada data grup yang valid untuk analisis korelasi Level 4. Menulis plot kosong.")
        # Create empty plot or skip saving
        plt.figure(figsize=(6, 5))
        plt.title("No Data Available for Positive Selection Rates vs Transmission Phenotype")
        os.makedirs(os.path.dirname(args.out_correlation) or ".", exist_ok=True)
        plt.savefig(args.out_correlation, dpi=300)
        plt.close()
    else:
        df_loo = pd.DataFrame(group_rates)
        print("Data untuk korelasi Spearman:")
        print(df_loo)
        
        # Hitung Spearman correlation observed
        obs_rho, obs_p = np.nan, np.nan
        if len(df_loo) >= 2:
            try:
                obs_rho, obs_p = spearmanr(df_loo['selection_rate'], df_loo['score'])
                print(f"\nObserved Spearman Correlation: rho={obs_rho:.4f}, p-value={obs_p:.4f}")
            except Exception as e:
                print(f"\nWarning: Gagal menghitung korelasi Spearman Observed: {e}")
        else:
            print("\nWarning: Data tidak mencukupi untuk menghitung korelasi Spearman Observed (butuh minimal 2 kelompok).")
            
        # Leave-One-Out Sensitivity Analysis
        print("\nLeave-One-Out Sensitivity Analysis:")
        loo_results = []
        for excluded in groups:
            # check if excluded is in df_loo
            if excluded not in df_loo['virus_group'].values:
                continue
            sub_df = df_loo[df_loo['virus_group'] != excluded]
            rho, p = np.nan, np.nan
            if len(sub_df) >= 2:
                try:
                    rho, p = spearmanr(sub_df['selection_rate'], sub_df['score'])
                except Exception as e:
                    pass
            loo_results.append({
                'excluded_virus': excluded,
                'spearman_rho': rho,
                'p_value': p
            })
        if loo_results:
            df_loo_res = pd.DataFrame(loo_results)
            print(df_loo_res.to_string(index=False))
        else:
            print("Tidak ada hasil Leave-One-Out karena data kurang.")
            
        # Plot Boxplot/Strip plot perbandingan rate H2H vs Spillover
        plt.figure(figsize=(6, 5))
        palette = {'H2H': '#ef4444', 'Spillover': '#f59e0b'}
        
        sns.boxplot(
            data=df_loo, 
            x='category', 
            y='selection_rate', 
            palette=palette,
            hue='category',
            legend=False,
            width=0.4
        )
        sns.stripplot(
            data=df_loo, 
            x='category', 
            y='selection_rate', 
            color='black', 
            alpha=0.8, 
            size=8
        )
        
        # Tampilkan teks Spearman di plot
        stats_text = f"Spearman Correlation:\nrho = {obs_rho:.3f}\np-value = {obs_p:.3f}"
        plt.text(
            0.05, 0.95, 
            stats_text, 
            transform=plt.gca().transAxes, 
            fontsize=10, 
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
        )
        
        plt.title("Positive Selection Rates vs Transmission Phenotype", fontsize=11, pad=15)
        plt.xlabel("Transmission Category (H2H score=2, Spillover score=1)", fontsize=10)
        plt.ylabel("Pooled Selection Rate (Entry + Replication)", fontsize=10)
        max_rate = df_loo['selection_rate'].max()
        plt.ylim(0, (max_rate if pd.notna(max_rate) and max_rate > 0 else 0.01) * 1.3)
        plt.tight_layout()
        
        os.makedirs(os.path.dirname(args.out_correlation) or ".", exist_ok=True)
        plt.savefig(args.out_correlation, dpi=300)
        plt.close()
        print(f"Plot korelasi disimpan di: {args.out_correlation}")

if __name__ == "__main__":
    main()
