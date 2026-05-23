import argparse
import json
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Parse HyPhy MEME, FEL, and FUBAR JSON outputs to extract selection sites")
    parser.add_argument("--meme", required=True, help="Input HyPhy MEME JSON file")
    parser.add_argument("--fel", required=True, help="Input HyPhy FEL JSON file")
    parser.add_argument("--fubar", required=True, help="Input HyPhy FUBAR JSON file")
    parser.add_argument("--pvalue", type=float, default=0.05, help="P-value threshold for MEME and FEL")
    parser.add_argument("--fubar-pp", type=float, default=0.90, help="Posterior probability threshold for FUBAR")
    parser.add_argument("--min-methods", type=int, default=2, help="Minimum number of methods for consensus")
    parser.add_argument("--out-consensus", required=True, help="Output TXT file for consensus sites")
    parser.add_argument("--out-meme", required=True, help="Output TXT file for MEME sites")
    parser.add_argument("--out-fel", required=True, help="Output TXT file for FEL sites")
    parser.add_argument("--out-fubar", required=True, help="Output TXT file for FUBAR sites")
    parser.add_argument("--out-complete", required=True, help="Output TSV file for all sites with detailed info")
    return parser.parse_args()

def load_json(filepath):
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return None
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
            if not data:
                return None
            return data
    except Exception as e:
        print(f"Warning: Gagal membaca {filepath}: {e}")
        return None

def parse_meme(data):
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    content = mle.get("content", [])
    
    p_idx = -1
    for i, h in enumerate(headers):
        h_str = h[0].lower() if isinstance(h, list) else h.lower()
        if "p-value" in h_str or "p value" in h_str or "pval" in h_str:
            p_idx = i
            break
    if p_idx == -1:
        p_idx = min(5, len(headers) - 1)
        
    results = {}
    for idx, row in enumerate(content):
        if len(row) > p_idx:
            try:
                results[idx + 1] = float(row[p_idx])
            except Exception:
                results[idx + 1] = 1.0
    return results

def parse_fel(data):
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    content = mle.get("content", [])
    
    p_idx = -1
    for i, h in enumerate(headers):
        h_str = h[0].lower() if isinstance(h, list) else h.lower()
        if "p-value" in h_str or "p value" in h_str or "pval" in h_str:
            p_idx = i
            break
    if p_idx == -1:
        p_idx = min(3, len(headers) - 1)
        
    results = {}
    for idx, row in enumerate(content):
        if len(row) > p_idx:
            try:
                results[idx + 1] = float(row[p_idx])
            except Exception:
                results[idx + 1] = 1.0
    return results

def parse_fubar(data):
    if not data:
        return {}
    
    fubar_data = None
    if "posterior" in data:
        fubar_data = data["posterior"]
    elif "MLE" in data:
        fubar_data = data["MLE"]
        
    if not fubar_data:
        return {}
        
    headers = fubar_data.get("headers", [])
    content = fubar_data.get("content", [])
    
    pp_idx = -1
    for i, h in enumerate(headers):
        h_str = h[0].lower() if isinstance(h, list) else h.lower()
        if "prob[alpha < beta]" in h_str or "prob[alpha<beta]" in h_str or "alpha<beta" in h_str or "alpha < beta" in h_str or "prob[beta > alpha]" in h_str or "beta > alpha" in h_str:
            pp_idx = i
            break
    if pp_idx == -1:
        pp_idx = min(4, len(headers) - 1)
        
    results = {}
    for idx, row in enumerate(content):
        if len(row) > pp_idx:
            try:
                results[idx + 1] = float(row[pp_idx])
            except Exception:
                results[idx + 1] = 0.0
    return results

def main():
    args = parse_args()
    
    meme_data = load_json(args.meme)
    fel_data = load_json(args.fel)
    fubar_data = load_json(args.fubar)
    
    meme_res = parse_meme(meme_data)
    fel_res = parse_fel(fel_data)
    fubar_res = parse_fubar(fubar_data)
    
    total_sites = 0
    if meme_res:
        total_sites = max(total_sites, max(meme_res.keys()))
    if fel_res:
        total_sites = max(total_sites, max(fel_res.keys()))
    if fubar_res:
        total_sites = max(total_sites, max(fubar_res.keys()))
        
    for out_file in [args.out_consensus, args.out_meme, args.out_fel, args.out_fubar, args.out_complete]:
        out_dir = os.path.dirname(out_file)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
            
    # Tulis hasil individual (Hanya yang signifikan untuk kompatibilitas hilir)
    meme_count = 0
    with open(args.out_meme, 'w') as f:
        f.write("Site\tP-value\n")
        for site in range(1, total_sites + 1):
            pval = meme_res.get(site, 1.0)
            if pval < args.pvalue:
                f.write(f"{site}\t{pval:.6f}\n")
                meme_count += 1
                
    fel_count = 0
    with open(args.out_fel, 'w') as f:
        f.write("Site\tP-value\n")
        for site in range(1, total_sites + 1):
            pval = fel_res.get(site, 1.0)
            if pval < args.pvalue:
                f.write(f"{site}\t{pval:.6f}\n")
                fel_count += 1
                
    fubar_count = 0
    with open(args.out_fubar, 'w') as f:
        f.write("Site\tPosteriorProb\n")
        for site in range(1, total_sites + 1):
            pp = fubar_res.get(site, 0.0)
            if pp >= args.fubar_pp:
                f.write(f"{site}\t{pp:.6f}\n")
                fubar_count += 1
                
    # Tulis file konsensus (Hanya yang signifikan untuk kompatibilitas hilir)
    consensus_count = 0
    with open(args.out_consensus, 'w') as f:
        f.write("Site\tP-value\n")
        for site in range(1, total_sites + 1):
            meme_p = meme_res.get(site, 1.0)
            fel_p = fel_res.get(site, 1.0)
            fubar_pp = fubar_res.get(site, 0.0)
            
            sig_meme = 1 if meme_p < args.pvalue else 0
            sig_fel = 1 if fel_p < args.pvalue else 0
            sig_fubar = 1 if fubar_pp >= args.fubar_pp else 0
            
            n_methods = sig_meme + sig_fel + sig_fubar
            
            if n_methods >= args.min_methods:
                # Gunakan p-value rata-rata dari MEME dan FEL sebagai representasi P-value
                avg_p = (meme_p + fel_p) / 2.0
                f.write(f"{site}\t{avg_p:.6f}\n")
                consensus_count += 1
                
    # Tulis hasil lengkap untuk semua site (Analisis lanjutan oleh user)
    with open(args.out_complete, 'w') as f:
        f.write("Site\tMEME_p\tFEL_p\tFUBAR_pp\tn_methods\tConsensus\n")
        for site in range(1, total_sites + 1):
            meme_p = meme_res.get(site, 1.0)
            fel_p = fel_res.get(site, 1.0)
            fubar_pp = fubar_res.get(site, 0.0)
            
            sig_meme = 1 if meme_p < args.pvalue else 0
            sig_fel = 1 if fel_p < args.pvalue else 0
            sig_fubar = 1 if fubar_pp >= args.fubar_pp else 0
            
            n_methods = sig_meme + sig_fel + sig_fubar
            is_consensus = 1 if n_methods >= args.min_methods else 0
            
            f.write(f"{site}\t{meme_p:.6f}\t{fel_p:.6f}\t{fubar_pp:.6f}\t{n_methods}\t{is_consensus}\n")
            
    print(f"Parser selesai. Total situs: {total_sites}")
    print(f"MEME signifikan: {meme_count}, FEL signifikan: {fel_count}, FUBAR signifikan: {fubar_count}")
    print(f"Konsensus signifikan (>= {args.min_methods} metode): {consensus_count}")

if __name__ == "__main__":
    main()
