import argparse
import json
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Parse HyPhy MEME JSON output to extract positive selection sites")
    parser.add_argument("--json", required=True, help="Input HyPhy MEME JSON file")
    parser.add_argument("--pvalue", type=float, default=0.05, help="P-value threshold")
    parser.add_argument("--out", required=True, help="Output TXT file")
    return parser.parse_args()

def main():
    args = parse_args()
    
    if not os.path.exists(args.json):
        print(f"Error: File {args.json} tidak ditemukan.")
        # Buat output kosong jika file JSON tidak terbuat (misal karena sequence homogen)
        # agar snakemake tidak error dan workflow tetap berlanjut
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, 'w') as f:
            f.write("Site\tP-value\n")
        sys.exit(0)
        
    try:
        with open(args.json, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error: Gagal membaca file JSON: {e}")
        # Fallback file kosong
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, 'w') as f:
            f.write("Site\tP-value\n")
        sys.exit(0)
        
    mle = data.get("MLE", {})
    headers = mle.get("headers", [])
    content = mle.get("content", [])
    
    # Cari indeks untuk p-value
    p_idx = -1
    for i, h in enumerate(headers):
        # Cari header yang merepresentasikan p-value
        h_lower = h[0].lower() if isinstance(h, list) else h.lower()
        if "p-value" in h_lower or "p value" in h_lower or "pval" in h_lower:
            p_idx = i
            break
            
    if p_idx == -1:
        # Fallback default: kolom ke-5 biasanya adalah p-value untuk MEME
        # Kolom typical: [alpha, beta-, beta+, proportion, LRT, p-value]
        p_idx = min(5, len(headers) - 1)
        
    selected_sites = []
    # Loop over all sites
    for site_idx, row in enumerate(content):
        if len(row) > p_idx:
            p_val = row[p_idx]
            try:
                p_val = float(p_val)
                if p_val < args.pvalue:
                    # 1-indexed site position untuk biologi
                    selected_sites.append((site_idx + 1, p_val))
            except Exception:
                continue
                
    # Tulis hasil
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as out_f:
        out_f.write("Site\tP-value\n")
        for site, pval in selected_sites:
            out_f.write(f"{site}\t{pval:.6f}\n")
            
    print(f"Parser berhasil: Mengekstrak {len(selected_sites)} situs terseleksi dengan p-value < {args.pvalue}")

if __name__ == "__main__":
    main()
