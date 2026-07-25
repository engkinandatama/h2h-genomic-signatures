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

def _extract_rows_from_content(content):
    """Handle HyPhy content, which is a dict {"0": [rows...]} or a plain list."""
    if isinstance(content, dict):
        # Standard HyPhy format: content = {"0": [row1, row2, ...]}
        key = list(content.keys())[0]
        return content[key]
    elif isinstance(content, list):
        return content
    return []


def _find_header_index(headers, keywords, default_idx):
    """Search headers (may be list of [name, desc] or plain strings) for a keyword."""
    for i, h in enumerate(headers):
        h_str = (h[0] + " " + h[1]).lower() if isinstance(h, list) and len(h) >= 2 else str(h).lower()
        if any(kw in h_str for kw in keywords):
            return i
    return default_idx


def parse_meme(data):
    if not data or "MLE" not in data:
        return {}
    mle = data["MLE"]
    headers = mle.get("headers", [])
    content = mle.get("content", {})

    # MEME p-value column: look for "p-value" in header description, default index 6
    p_idx = _find_header_index(headers, ["p-value", "p value", "pval"], 6)

    rows = _extract_rows_from_content(content)
    results = {}
    for idx, row in enumerate(rows):
        if len(row) > p_idx:
            try:
                results[idx + 1] = float(row[p_idx])
            except Exception:
                results[idx + 1] = 1.0
    return results

def parse_fel(data):
    """
    Return (p_values, positive_sites).

    FEL's p-value tests the two-sided hypothesis beta != alpha, so a low p-value
    alone does not mean positive selection: sites under strong purifying selection
    are equally significant. positive_sites holds only those sites where beta >
    alpha, and callers must require membership in it before counting a site as
    evidence of positive selection.
    """
    if not data or "MLE" not in data:
        return {}, set()
    mle = data["MLE"]
    headers = mle.get("headers", [])
    content = mle.get("content", {})

    # FEL p-value column: look for "p-value" in header description, default index 4
    p_idx = _find_header_index(headers, ["p-value", "p value", "pval"], 4)
    a_idx = _find_header_index(headers, ["alpha"], 0)
    b_idx = _find_header_index(headers, ["beta"], 1)

    rows = _extract_rows_from_content(content)
    results = {}
    positive_sites = set()
    for idx, row in enumerate(rows):
        site = idx + 1
        if len(row) > p_idx:
            try:
                results[site] = float(row[p_idx])
            except Exception:
                results[site] = 1.0
        if len(row) > max(a_idx, b_idx):
            try:
                if float(row[b_idx]) > float(row[a_idx]):
                    positive_sites.add(site)
            except Exception:
                pass
    return results, positive_sites

def parse_fubar(data):
    if not data:
        return {}

    # FUBAR stores results in MLE section
    fubar_data = None
    if "MLE" in data:
        fubar_data = data["MLE"]
    elif "posterior" in data:
        fubar_data = data["posterior"]

    if not fubar_data:
        return {}

    headers = fubar_data.get("headers", [])
    content = fubar_data.get("content", {})

    # FUBAR positive selection column: Prob[alpha<beta], default index 4
    pp_idx = _find_header_index(
        headers,
        ["prob[alpha<beta]", "prob[alpha < beta]", "alpha<beta", "alpha < beta",
         "prob[beta>alpha]", "prob[beta > alpha]", "positive selection"],
        4
    )

    rows = _extract_rows_from_content(content)
    results = {}
    for idx, row in enumerate(rows):
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
    fel_res, fel_pos = parse_fel(fel_data)
    fubar_res = parse_fubar(fubar_data)
    
    total_sites = 0
    if meme_res:
        total_sites = max(total_sites, max(meme_res.keys(), default=0))
    if fel_res:
        total_sites = max(total_sites, max(fel_res.keys(), default=0))
    if fubar_res:
        total_sites = max(total_sites, max(fubar_res.keys(), default=0))
        
    # Jika tidak ada situs sama sekali (semua JSON kosong), buat output kosong dan exit
    if total_sites == 0:
        for out_file in [args.out_consensus, args.out_meme, args.out_fel, args.out_fubar]:
            out_dir = os.path.dirname(out_file)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            with open(out_file, 'w') as f:
                f.write("Site\tP-value\n")
        out_dir = os.path.dirname(args.out_fubar)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out_fubar, 'w') as f:
            f.write("Site\tPosteriorProb\n")
        out_dir = os.path.dirname(args.out_complete)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out_complete, 'w') as f:
            f.write("Site\tMEME_p\tFEL_p\tFUBAR_pp\tn_methods\tConsensus\n")
        print("Warning: Semua JSON input kosong atau tidak valid. File output kosong dibuat.")
        sys.exit(0)
        
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
            if pval < args.pvalue and site in fel_pos:
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
            sig_fel = 1 if (fel_p < args.pvalue and site in fel_pos) else 0
            sig_fubar = 1 if fubar_pp >= args.fubar_pp else 0

            n_methods = sig_meme + sig_fel + sig_fubar

            if n_methods >= args.min_methods:
                # Gunakan p-value rata-rata dari MEME dan FEL sebagai representasi P-value
                avg_p = (meme_p + fel_p) / 2.0
                f.write(f"{site}\t{avg_p:.6f}\n")
                consensus_count += 1
                
    # Tulis hasil lengkap untuk semua site (Analisis lanjutan oleh user)
    with open(args.out_complete, 'w') as f:
        f.write("Site\tMEME_p\tFEL_p\tFEL_positive\tFUBAR_pp\tn_methods\tConsensus\n")
        for site in range(1, total_sites + 1):
            meme_p = meme_res.get(site, 1.0)
            fel_p = fel_res.get(site, 1.0)
            fel_is_pos = site in fel_pos
            fubar_pp = fubar_res.get(site, 0.0)

            sig_meme = 1 if meme_p < args.pvalue else 0
            sig_fel = 1 if (fel_p < args.pvalue and fel_is_pos) else 0
            sig_fubar = 1 if fubar_pp >= args.fubar_pp else 0

            n_methods = sig_meme + sig_fel + sig_fubar
            is_consensus = 1 if n_methods >= args.min_methods else 0

            f.write(f"{site}\t{meme_p:.6f}\t{fel_p:.6f}\t{int(fel_is_pos)}"
                    f"\t{fubar_pp:.6f}\t{n_methods}\t{is_consensus}\n")
            
    print(f"Parser selesai. Total situs: {total_sites}")
    print(f"MEME signifikan: {meme_count}, FEL signifikan: {fel_count}, FUBAR signifikan: {fubar_count}")
    print(f"Konsensus signifikan (>= {args.min_methods} metode): {consensus_count}")

if __name__ == "__main__":
    main()
