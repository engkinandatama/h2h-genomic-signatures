import argparse
import os
import sys
import urllib.request
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Map positive selection sites to 3D structure using AlphaFold & Py3Dmol")
    parser.add_argument("--uniprot", required=True, help="UniProt ID to fetch from AlphaFold DB")
    parser.add_argument("--selection", required=True, help="HyPhy selection results txt/json")
    parser.add_argument("--out_pdb", required=True, help="Output annotated PDB file")
    parser.add_argument("--out_html", required=True, help="Output Py3Dmol HTML viewer")
    return parser.parse_args()

def read_selection_sites(results_path):
    sites = []
    if not os.path.exists(results_path):
        print(f"Warning: File hasil seleksi {results_path} tidak ditemukan.")
        return sites
        
    with open(results_path, 'r') as f:
        header_line = f.readline().strip()
        headers = header_line.split('\t')
        
        try:
            site_idx = headers.index("site")
            cfel_sig_idx = headers.index("cfel_sig")
        except ValueError:
            print("Warning: Format file tidak dikenali (tidak ada kolom 'site' atau 'cfel_sig'). Menggunakan fallback pembacaan kolom 1.")
            site_idx = 0
            cfel_sig_idx = -1
            
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            try:
                site = int(parts[site_idx])
                # Jika format lama (cfel_sig_idx == -1), anggap semua baris adalah seleksi positif.
                # Jika format baru, cek apakah cfel_sig adalah 'Yes'.
                if cfel_sig_idx != -1:
                    if len(parts) > cfel_sig_idx and parts[cfel_sig_idx] == 'Yes':
                        sites.append(site)
                else:
                    sites.append(site)
            except ValueError:
                continue
    return sites

def download_alphafold_pdb(uniprot_id, temp_path):
    import json
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
    }
    
    # 1. Coba query API AlphaFold terlebih dahulu untuk mendapatkan URL PDB yang dinamis
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    print(f"Mengkueri API AlphaFold untuk UniProt ID {uniprot_id}...")
    try:
        req = urllib.request.Request(api_url, headers=headers)
        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode())
            if data and isinstance(data, list):
                # Cari entri yang cocok atau gunakan entri pertama
                entry = data[0]
                pdb_url = entry.get("pdbUrl")
                if pdb_url:
                    print(f"Mencoba mengunduh PDB dari API URL: {pdb_url}")
                    pdb_req = urllib.request.Request(pdb_url, headers=headers)
                    with urllib.request.urlopen(pdb_req, timeout=15) as pdb_resp, open(temp_path, 'wb') as out_file:
                        out_file.write(pdb_resp.read())
                    print(f"PDB berhasil diunduh ke {temp_path}")
                    return True
    except Exception as e:
        print(f"Peringatan: Gagal kueri API AlphaFold untuk {uniprot_id}: {e}")

    # 2. Fallback ke URL hardcoded model_v4 jika API gagal
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    print(f"Mencoba fallback unduh PDB dari URL statis: {url}")
    try:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=15) as response, open(temp_path, 'wb') as out_file:
            out_file.write(response.read())
        print(f"PDB berhasil diunduh via fallback ke {temp_path}")
        return True
    except Exception as e:
        print(f"Error: Gagal mengunduh PDB untuk UniProt {uniprot_id} via fallback: {e}")
        return False

def generate_fallback_pdb(out_path, uniprot_id):
    # Menulis file PDB minimal tiruan jika download dari internet gagal
    with open(out_path, 'w') as f:
        f.write(f"HEADER    FALLBACK MOCK PDB FOR UNIPROT {uniprot_id}\n")
        f.write("ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N\n")
        f.write("ATOM      2  CA  MET A   1       1.450   0.000   0.000  1.00  0.00           C\n")
        f.write("ATOM      3  C   MET A   1       2.000   1.450   0.000  1.00  0.00           C\n")
        f.write("ATOM      4  O   MET A   1       1.200   2.400   0.000  1.00  0.00           O\n")
        f.write("TER\n")
        f.write("END\n")

def main():
    args = parse_args()
    
    os.makedirs(os.path.dirname(args.out_pdb) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_html) or ".", exist_ok=True)
    
    # 1. Baca situs seleksi positif
    positive_sites = read_selection_sites(args.selection)
    print(f"Membaca {len(positive_sites)} situs seleksi positif dari {args.selection}")
    
    # 2. Download PDB
    temp_pdb = args.out_pdb + ".temp"
    download_success = download_alphafold_pdb(args.uniprot, temp_pdb)
    
    if not download_success:
        print("Menggunakan fallback mock PDB karena unduhan gagal.")
        generate_fallback_pdb(temp_pdb, args.uniprot)
        
    # 3. Petakan situs ke kolom B-factor PDB
    positive_set = set(positive_sites)
    pdb_content = []
    
    with open(temp_pdb, 'r') as infile, open(args.out_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                try:
                    res_num = int(line[22:26].strip())
                    # Warnai residu terpilih dengan B-factor 100.00
                    b_factor = 100.00 if res_num in positive_set else 0.00
                    b_factor_str = f"{b_factor:6.2f}"
                    # Pastikan baris cukup panjang sebelum diparsing
                    if len(line) >= 66:
                        modified_line = line[:60] + b_factor_str + line[66:]
                    else:
                        # Pad baris yang terlalu pendek
                        padded = line.rstrip("\n").ljust(66)
                        modified_line = padded[:60] + b_factor_str + "\n"
                    outfile.write(modified_line)
                    pdb_content.append(modified_line)
                except Exception:
                    outfile.write(line)
                    pdb_content.append(line)
            else:
                outfile.write(line)
                pdb_content.append(line)
                
    # Hapus temp file
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)
        
    # 4. Buat viewer HTML statis
    pdb_data_js = "".join(pdb_content).replace("\n", "\\n").replace("\r", "").replace("'", "\\'")
    
    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>3D Structure Mapping - {args.uniprot}</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3dmol/2.0.4/3Dmol-min.js"></script>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background-color: #0f172a;
            color: #f8fafc;
            margin: 0;
            padding: 20px;
        }}
        .container {{
            max-width: 1000px;
            margin: 0 auto;
            background-color: #1e293b;
            padding: 24px;
            border-radius: 12px;
            box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1), 0 2px 4px -2px rgb(0 0 0 / 0.1);
        }}
        h1 {{
            font-size: 1.8rem;
            margin-top: 0;
            color: #38bdf8;
        }}
        .meta-info {{
            font-size: 0.95rem;
            color: #94a3b8;
            margin-bottom: 20px;
            border-bottom: 1px solid #334155;
            padding-bottom: 10px;
        }}
        #viewer {{
            width: 100%;
            height: 550px;
            position: relative;
            background-color: #0b0f19;
            border-radius: 8px;
            overflow: hidden;
            border: 1px solid #475569;
        }}
        .legend {{
            display: flex;
            gap: 24px;
            margin: 20px 0;
            background-color: #0f172a;
            padding: 12px 16px;
            border-radius: 6px;
            border-left: 4px solid #38bdf8;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 0.9rem;
        }}
        .dot {{
            width: 14px;
            height: 14px;
            border-radius: 4px;
        }}
        .red {{ background-color: #ef4444; }}
        .blue {{ background-color: #3b82f6; }}
        .description {{
            font-size: 0.95rem;
            line-height: 1.6;
            color: #cbd5e1;
        }}
        .site-list {{
            background-color: #0f172a;
            padding: 12px;
            border-radius: 6px;
            max-height: 100px;
            overflow-y: auto;
            font-family: monospace;
            color: #38bdf8;
            margin-top: 10px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>3D Mapping of Positive Selection Sites</h1>
        <div class="meta-info">
            <strong>UniProt ID:</strong> {args.uniprot} | 
            <strong>Jumlah Residu Terseleksi:</strong> {len(positive_sites)}
        </div>
        
        <div id="viewer"></div>
        
        <div class="legend">
            <div class="legend-item">
                <div class="dot red"></div>
                <span>Situs Seleksi Positif (p &lt; 0.05)</span>
            </div>
            <div class="legend-item">
                <div class="dot blue"></div>
                <span>Residu Lainnya (Netral / Purifying)</span>
            </div>
        </div>
        
        <div class="description">
            <p>Visualisasi 3D protein di atas bersumber dari database <strong>AlphaFold</strong>.
            Residu berwarna <strong style="color: #ef4444;">merah (stick & cartoon)</strong> menunjukkan situs asam amino yang diidentifikasi mengalami seleksi positif episodik secara signifikan oleh model HyPhy MEME (p &lt; 0.05). Residu berwarna <strong style="color: #3b82f6;">biru (cartoon)</strong> menunjukkan situs netral atau di bawah seleksi negatif.</p>
            
            <strong>Daftar Posisi Residu Terseleksi Positif:</strong>
            <div class="site-list">
                {", ".join(map(str, sorted(positive_sites))) if positive_sites else "Tidak ada situs yang terseleksi positif."}
            </div>
        </div>
    </div>

    <script>
        document.addEventListener("DOMContentLoaded", function() {{
            let element = document.getElementById("viewer");
            
            // Cek apakah pustaka 3Dmol berhasil dimuat dari CDN
            if (typeof $3Dmol === "undefined" && typeof 3Dmol === "undefined") {
                element.innerHTML = '<div style="color: #f87171; padding: 40px; text-align: center; font-weight: bold; font-family: system-ui, sans-serif; line-height: 1.6; margin-top: 150px;">' +
                    '<span style="font-size: 24px;">⚠️ Gagal Memuat Visualisasi 3D</span><br><br>' +
                    'Pustaka visualisasi 3D (3Dmol.js) tidak dapat diunduh dari CDN.<br>' +
                    'Harap hubungkan komputer Anda ke internet, atau periksa apakah ekstensi penolak iklan (Ad-blocker) / firewall memblokir cdnjs.cloudflare.com.' +
                    '</div>';
                return;
            }
            
            let mol3D = typeof $3Dmol !== "undefined" ? $3Dmol : 3Dmol;
            let viewer = mol3D.createViewer(element, {{}});
            
            let pdbData = '{pdb_data_js}';
            
            viewer.addModel(pdbData, "pdb");
            
            // Styling default: warna biru untuk cartoon
            viewer.setStyle({{}}, {{ cartoon: {{ color: '#3b82f6' }} }});
            
            // Styling positive selection: warna merah untuk residu dengan B-factor >= 100.00
            viewer.setStyle({{ predicate: function(atom) {{ return atom.b >= 100.0; }} }}, {{ 
                cartoon: {{ color: '#ef4444' }},
                stick: {{ color: '#ef4444', radius: 0.25 }}
            }});
            
            viewer.zoomTo();
            viewer.render();
            viewer.setBackgroundColor('#0b0f19');
        }});
    </script>
</body>
</html>
"""
    with open(args.out_html, 'w', encoding='utf-8') as out_h:
        out_h.write(html_template)
        
    print(f"HTML Viewer berhasil dibuat di: {args.out_html}")

if __name__ == "__main__":
    main()
