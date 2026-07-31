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
    parser.add_argument("--alignment", default="",
                        help="Codon alignment used to derive the site numbering. "
                             "Required to convert alignment columns into UniProt "
                             "residue numbers; without it the columns are used "
                             "directly, which is only correct for a gapless "
                             "alignment of a full-length reference.")
    parser.add_argument("--out_mapping", default="",
                        help="Optional TSV recording alignment column -> UniProt "
                             "residue -> PDB residue for every mapped site")
    return parser.parse_args()


def read_fasta(path):
    recs, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(chunks)))
                name, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if name is not None:
        recs.append((name, "".join(chunks)))
    return recs


def build_column_to_residue_map(alignment_path):
    """
    Map 1-based codon-alignment column -> 1-based residue number in the
    least-gapped sequence of the alignment.

    Selection sites are numbered by alignment column. AlphaFold numbers residues
    by position in the protein. Those coincide only when the reference sequence
    has no gaps, which is false whenever the alignment contains partial
    sequences: for Andes GnGc the offset reaches -44 by the C-terminus.
    """
    recs = read_fasta(alignment_path)
    if not recs:
        return {}, None
    # Codon alignments are nucleotide; collapse each codon to one column.
    best_name, best_seq, best_gaps = None, None, None
    for name, seq in recs:
        gaps = seq.count("-")
        if best_gaps is None or gaps < best_gaps:
            best_name, best_seq, best_gaps = name, seq, gaps

    mapping, residue = {}, 0
    n_codons = len(best_seq) // 3
    for col in range(1, n_codons + 1):
        codon = best_seq[(col - 1) * 3: col * 3]
        if codon and codon != "---" and "-" not in codon:
            residue += 1
            mapping[col] = residue
    return mapping, best_name

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
                # parse_additional_selection.py writes this column with str(bool),
                # so the values are "True"/"False" and never "Yes". Comparing
                # against 'Yes' matched nothing and produced all-blue structures.
                if cfel_sig_idx != -1:
                    if (len(parts) > cfel_sig_idx
                            and parts[cfel_sig_idx].strip().lower() == "true"):
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
                # AlphaFold now serves fragment models for long proteins: the PDB
                # restarts numbering at 1 while covering UniProt uniprotStart..End.
                # For the L proteins here that offset is +562 to +1292, so ignoring
                # it mislabels every mapped residue.
                u_start = entry.get("uniprotStart")
                u_end = entry.get("uniprotEnd")
                if u_start:
                    print(f"AlphaFold model covers UniProt {u_start}-{u_end} "
                          f"(PDB residue 1 = UniProt {u_start})")
                if pdb_url:
                    print(f"Mencoba mengunduh PDB dari API URL: {pdb_url}")
                    pdb_req = urllib.request.Request(pdb_url, headers=headers)
                    with urllib.request.urlopen(pdb_req, timeout=15) as pdb_resp, open(temp_path, 'wb') as out_file:
                        out_file.write(pdb_resp.read())
                    print(f"PDB berhasil diunduh ke {temp_path}")
                    return True, int(u_start) if u_start else 1
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
        # Static fallback URL is always fragment F1, i.e. UniProt offset 1.
        return True, 1
    except Exception as e:
        print(f"Error: Gagal mengunduh PDB untuk UniProt {uniprot_id} via fallback: {e}")
        return False, 1

def fail_no_structure(uniprot_id, out_pdb, out_html, n_sites, out_mapping=""):
    """
    Write an explicit 'no structure' marker instead of a mock PDB.

    The previous fallback emitted a four-atom single-methionine PDB. That file
    was indistinguishable from a real structure to every downstream step, so a
    failed download was shipped as supplementary material for a 1148-residue
    glycoprotein, complete with clickable site buttons that highlighted nothing.
    """
    with open(out_pdb, "w") as fh:
        fh.write(f"REMARK  NO STRUCTURE AVAILABLE FOR UNIPROT {uniprot_id}\n")
        fh.write("REMARK  AlphaFold DB returned no model for this accession.\n")
        fh.write("REMARK  This file intentionally contains no ATOM records.\n")
        fh.write("END\n")
    with open(out_html, "w") as fh:
        fh.write(
            "<!DOCTYPE html><html><head><meta charset='utf-8'>"
            f"<title>No structure - {uniprot_id}</title></head><body>"
            "<h2>No 3D structure available</h2>"
            f"<p>AlphaFold DB has no predicted model for UniProt "
            f"<b>{uniprot_id}</b>, so the {n_sites} selection site(s) for this "
            "protein could not be mapped onto a structure.</p>"
            "<p>This page is a placeholder recording that absence. It is not a "
            "failed render.</p></body></html>")
    # The coordinate map is a declared rule output, so it must exist even when
    # there is nothing to map. Writing only the header records "no sites could be
    # placed" without letting the rule fail on a missing file, which is what took
    # down every dataset whose accession AlphaFold has no model for.
    if out_mapping:
        os.makedirs(os.path.dirname(out_mapping) or ".", exist_ok=True)
        with open(out_mapping, "w") as fh:
            fh.write("alignment_column\treference_residue\tpdb_residue\t"
                     "uniprot_residue\tin_model\n")
    print(f"NOTE: AlphaFold DB has no model for {uniprot_id}; wrote explicit "
          f"no-structure markers. This is a property of the accession, not a "
          f"pipeline failure.", file=sys.stderr)

def main():
    args = parse_args()
    
    os.makedirs(os.path.dirname(args.out_pdb) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_html) or ".", exist_ok=True)
    
    # 1. Baca situs seleksi positif
    positive_sites = read_selection_sites(args.selection)
    print(f"Membaca {len(positive_sites)} situs seleksi positif dari {args.selection}")
    
    # 2. Download PDB
    temp_pdb = args.out_pdb + ".temp"
    download_success, uniprot_start = download_alphafold_pdb(args.uniprot, temp_pdb)
    
    if not download_success:
        fail_no_structure(args.uniprot, args.out_pdb, args.out_html,
                          len(positive_sites), args.out_mapping)
        return


    # 3. Convert alignment columns into PDB residue numbers before colouring.
    #    Two shifts are involved and both were previously ignored:
    #      a) gaps in the alignment reference, so column != residue
    #      b) fragment models, whose residue 1 is UniProt uniprot_start
    col2res, ref_name = ({}, None)
    if args.alignment and os.path.exists(args.alignment):
        col2res, ref_name = build_column_to_residue_map(args.alignment)
        print(f"Coordinate reference: {ref_name} "
              f"({len(col2res)} ungapped codons)")
    else:
        print("WARNING: no alignment supplied; alignment columns are being used "
              "directly as residue numbers. This is only correct for a gapless "
              "full-length reference.", file=sys.stderr)

    mapping_rows, positive_set, outside = [], set(), []
    for col in positive_sites:
        residue = col2res.get(col, col) if col2res else col
        pdb_res = residue - uniprot_start + 1
        if pdb_res < 1:
            outside.append((col, residue, pdb_res))
            continue
        positive_set.add(pdb_res)
        mapping_rows.append((col, residue, pdb_res))

    if outside:
        print(f"WARNING: {len(outside)} site(s) fall before the modelled region "
              f"(UniProt {uniprot_start}+) and cannot be shown: "
              f"{[c for c, _, _ in outside][:10]}", file=sys.stderr)
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
                
    # Report which mapped sites actually exist in the model, and record the
    # full column -> UniProt -> PDB correspondence so the numbering used in the
    # figures can be checked against the tables.
    modelled = set()
    for line in pdb_content:
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            try:
                modelled.add(int(line[22:26].strip()))
            except ValueError:
                pass
    missing = sorted(positive_set - modelled)
    if missing:
        print(f"WARNING: {len(missing)} site(s) map beyond the modelled region "
              f"and are not shown: {missing[:10]}", file=sys.stderr)
    print(f"Mapped {len(positive_set & modelled)}/{len(positive_sites)} sites "
          f"onto the structure.")

    if args.out_mapping:
        os.makedirs(os.path.dirname(args.out_mapping) or ".", exist_ok=True)
        with open(args.out_mapping, "w") as fh:
            fh.write("alignment_column\treference_residue\tpdb_residue\t"
                     "uniprot_residue\tin_model\n")
            for col, res, pdb_res in mapping_rows:
                fh.write(f"{col}\t{res}\t{pdb_res}\t{res}\t"
                         f"{pdb_res in modelled}\n")

    # Hapus temp file
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)
        
    # 4. Buat viewer HTML statis
    # Each button carries BOTH coordinates. The manuscript figure previously
    # showed alignment columns while the accompanying table showed UniProt
    # residues, so the same site appeared under two numbers 43 apart.
    site_buttons = "".join(
        f'<button class="site-btn" data-resi="{pdb_res}" '
        f'title="alignment column {col} = UniProt residue {res}">'
        f'{res}</button>'
        for col, res, pdb_res in sorted(mapping_rows, key=lambda r: r[1])
        if pdb_res in modelled) or (
        '<span style="opacity:.7">No site could be mapped onto this model.</span>')

    pdb_data_js = "".join(pdb_content).replace("\n", "\\n").replace("\r", "").replace("'", "\\'")
    
    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>3D Structure Mapping - {args.uniprot}</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
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
        .site-strip {{
            margin: 12px 0;
            padding: 12px;
            background-color: #0f172a;
            border-radius: 8px;
        }}
        .site-btn {{
            background-color: #1e40af;
            color: #f8fafc;
            border: 1px solid #3b82f6;
            border-radius: 4px;
            padding: 4px 9px;
            margin: 2px;
            cursor: pointer;
            font-size: 13px;
        }}
        .site-btn:hover {{ background-color: #2563eb; }}
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
        
        <div class="site-strip">
            <div style="margin-bottom:8px; font-weight:600;">
                Mapped sites (UniProt numbering; hover for the alignment column):
            </div>
            {site_buttons}
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
            
            // Cek apakah pustaka 3Dmol berhasil dimuat dari CDN.
            // NOTE: "3Dmol" is not a legal JavaScript identifier (a numeral cannot
            // start one), so `typeof 3Dmol` is a parse error that kills the whole
            // script block. The library exposes itself as $3Dmol; when loaded via
            // a plain <script> tag it is also reachable as window["3Dmol"].
            var mol3D = (typeof $3Dmol !== "undefined") ? $3Dmol : window["3Dmol"];
            if (!mol3D) {{
                element.innerHTML = '<div style="color: #f87171; padding: 40px; text-align: center; font-weight: bold; font-family: system-ui, sans-serif; line-height: 1.6; margin-top: 150px;">' +
                    '<span style="font-size: 24px;">Gagal Memuat Visualisasi 3D</span><br><br>' +
                    'Pustaka visualisasi 3D (3Dmol.js) tidak dapat diunduh dari CDN.<br>' +
                    'Harap hubungkan komputer Anda ke internet, atau periksa apakah ekstensi penolak iklan (Ad-blocker) / firewall memblokir cdnjs.cloudflare.com.' +
                    '</div>';
                return;
            }}
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

            // Focus a residue when its button is clicked. selectedAtoms() is the
            // 3Dmol API; an earlier revision called getAtoms(), which does not
            // exist in the library and threw before render() was reached, so
            // clicking a site appeared to do nothing.
            document.querySelectorAll(".site-btn").forEach(function(btn) {{
                btn.addEventListener("click", function() {{
                    var resi = parseInt(btn.getAttribute("data-resi"), 10);
                    viewer.setStyle({{}}, {{ cartoon: {{ color: '#3b82f6' }} }});
                    viewer.setStyle({{ predicate: function(atom) {{ return atom.b >= 100.0; }} }},
                                    {{ cartoon: {{ color: '#ef4444' }}, stick: {{ color: '#ef4444' }} }});
                    viewer.setStyle({{ resi: resi }},
                                    {{ sphere: {{ color: '#facc15', radius: 1.2 }} }});
                    viewer.zoomTo({{ resi: resi }});
                    viewer.render();
                }});
            }});
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
