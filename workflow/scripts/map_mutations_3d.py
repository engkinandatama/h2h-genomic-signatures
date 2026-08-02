import argparse
import collections
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
    parser.add_argument("--site-columns", dest="site_columns",
                        default="episodic_only,consensus_pos,cfel_sig",
                        help="Comma-separated boolean columns of the all_sites "
                             "table whose True rows are highlighted. Later names "
                             "win when a site is flagged by more than one.")
    parser.add_argument("--virus-group", dest="virus_group", default="",
                        help="Group name, used in the page title and header.")
    parser.add_argument("--protein", dest="protein", default="",
                        help="Protein name, used in the page title and header.")
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

def read_selection_sites(results_path, columns):
    """
    Return {alignment_column: category} for every site flagged in `columns`.

    The viewer used to key on cfel_sig alone. Contrast-FEL significance is empty
    in all but one dataset here, and that dataset has no AlphaFold model, so every
    published viewer highlighted nothing at all while the tables listed sites. The
    default set now carries all three evidence tiers, weakest first so that the
    strongest claim wins a site flagged by more than one: episodic_only (MEME
    alone, no pervasive support), consensus_pos (agreed by at least min_methods
    site models), and cfel_sig (rates differ between branch sets). Each is drawn
    in its own colour rather than merged, so the viewer never implies more
    evidence for a site than the tables hold.
    """
    sites = {}
    if not os.path.exists(results_path):
        print(f"Warning: File hasil seleksi {results_path} tidak ditemukan.")
        return sites

    with open(results_path, "r") as f:
        headers = f.readline().rstrip("\n").split("\t")
        if "site" not in headers:
            print("Warning: kolom 'site' tidak ada; tidak ada situs yang dipetakan.",
                  file=sys.stderr)
            return sites
        wanted = [c for c in columns if c in headers]
        missing = [c for c in columns if c not in headers]
        if missing:
            print(f"Warning: kolom {missing} tidak ada di {results_path}.",
                  file=sys.stderr)
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            row = dict(zip(headers, line.split("\t")))
            try:
                site = int(row["site"])
            except (KeyError, ValueError):
                continue
            for col in wanted:
                if str(row.get(col, "")).strip().lower() == "true":
                    sites[site] = col
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
    site_columns = [c.strip() for c in args.site_columns.split(",") if c.strip()]
    positive_sites = read_selection_sites(args.selection, site_columns)
    print(f"Membaca {len(positive_sites)} situs terseleksi dari {args.selection} "
          f"(kolom: {', '.join(site_columns)})")
    
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

    mapping_rows, category_of, outside = [], {}, []
    for col, category in sorted(positive_sites.items()):
        residue = col2res.get(col, col) if col2res else col
        pdb_res = residue - uniprot_start + 1
        if pdb_res < 1:
            outside.append((col, residue, pdb_res))
            continue
        category_of[pdb_res] = category
        mapping_rows.append((col, residue, pdb_res, category))
    positive_set = set(category_of)

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
                    # Three tiers so the viewer never conflates the claims:
                    # 100 = differential between branch sets (Contrast-FEL),
                    #  60 = consensus positive selection across site models,
                    #  30 = episodic only, MEME with no pervasive support.
                    b_factor = {"cfel_sig": 100.00,
                                "consensus_pos": 60.00,
                                "episodic_only": 30.00}.get(
                                    category_of.get(res_num), 0.00)
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
                     "uniprot_residue\tcategory\tin_model\n")
            for col, res, pdb_res, category in mapping_rows:
                fh.write(f"{col}\t{res}\t{pdb_res}\t{res}\t{category}\t"
                         f"{pdb_res in modelled}\n")

    # Hapus temp file
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)

    # 4. Static viewer, matching the layout agreed for the supplementary material.
    #    The PDB is embedded in a <script type="text/plain"> block and read back
    #    with .text(), rather than escaped into a JavaScript string literal. A
    #    coordinate file contains backslashes and quotes; escaping it by hand is
    #    a source of silent corruption that only shows as an empty viewer.
    CATEGORY_STYLE = {
        "cfel_sig":      ("#ef4444", 100.0, "Differential selection between clades (Contrast-FEL)"),
        "consensus_pos": ("#f59e0b",  60.0, "Positive selection agreed by two or more site models"),
        "episodic_only": ("#a855f7",  30.0, "Episodic only (MEME, no pervasive support)"),
    }

    shown = [r for r in sorted(mapping_rows, key=lambda r: r[1]) if r[2] in modelled]
    site_buttons = "".join(
        f'<button class="site-btn" data-resi="{pdb_res}" data-category="{category}" '
        f'title="alignment column {col} = UniProt residue {res} ({category})">'
        f'Site {res}</button>'
        for col, res, pdb_res, category in shown) or (
        '<span style="opacity:.7">No site could be mapped onto this model.</span>')

    counts = collections.Counter(c for _, _, _, c in shown)
    legend_items = "".join(
        f'''
            <div class="legend-item">
                <span class="dot" style="background: {colour};"></span>
                <span><b>{label}:</b> {desc} &mdash; {counts.get(cat, 0)} site(s).</span>
            </div>'''
        for cat, (colour, _b, desc), label in (
            ("cfel_sig", CATEGORY_STYLE["cfel_sig"], "Differential"),
            ("consensus_pos", CATEGORY_STYLE["consensus_pos"], "Consensus"),
            ("episodic_only", CATEGORY_STYLE["episodic_only"], "Episodic only")))

    style_rules = "\n".join(
        f"            viewer.addStyle({{ b: {b} }}, {{ cartoon: {{ color: '{colour}' }}, "
        f"stick: {{ color: '{colour}', radius: 0.3 }} }});"
        for _cat, (colour, b, _d) in CATEGORY_STYLE.items())
    style_rules_focus = "\n".join(
        f"            viewer.addStyle({{ b: {b} }}, {{ cartoon: {{ color: '{colour}' }}, "
        f"stick: {{ color: '{colour}', radius: 0.2 }} }});"
        for _cat, (colour, b, _d) in CATEGORY_STYLE.items())

    group_label = (args.virus_group or "").replace("_", " ") or "Virus group"
    protein_label = (args.protein or "").replace("_", " ") or "Protein"
    pdb_block = "".join(pdb_content).replace("</", "<\\/")

    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>3D Mapping - {args.virus_group} {args.protein}</title>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background-color: #0f172a; color: #f8fafc; padding: 20px; line-height: 1.6; }}
        .container {{ max-width: 1000px; margin: 0 auto; background-color: #1e293b; padding: 30px; border-radius: 12px; box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.3); }}
        h1 {{ color: #38bdf8; margin-top: 0; font-size: 1.8rem; border-bottom: 2px solid #334155; padding-bottom: 10px; }}

        #info_panel {{
            background: #0f172a; border-left: 5px solid #ef4444; padding: 20px; margin: 20px 0; border-radius: 8px;
            font-size: 1.2rem; min-height: 60px; display: flex; align-items: center; justify-content: center;
            box-shadow: inset 0 2px 4px 0 rgba(0, 0, 0, 0.06);
        }}

        #viewer_container {{ position: relative; width: 100%; height: 550px; background-color: #0b0f19; border-radius: 10px; border: 1px solid #475569; margin-bottom: 25px; }}
        #viewer {{ width: 100%; height: 100%; }}

        .legend {{ display: flex; flex-direction: column; gap: 10px; margin: 20px 0; padding: 15px; background: #0f172a; border-radius: 8px; }}
        .legend-item {{ display: flex; align-items: center; gap: 10px; font-size: 0.95rem; }}
        .dot {{ width: 16px; height: 16px; border-radius: 4px; display: inline-block; flex: none; }}

        .section-title {{ color: #38bdf8; font-weight: bold; font-size: 1.1rem; margin-bottom: 10px; display: block; }}
        .description-box {{ background: #0f172a; padding: 20px; border-radius: 8px; margin-top: 25px; color: #cbd5e1; font-size: 0.95rem; }}

        .site-list {{ background-color: #0f172a; padding: 15px; border-radius: 8px; max-height: 180px; overflow-y: auto; border: 1px solid #334155; display: flex; flex-wrap: wrap; gap: 8px; }}
        .site-btn {{
            background: #1e293b; color: #38bdf8; border: 1px solid #38bdf8; padding: 6px 14px; border-radius: 6px;
            cursor: pointer; font-weight: bold; transition: all 0.2s; font-size: 0.9rem;
        }}
        .site-btn:hover {{ background: #38bdf8; color: #0f172a; transform: translateY(-2px); }}
        .site-btn.active {{ background: #facc15; color: #0f172a; border-color: #facc15; box-shadow: 0 0 10px rgba(250, 204, 21, 0.4); }}

        b {{ color: #f8fafc; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>3D Structural Mapping: {group_label}</h1>

        <div style="margin-bottom: 20px; font-size: 1rem; color: #94a3b8;">
            <b>Virus Group:</b> {group_label} &nbsp;|&nbsp;
            <b>Protein:</b> {protein_label} &nbsp;|&nbsp;
            <b>UniProt ID:</b> <a href="https://www.uniprot.org/uniprotkb/{args.uniprot}" target="_blank" style="color: #38bdf8;">{args.uniprot}</a>
        </div>

        <div id="info_panel">
            <span id="info_content" style="color: #94a3b8;"><i>Select a site below to begin structural inspection</i></span>
        </div>

        <div id="viewer_container">
            <div id="viewer"></div>
        </div>

        <div class="legend">{legend_items}
            <div class="legend-item">
                <span class="dot" style="background: #3b82f6;"></span>
                <span><b>Other residues:</b> protein backbone, neutral or under purifying selection.</span>
            </div>
        </div>

        <span class="section-title">Interactive Site Inspection:</span>
        <div class="site-list">{site_buttons}</div>

        <div class="description-box">
            <span class="section-title" style="color: #f8fafc;">About this Visualization:</span>
            <p>This interactive viewer maps sites under positive selection onto a 3D protein
            structure. Three tiers of evidence are drawn in separate colours rather than merged,
            so a site is never shown as better supported than the underlying tables allow. A site
            flagged by more than one test takes the colour of its strongest claim.</p>

            <p><b>Differential</b> sites have significantly different non-synonymous rates between
            the foreground and reference clades. <b>Consensus</b> sites are called under positive
            selection by at least two independent site models. <b>Episodic only</b> sites are
            significant under MEME alone, with no support from the pervasive models, and are the
            weakest of the three.</p>

            <p>Residue numbering follows the <b>UniProt</b> entry, not the alignment column, and the
            two differ wherever the reference carries a gap. Each button's tooltip gives both.
            Where the deposited model covers only part of the protein, sites outside the modelled
            region cannot be shown and are listed in the accompanying coordinate map.</p>

            <p style="font-size: 0.85rem; border-top: 1px solid #334155; padding-top: 10px; margin-top: 15px;">
                <b>Controls:</b> Left-click to Rotate | Right-click to Pan | Scroll to Zoom | Click buttons above to auto-focus.
            </p>
        </div>
    </div>

    <script type="text/plain" id="pdb_data">{pdb_block}</script>
    <script>
        var viewer;

        function baseStyle(opacity) {{
            viewer.setStyle({{}}, {{ cartoon: {{ color: '#3b82f6', opacity: opacity }} }});
        }}

        function zoomToSite(resi) {{
            let spec = {{ resi: parseInt(resi) }};
            viewer.zoomTo(spec, 1000);

            baseStyle(0.5);
{style_rules_focus}
            viewer.addStyle(spec, {{ stick: {{ color: '#facc15', radius: 0.5 }}, sphere: {{ color: '#facc15', radius: 1.0 }} }});

            let atoms = viewer.getModel().selectedAtoms(spec);
            if (atoms.length > 0) {{
                let a = atoms.find(x => x.name === 'CA') || atoms[0];
                let btn = $(".site-btn[data-resi='" + resi + "']");
                let cat = btn.data("category") || "";
                let infoHtml = "<b style='color: #facc15;'>Inspecting Site " + resi + "</b> &nbsp;|&nbsp; " +
                               "<b>Residue:</b> " + a.resn + " &nbsp;|&nbsp; " +
                               "<b>Evidence:</b> " + cat;
                $("#info_content").html(infoHtml).css("color", "#f8fafc");

                viewer.removeAllLabels();
                viewer.addLabel("Site " + resi, {{ fontSize: 14, fontColor: '#0f172a', backgroundColor: '#facc15', position: {{ x: a.x, y: a.y, z: a.z }} }});
            }}
            viewer.render();
        }}

        $(function() {{
            if (typeof window["$3Dmol"] === "undefined") {{
                $("#info_content").html("<b>3Dmol.js could not be loaded from the CDN. " +
                    "The structure cannot be displayed offline.</b>").css("color", "#ef4444");
                return;
            }}
            viewer = $3Dmol.createViewer($("#viewer"), {{ backgroundColor: '#0b0f19' }});
            viewer.addModel($("#pdb_data").text(), "pdb");
            baseStyle(0.8);
{style_rules}
            viewer.zoomTo();
            viewer.render();

            $(".site-btn").click(function() {{
                $(".site-btn").removeClass("active");
                $(this).addClass("active");
                zoomToSite($(this).data("resi"));
            }});
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
