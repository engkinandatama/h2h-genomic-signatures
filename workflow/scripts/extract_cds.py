import argparse
import zipfile
import json
import random
import os
import re
import sys

CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate_dna(dna_seq):
    dna_seq = dna_seq.upper().replace('-', '').replace('\n', '').replace('\r', '').strip()
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        amino_acid = CODON_TABLE.get(codon, 'X')
        protein.append(amino_acid)
    return "".join(protein)

def matches_any(value, pattern):
    """
    True when `value` contains any of the pipe-separated keywords in `pattern`.

    Config filters are alternations: geo_filter "Bangladesh|India" means either
    country. Testing the whole pattern as one substring can never match a field
    holding a single country, so every record is dropped and the dataset comes
    back empty. That is what emptied Nipah_NiVB_H2H the moment India was merged
    into the NiV-B clade; while the value was the single word "Bangladesh" the
    substring test happened to agree with the intent.
    """
    if not pattern:
        return True
    haystack = (value or "").lower()
    return any(k.strip().lower() in haystack
               for k in pattern.split("|") if k.strip())


def load_reference_protein(path):
    """Residues of the UniProt reference for this protein, or '' when absent."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return ""
    residues = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                residues.append(line.strip())
    return "".join(residues).replace("*", "")


def _kmers(seq, k=5):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def reference_gate(nuc_seqs, prot_seqs, reference, min_frac, max_frac,
                   min_identity, tag):
    """
    Indices to keep after gating candidates against the UniProt reference.

    Two independent checks, because neither alone is sufficient.

    Length, measured against the reference rather than against the candidates
    themselves. The previous gate took the 90th percentile of whatever came back,
    which is circular: when BV-BRC returned 135 short CDS for Puumala L, p90 was
    273 bp, the floor became 191 bp, every fragment passed, and a genuine 6471 bp
    L protein would have been rejected for exceeding the ceiling. A gate derived
    from the contaminated set blesses the contamination.

    Identity, as the fraction of a candidate's 5-mers found in the reference.
    Length cannot separate neighbouring paralogs: a Marburg nucleoprotein is
    695 aa against a 681 aa glycoprotein and reached the GP dataset at 102% of
    the reference length. Across all 1033 sequences in this study the two
    populations do not overlap -- every on-target record scores at least 0.459
    and every off-target one 0.00 -- so the default cut sits between them.
    """
    n = len(nuc_seqs)
    notes = []
    if n == 0:
        return [], notes

    if reference:
        expected = (len(reference) + 1) * 3      # +1 codon for the stop
        floor, ceiling = int(min_frac * expected), int(max_frac * expected)
        basis = f"UniProt reference {len(reference)} aa = {expected} bp"
        ref_kmers = _kmers(reference)
    elif n >= 5:
        lengths = sorted(len(s) for s in nuc_seqs)
        p90 = lengths[int(0.9 * (n - 1))]
        floor, ceiling = int(min_frac * p90), int(max_frac * p90)
        basis = f"90th percentile of candidates {p90} bp"
        ref_kmers = set()
        notes.append(f"[{tag}] WARNING: no reference protein supplied; the length "
                     f"gate falls back to the candidate distribution and cannot "
                     f"detect a dataset that is contaminated throughout.")
    else:
        return list(range(n)), notes

    keep, short, long_, off = [], 0, 0, []
    for i, nuc in enumerate(nuc_seqs):
        if len(nuc) < floor:
            short += 1
            continue
        if len(nuc) > ceiling:
            long_ += 1
            continue
        if ref_kmers and min_identity > 0:
            prot = prot_seqs[i] if i < len(prot_seqs) else ""
            q = _kmers((prot or "").replace("*", ""))
            ident = len(q & ref_kmers) / len(q) if q else 1.0
            if ident < min_identity:
                off.append((i, round(ident, 3), len(prot or "")))
                continue
        keep.append(i)

    notes.append(f"[{tag}] Reference gate ({basis}): keep {floor}-{ceiling} bp and "
                 f"k-mer identity >= {min_identity}.")
    if short or long_:
        notes.append(f"[{tag}]   dropped {short} too short, {long_} too long.")
    if off:
        notes.append(f"[{tag}]   dropped {len(off)} off-target (wrong gene): "
                     f"{[(i, s, f'{l}aa') for i, s, l in off[:5]]}")
    notes.append(f"[{tag}]   kept {len(keep)} of {n}.")
    return keep, notes

def parse_args():
    parser = argparse.ArgumentParser(description="Filter virus genomes and extract CDS from NCBI Datasets ZIP")
    parser.add_argument("--zip", required=True, help="Input NCBI datasets zip file")
    parser.add_argument("--protein", required=True, help="Target protein name to extract")
    parser.add_argument("--geo", default="", help="Geographic filter (e.g. 'Bangladesh')")
    parser.add_argument("--host", default="", help="Host filter (e.g. 'Homo sapiens' or 'Pteropus')")
    parser.add_argument("--max", type=int, default=100, help="Maximum number of sequences to keep")
    parser.add_argument("--min-len", type=int, default=1500, help="Minimum sequence length (bp) to avoid partial fragments")
    parser.add_argument("--min-len-fraction", type=float, default=0.70,
                        dest="min_len_fraction",
                        help="Minimum CDS length as a fraction of the 90th-percentile "
                             "candidate length for this protein; 0 disables the gate")
    parser.add_argument("--reference", default="",
                        dest="reference",
                        help="UniProt reference protein FASTA for this virus and "
                             "protein. Defines the expected CDS length and the "
                             "k-mer identity a candidate must reach.")
    parser.add_argument("--min-reference-identity", type=float, default=0.30,
                        dest="min_reference_identity",
                        help="Minimum fraction of a candidate's 5-mers that must "
                             "occur in the reference protein. Default 0.30.")
    parser.add_argument("--max-len-fraction", type=float, default=1.30,
                        dest="max_len_fraction",
                        help="Maximum CDS length as a fraction of the same reference; "
                             "rejects a longer paralog that clears the floor from above")
    parser.add_argument("--out-nuc", required=True, help="Output FASTA file for CDS nucleotides")
    parser.add_argument("--out-prot", required=True, help="Output FASTA file for protein translation")
    return parser.parse_args()

def parse_fasta_header(header):
    main_id = header.split()[0].lstrip('>')
    attrs = {}
    matches = re.findall(r'\[([a-zA-Z0-9_]+)=([^\]]+)\]', header)
    for key, val in matches:
        attrs[key] = val.strip()
    
    # Extract protein name outside brackets
    desc_match = re.search(r'^[^\s]+\s+([^\[]+)', header)
    if desc_match:
        attrs['desc_protein'] = desc_match.group(1).strip()
        
    return main_id, attrs

def extract_genomic_accession(main_id):
    # main_id format: e.g. "AY228237.1:12-1298" atau "lcl|AY228237.1_cds_1"
    base = main_id.split(':')[0]
    
    # Hapus prefix lcl| jika ada dari NCBI datasets
    if base.startswith('lcl|'):
        base = base[4:]
        
    # Hapus suffix _cds_xxx atau _prot_xxx jika ada
    base = re.sub(r'_(cds|prot)_.*$', '', base)
    
    return base

def _has_keyword(text, keyword):
    """
    Substring match bounded to whole tokens.

    A bare `in` test made "l protein" match "nonstructural protein short": the
    keyword sits inside "nonstructura|l protein". That is how 135 copies of the
    273 bp NSs gene entered the Puumala L dataset and outnumbered the real
    6471 bp polymerase 100 to 1, which in turn made the length gate calibrate
    itself on the contamination. The same trap is open for every short keyword
    here -- gn, gc, gp, l.
    """
    if not text or not keyword:
        return False
    return re.search(r"(?<![a-z0-9])" + re.escape(keyword) + r"(?![a-z0-9])",
                     text) is not None


def _any_keyword(text, keywords):
    return any(_has_keyword(text, k) for k in keywords)


def protein_matches(target, gene_val, protein_val):
    target = target.lower().strip()
    gene_val = gene_val.lower().strip() if gene_val else ""
    protein_val = protein_val.lower().strip() if protein_val else ""
    
    # Annotations that carry the target's keywords but name a different gene.
    # Without these, "polymerase cofactor" (VP35) is collected as L, and the
    # "small secreted glycoprotein" (sGP) and "super small secreted glycoprotein"
    # (ssGP) products of filovirus transcriptional editing are collected as GP.
    # sGP and ssGP share only the N-terminal ~295 codons with GP before shifting
    # frame, so mixing them destroys the codon alignment.
    exclusions = {
        "gngc":      ["polymerase", "nucleoprotein", "nucleocapsid", "cofactor"],
        "l_protein": ["nonstructural", "nss", "cofactor", "complex protein", "vp35", "vp30", "vp24",
                      "vp40", "nucleoprotein", "matrix", "phosphoprotein"],
        "l":         ["cofactor", "complex protein", "vp35", "vp30", "vp24",
                      "vp40", "nucleoprotein", "matrix", "phosphoprotein"],
        "g_protein": ["polymerase", "fusion", "nucleoprotein", "phosphoprotein"],
        "f_protein": ["polymerase", "nucleoprotein", "phosphoprotein"],
        "gp":        ["secreted", "sgp", "ssgp", "soluble", "delta peptide",
                      "polymerase", "nucleoprotein", "cofactor"],
    }
    for bad in exclusions.get(target, []):
        if bad in protein_val or bad in gene_val:
            return False

    if target == "gngc":
        # Hantavirus GnGc (Glycoprotein precursor / M segment glycoprotein)
        return (_any_keyword(protein_val, ["glycoprotein", "gpc", "gn", "gc"])
                or gene_val in ("m", "gpc", "gn", "gc"))
    elif target in ["l_protein", "l"]:
        # RNA-dependent RNA polymerase / L protein
        return (_any_keyword(protein_val, ["polymerase", "large protein",
                                           "l protein", "rdrp", "transcriptase"])
                or gene_val == "l" or _has_keyword(gene_val, "polymerase"))
    elif target == "g_protein":
        # Nipah G protein (Attachment glycoprotein)
        return (_any_keyword(protein_val, ["glycoprotein g", "attachment",
                                           "g protein", "g-protein",
                                           "receptor-binding"])
                or gene_val == "g" or protein_val == "glycoprotein")
    elif target == "f_protein":
        # Nipah F protein (Fusion glycoprotein)
        return (_any_keyword(protein_val, ["fusion", "f protein", "f-protein"])
                or gene_val == "f")
    elif target == "gp":
        # Filovirus GP (full-length envelope glycoprotein only)
        if gene_val == "l":
            return False
        return _any_keyword(protein_val, ["glycoprotein", "gp"]) or gene_val == "gp"
    else:
        return _has_keyword(protein_val, target) or _has_keyword(gene_val, target)

def write_empty_outputs(args):
    # Buat file kosong agar Snakemake tidak mengeluh "Missing output files"
    nuc_dir = os.path.dirname(args.out_nuc)
    prot_dir = os.path.dirname(args.out_prot)
    if nuc_dir:
        os.makedirs(nuc_dir, exist_ok=True)
    if prot_dir:
        os.makedirs(prot_dir, exist_ok=True)
    with open(args.out_nuc, 'w') as f: pass
    with open(args.out_prot, 'w') as f: pass

def main():
    args = parse_args()
    
    # Debug info: Cek ukuran file zip
    try:
        zip_size = os.path.getsize(args.zip)
        print(f"DEBUG: Membaca file ZIP: {args.zip} ({zip_size} bytes)")
    except Exception as e:
        print(f"DEBUG: Gagal membaca ukuran file ZIP: {e}")
        write_empty_outputs(args)
        return

    try:
        with zipfile.ZipFile(args.zip, 'r') as z:
            # 1. Parsing data_report.jsonl untuk filtering metadata
            metadata_file = [f for f in z.namelist() if f.endswith('data_report.jsonl')]
            if not metadata_file:
                print(f"Error: data_report.jsonl tidak ditemukan di dalam {args.zip}.")
                print("NCBI datasets mungkin mengembalikan hasil kosong (tidak ada sekuens yang cocok).")
                print(f"DEBUG: Isi file ZIP ({len(z.namelist())} file):")
                for name in z.namelist()[:30]:
                    print(f" - {name}")
                if len(z.namelist()) > 30:
                    print(" - ... (dan file lainnya)")
                write_empty_outputs(args)
                return
            
            valid_accessions = []
            with z.open(metadata_file[0]) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        record = json.loads(line.decode('utf-8'))
                    except Exception as e:
                        print(f"Warning: gagal parse JSON line: {e}")
                        continue
                    acc = record.get('accession')
                    
                    # Menggunakan OR {} untuk menghindari AttributeError jika objek bernilai None di JSON
                    loc_obj = record.get('location')
                    if not isinstance(loc_obj, dict):
                        loc_obj = {}
                    location = loc_obj.get('geographicLocation') or ''
                    
                    host_obj = record.get('host')
                    if not isinstance(host_obj, dict):
                        host_obj = {}
                    host_name = host_obj.get('organismName') or host_obj.get('name') or ''
                    
                    # Apply Geographic Filter jika ada
                    if not matches_any(location, args.geo):
                        continue

                    # Apply Host Filter jika ada (H2H vs Reservoir)
                    if not matches_any(host_name, args.host):
                        continue
                    
                    valid_accessions.append(acc)
            
            valid_set = set(valid_accessions)
            print(f"Ditemukan {len(valid_set)} accession genom yang cocok setelah filter Host='{args.host}' dan Geo='{args.geo}'")
            
            # 2. Cari file cds.fna di dalam ZIP
            cds_files = [f for f in z.namelist() if f.endswith('cds.fna')]
            if not cds_files:
                print(f"Error: cds.fna tidak ditemukan di dalam {args.zip}.")
                write_empty_outputs(args)
                return
                
            # Parse file cds.fna
            records_to_process = []
            current_header = None
            current_seq = []
            
            with z.open(cds_files[0], 'r') as f:
                for line in f:
                    line_str = line.decode('utf-8').strip()
                    if not line_str:
                        continue
                    if line_str.startswith('>'):
                        if current_header:
                            records_to_process.append((current_header, "".join(current_seq)))
                        current_header = line_str
                        current_seq = []
                    else:
                        current_seq.append(line_str)
                if current_header:
                    records_to_process.append((current_header, "".join(current_seq)))
                    
            # 3. Filter sekuens berdasarkan protein target dan kualitas
            passed_nucleotides = []
            passed_proteins = []
            
            dropped_by_length = 0
            dropped_by_codon = 0
            dropped_by_stop = 0
            dropped_by_host = 0
            dropped_by_protein = 0
            
            seen_ids = set()
            
            for header, seq in records_to_process:
                main_id, attrs = parse_fasta_header(header)
                genomic_acc = extract_genomic_accession(main_id)
                
                # Filter 0: Cek duplikasi main_id
                if main_id in seen_ids:
                    continue
                seen_ids.add(main_id)
                
                # Filter 1: Cek apakah genomic accession ada di set valid
                if genomic_acc not in valid_set:
                    dropped_by_host += 1
                    continue
                    
                # Filter 2: Pencocokan nama protein/gene
                gene_val = attrs.get('gene', '')
                protein_val = attrs.get('protein', '') or attrs.get('desc_protein', '')
                if not protein_matches(args.protein, gene_val, protein_val):
                    dropped_by_protein += 1
                    continue
                    
                # Filter 3: Filter panjang minimum sekuens (hindari parsial)
                if args.min_len and len(seq) < args.min_len:
                    dropped_by_length += 1
                    continue
                    
                # Filter 4: Kelipatan 3 (codon)
                if len(seq) % 3 != 0:
                    dropped_by_codon += 1
                    continue
                    
                # Filter 5: Cek internal stop codon
                translated_seq = translate_dna(seq)
                if '*' in translated_seq[:-1]: # Abaikan stop codon di ujung akhir sekuens
                    dropped_by_stop += 1
                    continue
                    
                # Sekuens lolos seleksi
                # Bersihkan sekuens dari spasi/new line jika ada
                clean_seq = seq.upper().replace('\n', '').replace('\r', '').strip()
                passed_nucleotides.append((main_id, clean_seq))
                passed_proteins.append((main_id, translated_seq))
                
            print(f"[{args.protein}] Filter Stats -> Host/Geo: -{dropped_by_host}, Protein: -{dropped_by_protein}, Length: -{dropped_by_length}, Codon: -{dropped_by_codon}, StopCodon: -{dropped_by_stop}")
            print(f"[{args.protein}] Berhasil meloloskan {len(passed_nucleotides)} sekuens CDS berkualitas tinggi.")
            
            # Jika kosong, tetap tulis file kosong agar pipeline tidak crash karena missing file
            if not passed_nucleotides:
                write_empty_outputs(args)
                return

            # Relative length gate. A flat --min-len cannot distinguish a 265 bp
            # surveillance fragment from a 2031 bp full-length CDS, so partial
            # sequences reach the alignment and reduce per-site occupancy to a few
            # percent, which is what happened to the Puumala L alignment. The
            # reference is the 90th percentile of candidate lengths rather than the
            # maximum, so one over-long mis-annotation cannot raise the bar for all.
            if args.min_len_fraction > 0:
                reference = load_reference_protein(args.reference)
                keep, notes = reference_gate(
                    [sq for _, sq in passed_nucleotides],
                    [sq for _, sq in passed_proteins],
                    reference, args.min_len_fraction, args.max_len_fraction,
                    args.min_reference_identity, args.protein)
                for line in notes:
                    print(line)
                passed_nucleotides = [passed_nucleotides[i] for i in keep]
                passed_proteins = [passed_proteins[i] for i in keep]
                if not passed_nucleotides:
                    print(f"[{args.protein}] WARNING: every sequence fell outside "
                          f"the relative length gate.", file=sys.stderr)
                    write_empty_outputs(args)
                    return

            # Deduplicate BEFORE applying the cap. A single genome contributes one
            # record per annotated product, and the same isolate is often deposited
            # more than once, so capping first spends the quota on duplicates that
            # the downstream merge then removes. In the 2026-05-24 run this reduced
            # the Ebola non-human group from 105 candidates to 8 sequences.
            seen_seq = set()
            dedup_nuc, dedup_prot = [], []
            for (nid, nseq), (pid, pseq) in zip(passed_nucleotides, passed_proteins):
                if nseq in seen_seq:
                    continue
                seen_seq.add(nseq)
                dedup_nuc.append((nid, nseq))
                dedup_prot.append((pid, pseq))
            n_dupes = len(passed_nucleotides) - len(dedup_nuc)
            if n_dupes:
                print(f"[{args.protein}] Removed {n_dupes} duplicate sequences before capping.")
            passed_nucleotides, passed_proteins = dedup_nuc, dedup_prot

            # Downsample jika jumlah melebihi --max
            if len(passed_nucleotides) > args.max:
                random.seed(42)
                indices = sorted(random.sample(range(len(passed_nucleotides)), args.max))
                passed_nucleotides = [passed_nucleotides[i] for i in indices]
                passed_proteins = [passed_proteins[i] for i in indices]
                print(f"[{args.protein}] Downsampled menjadi {args.max} sekuens sesuai batas maksimal.")
                
            # 4. Tulis hasil output
            nuc_dir = os.path.dirname(args.out_nuc)
            prot_dir = os.path.dirname(args.out_prot)
            if nuc_dir:
                os.makedirs(nuc_dir, exist_ok=True)
            if prot_dir:
                os.makedirs(prot_dir, exist_ok=True)
            
            # Tulis nukleotida
            with open(args.out_nuc, 'w') as out_n:
                for main_id, seq in passed_nucleotides:
                    out_n.write(f">{main_id}\n{seq}\n")
                    
            # Tulis protein
            with open(args.out_prot, 'w') as out_p:
                for main_id, seq in passed_proteins:
                    out_p.write(f">{main_id}\n{seq}\n")
                    
            print(f"[{args.protein}] Hasil nukleotida disimpan di: {args.out_nuc}")
            print(f"[{args.protein}] Hasil protein disimpan di: {args.out_prot}")

    except zipfile.BadZipFile:
        print(f"Error: {args.zip} bukan file ZIP yang valid atau korup.")
        print("Ini biasanya terjadi jika NCBI Datasets memblokir request karena rate limit atau koneksi terputus.")
        print("Pastikan NCBI_API_KEY diset dan coba batasi concurrency Snakemake.")
        import sys; sys.exit(1)

if __name__ == "__main__":
    main()
