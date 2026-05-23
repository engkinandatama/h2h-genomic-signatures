import argparse
import zipfile
import json
import random
import os
import re

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

def parse_args():
    parser = argparse.ArgumentParser(description="Filter virus genomes and extract CDS from NCBI Datasets ZIP")
    parser.add_argument("--zip", required=True, help="Input NCBI datasets zip file")
    parser.add_argument("--protein", required=True, help="Target protein name to extract")
    parser.add_argument("--geo", default="", help="Geographic filter (e.g. 'Bangladesh')")
    parser.add_argument("--host", default="", help="Host filter (e.g. 'Homo sapiens' or 'Pteropus')")
    parser.add_argument("--max", type=int, default=100, help="Maximum number of sequences to keep")
    parser.add_argument("--min-len", type=int, default=1500, help="Minimum sequence length (bp) to avoid partial fragments")
    parser.add_argument("--out-nuc", required=True, help="Output FASTA file for CDS nucleotides")
    parser.add_argument("--out-prot", required=True, help="Output FASTA file for protein translation")
    return parser.parse_args()

def parse_fasta_header(header):
    main_id = header.split()[0].lstrip('>')
    attrs = {}
    matches = re.findall(r'\[([a-zA-Z0-9_]+)=([^\]]+)\]', header)
    for key, val in matches:
        attrs[key] = val.strip()
    return main_id, attrs

def extract_genomic_accession(main_id):
    match = re.search(r'lcl\|([A-Z0-9_.]+)_cds', main_id)
    if match:
        return match.group(1)
    if '_cds' in main_id:
        return main_id.split('_cds')[0]
    return main_id

def protein_matches(target, gene_val, protein_val):
    target = target.lower().strip()
    gene_val = gene_val.lower().strip() if gene_val else ""
    protein_val = protein_val.lower().strip() if protein_val else ""
    
    if target == "gngc":
        # Hantavirus GnGc (Glycoprotein precursor / M segment glycoprotein)
        return any(x in protein_val for x in ["glycoprotein", "gpc", "gn", "gc"]) or any(x in gene_val for x in ["m", "gpc", "gn", "gc"])
    elif target in ["l_protein", "l"]:
        # RNA-dependent RNA polymerase / L protein
        return any(x in protein_val for x in ["polymerase", "large protein", "l protein", "rdrp", "transcriptase"]) or gene_val == "l" or "polymerase" in gene_val
    elif target == "g_protein":
        # Nipah G protein (Attachment glycoprotein)
        return any(x in protein_val for x in ["glycoprotein g", "attachment", "g protein", "g-protein"]) or gene_val == "g"
    elif target == "f_protein":
        # Nipah F protein (Fusion glycoprotein)
        return any(x in protein_val for x in ["fusion", "f protein", "f-protein"]) or gene_val == "f"
    elif target == "gp":
        # Ebola GP (Glycoprotein) - exclude L protein polymerases
        if "polymerase" in protein_val or gene_val == "l":
            return False
        return any(x in protein_val for x in ["glycoprotein", "gp"]) or gene_val == "gp"
    else:
        return target in protein_val or target in gene_val

def main():
    args = parse_args()
    
    with zipfile.ZipFile(args.zip, 'r') as z:
        # 1. Parsing data_report.jsonl untuk filtering metadata
        metadata_file = [f for f in z.namelist() if f.endswith('data_report.jsonl')]
        if not metadata_file:
            print(f"Error: data_report.jsonl tidak ditemukan di dalam {args.zip}.")
            print("NCBI datasets mungkin mengembalikan hasil kosong (tidak ada sekuens yang cocok).")
            import sys; sys.exit(1)
        
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
                host_name = host_obj.get('name') or ''
                
                # Apply Geographic Filter jika ada
                if args.geo and args.geo.lower() not in location.lower():
                    continue
                    
                # Apply Host Filter jika ada (H2H vs Reservoir)
                if args.host and args.host.lower() not in host_name.lower():
                    continue
                
                valid_accessions.append(acc)
        
        valid_set = set(valid_accessions)
        print(f"Ditemukan {len(valid_set)} accession genom yang cocok setelah filter Host='{args.host}' dan Geo='{args.geo}'")
        
        # 2. Cari file cds.fna di dalam ZIP
        cds_files = [f for f in z.namelist() if f.endswith('cds.fna')]
        if not cds_files:
            print(f"Error: cds.fna tidak ditemukan di dalam {args.zip}.")
            import sys; sys.exit(1)
            
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
        
        for header, seq in records_to_process:
            main_id, attrs = parse_fasta_header(header)
            genomic_acc = extract_genomic_accession(main_id)
            
            # Filter 1: Cek apakah genomic accession ada di set valid
            if genomic_acc not in valid_set:
                continue
                
            # Filter 2: Pencocokan nama protein/gene
            gene_val = attrs.get('gene', '')
            protein_val = attrs.get('protein', '')
            if not protein_matches(args.protein, gene_val, protein_val):
                continue
                
            # Filter 3: Filter panjang minimum sekuens (hindari parsial)
            if len(seq) < args.min_len:
                # print(f"DEBUG: Dropped {main_id} - Terlahu pendek ({len(seq)} bp)")
                continue
                
            # Filter 4: Kelipatan 3 (codon)
            if len(seq) % 3 != 0:
                # print(f"DEBUG: Dropped {main_id} - Panjang bukan kelipatan 3 ({len(seq)} bp)")
                continue
                
            # Filter 5: Cek internal stop codon
            translated_seq = translate_dna(seq)
            if '*' in translated_seq[:-1]: # Abaikan stop codon di ujung akhir sekuens
                # print(f"DEBUG: Dropped {main_id} - Memiliki internal stop codon")
                continue
                
            # Sekuens lolos seleksi
            # Bersihkan sekuens dari spasi/new line jika ada
            clean_seq = seq.upper().replace('\n', '').replace('\r', '').strip()
            passed_nucleotides.append((main_id, clean_seq))
            passed_proteins.append((main_id, translated_seq))
            
        print(f"[{args.protein}] Berhasil meloloskan {len(passed_nucleotides)} sekuens CDS berkualitas tinggi.")
        
        # Downsample jika jumlah melebihi --max
        if len(passed_nucleotides) > args.max:
            random.seed(42)
            indices = random.sample(range(len(passed_nucleotides)), args.max)
            passed_nucleotides = [passed_nucleotides[i] for i in indices]
            passed_proteins = [passed_proteins[i] for i in indices]
            print(f"[{args.protein}] Downsampled menjadi {args.max} sekuens sesuai batas maksimal.")
            
        # 4. Tulis hasil output
        os.makedirs(os.path.dirname(args.out_nuc), exist_ok=True)
        os.makedirs(os.path.dirname(args.out_prot), exist_ok=True)
        
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

if __name__ == "__main__":
    main()
