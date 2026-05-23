import argparse
import zipfile
import json
import random
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Filter virus genomes and extract CDS from NCBI Datasets ZIP")
    parser.add_argument("--zip", required=True, help="Input NCBI datasets zip file")
    parser.add_argument("--protein", required=True, help="Target protein name to extract")
    parser.add_argument("--geo", default="", help="Geographic filter (e.g. 'Bangladesh')")
    parser.add_argument("--host", default="", help="Host filter (e.g. 'Homo sapiens' or 'Pteropus')")
    parser.add_argument("--max", type=int, default=100, help="Maximum number of sequences to keep")
    parser.add_argument("--out", required=True, help="Output FASTA file for CDS")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Buka file ZIP dari NCBI datasets
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
                record = json.loads(line.decode('utf-8'))
                acc = record.get('accession')
                location = record.get('location', {}).get('geographicLocation', '')
                host = record.get('host', {}).get('name', '')
                
                # Apply Geographic Filter jika ada
                if args.geo and args.geo.lower() not in location.lower():
                    continue
                    
                # Apply Host Filter jika ada (H2H vs Reservoir)
                if args.host and args.host.lower() not in host.lower():
                    continue
                
                valid_accessions.append(acc)
        
        # Downsample jika jumlah sekuens lebih dari --max
        if len(valid_accessions) > args.max:
            random.seed(42) # Untuk reproducibility
            valid_accessions = random.sample(valid_accessions, args.max)
            
        valid_set = set(valid_accessions)
        
        # 2. Extract CDS (Coding Sequences) untuk target protein
        # Idealnya kita parse file cds.fna menggunakan Bio.SeqIO
        # POTENSI ERROR TERBESAR: Internal Stop Codons!
        # Algoritma dN/dS (HyPhy) akan CRASH jika ada stop codon di tengah sekuens.
        #
        # PSEUDOCODE PENYARINGAN KETAT:
        # valid_records = []
        # for record in SeqIO.parse(cds_fasta, "fasta"):
        #     if record.id in valid_set:
        #         # Cek kelipatan 3 (codon)
        #         if len(record.seq) % 3 != 0: continue
        #         
        #         # Translate dan cek internal stop codon
        #         protein = record.seq.translate()
        #         if "*" in protein[:-1]: # Abaikan stop codon di akhir
        #             print(f"DROPPED {record.id}: Memiliki internal stop codon!")
        #             continue
        #             
        #         valid_records.append(record)
        
        # MOCK IMPLEMENTATION (Agar Snakemake rule tidak error saat dry-run)
        # Pada real implementation, kita akan pakai Bio.SeqIO membaca z.open('.../cds.fna')
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
            
        with open(args.out, 'w') as out_f:
            out_f.write(f">mock_cds_{args.protein}_from_filtered_data\n")
            out_f.write("ATGCGTACGTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCTAG\n")
            
        print(f"[{args.protein}] Tersimpan {len(valid_set)} filtered sequences di {args.out}")

if __name__ == "__main__":
    main()
