import argparse
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Merge multiple FASTA files and prefix sequence headers with virus names")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input FASTA files to merge")
    parser.add_argument("--virus-names", nargs="+", required=True, help="Virus names corresponding to each input file")
    parser.add_argument("--output", required=True, help="Output merged FASTA file")
    return parser.parse_args()

def main():
    args = parse_args()
    
    if len(args.inputs) != len(args.virus_names):
        print("Error: Jumlah file input harus sama dengan jumlah nama virus.", flush=True)
        sys.exit(1)
        
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        
    total_seqs = 0
    with open(args.output, 'w') as out_f:
        for filepath, virus_name in zip(args.inputs, args.virus_names):
            if not os.path.exists(filepath):
                print(f"Warning: File {filepath} tidak ditemukan. Dilewati.", flush=True)
                continue
                
            # Cek jika file kosong
            if os.path.getsize(filepath) == 0:
                print(f"Warning: File {filepath} kosong. Dilewati.", flush=True)
                continue
                
            print(f"Membaca {filepath} untuk virus {virus_name}...", flush=True)
            with open(filepath, 'r') as in_f:
                current_header = None
                current_seq = []
                
                for line in in_f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        if current_header:
                            # Tulis record sebelumnya
                            seq_str = "".join(current_seq)
                            # Bersihkan header dari spasi dan karakter ilegal jika ada
                            clean_header = current_header.replace(" ", "_").replace(",", "_").replace(".", "_").replace("-", "_").replace(":", "_").replace("|", "_")
                            out_f.write(f">{virus_name}_{clean_header}\n{seq_str}\n")
                            total_seqs += 1
                        current_header = line.lstrip('>')
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                if current_header:
                    seq_str = "".join(current_seq)
                    clean_header = current_header.replace(" ", "_").replace(",", "_").replace(".", "_").replace("-", "_").replace(":", "_").replace("|", "_")
                    out_f.write(f">{virus_name}_{clean_header}\n{seq_str}\n")
                    total_seqs += 1
                    
    print(f"Selesai menggabungkan. Total sekuens ditulis: {total_seqs} ke {args.output}", flush=True)

if __name__ == "__main__":
    main()
