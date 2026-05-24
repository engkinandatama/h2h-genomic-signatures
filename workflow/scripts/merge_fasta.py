#!/usr/bin/env python3
import argparse
import os

def parse_fasta(file_path):
    records = []
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return records
        
    current_header = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    records.append((current_header, "".join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
                
    if current_header:
        records.append((current_header, "".join(current_seq)))
        
    return records

def extract_accession(header):
    # Try to extract a clean accession number
    # E.g. ">AY228237.1" -> "AY228237.1"
    # ">lcl|AY228237.1_cds_1" -> "AY228237.1"
    header = header.split()[0].lstrip('>')
    if header.startswith('lcl|'):
        header = header[4:]
    if '_cds_' in header:
        header = header.split('_cds_')[0]
    return header

def merge_and_deduplicate(ncbi_file, bvbrc_file, output_file, max_seq):
    ncbi_records = parse_fasta(ncbi_file)
    bvbrc_records = parse_fasta(bvbrc_file)
    
    seen_accessions = set()
    seen_sequences = set()
    final_records = []
    
    # Process NCBI first (usually our primary trusted source)
    for header, seq in ncbi_records:
        acc = extract_accession(header)
        if acc not in seen_accessions and seq not in seen_sequences:
            final_records.append((header, seq))
            seen_accessions.add(acc)
            seen_sequences.add(seq)
            
    # Process BV-BRC
    bvbrc_added = 0
    for header, seq in bvbrc_records:
        acc = extract_accession(header)
        if acc not in seen_accessions and seq not in seen_sequences:
            final_records.append((header, seq))
            seen_accessions.add(acc)
            seen_sequences.add(seq)
            bvbrc_added += 1
            
    print(f"Merged {len(ncbi_records)} NCBI and {len(bvbrc_records)} BV-BRC records.")
    print(f"Added {bvbrc_added} unique records from BV-BRC.")
    print(f"Total unique records before sampling: {len(final_records)}")
    
    # Truncate to max_seq if needed (we shouldn't randomly downsample here unless required, 
    # but the inputs are already downsampled to max_seq individually. If the sum exceeds max_seq, we cap it).
    if len(final_records) > max_seq:
        # Keep all from NCBI (up to max) and fill the rest with BV-BRC
        final_records = final_records[:max_seq]
        print(f"Capped total records to {max_seq}.")
        
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        for header, seq in final_records:
            f.write(f"{header}\n{seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge and deduplicate FASTA files from NCBI and BV-BRC")
    parser.add_argument("--ncbi", required=True)
    parser.add_argument("--bvbrc", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max", type=int, default=100)
    args = parser.parse_args()
    
    merge_and_deduplicate(args.ncbi, args.bvbrc, args.out, args.max)
