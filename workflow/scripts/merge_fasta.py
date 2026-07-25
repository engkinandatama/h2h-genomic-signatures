#!/usr/bin/env python3
import argparse
import os
import random
import re

def parse_fasta(file_path):
    records = {}
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
                    seq_id = current_header.split()[0].lstrip('>')
                    records[seq_id] = (current_header, "".join(current_seq))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
                
    if current_header:
        seq_id = current_header.split()[0].lstrip('>')
        records[seq_id] = (current_header, "".join(current_seq))
        
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

def merge_and_deduplicate_paired(ncbi_nuc_file, ncbi_prot_file, bvbrc_nuc_file, bvbrc_prot_file, out_nuc_file, out_prot_file, max_seq,
                                 extra_nuc_file=None, extra_prot_file=None):
    # Parse paired records for NCBI
    ncbi_nuc_records = parse_fasta(ncbi_nuc_file)
    ncbi_prot_records = parse_fasta(ncbi_prot_file)

    # Parse paired records for BV-BRC
    bvbrc_nuc_records = parse_fasta(bvbrc_nuc_file)
    bvbrc_prot_records = parse_fasta(bvbrc_prot_file)

    # Optional third source (GISAID). Its sequences are deduplicated against the
    # public archives by exact nucleotide match, so records already present in
    # GenBank contribute nothing and only genuinely new isolates are added.
    extra_nuc_records = parse_fasta(extra_nuc_file) if extra_nuc_file and os.path.exists(extra_nuc_file) else {}
    extra_prot_records = parse_fasta(extra_prot_file) if extra_prot_file and os.path.exists(extra_prot_file) else {}

    # Pair by ID. Sort the intersections: iterating a set of strings follows an
    # order that depends on PYTHONHASHSEED, which would make the capped subset
    # below differ between runs on identical input.
    ncbi_ids = sorted(set(ncbi_nuc_records.keys()) & set(ncbi_prot_records.keys()))
    bvbrc_ids = sorted(set(bvbrc_nuc_records.keys()) & set(bvbrc_prot_records.keys()))
    
    seen_accessions = set()
    seen_nuc_sequences = set()
    final_ids = []
    final_nuc_records = {}
    final_prot_records = {}
    
    # Process NCBI first
    for seq_id in ncbi_ids:
        nuc_header, nuc_seq = ncbi_nuc_records[seq_id]
        prot_header, prot_seq = ncbi_prot_records[seq_id]
        
        acc = extract_accession(nuc_header)
        if acc not in seen_accessions and nuc_seq not in seen_nuc_sequences:
            final_ids.append(seq_id)
            final_nuc_records[seq_id] = (nuc_header, nuc_seq)
            final_prot_records[seq_id] = (prot_header, prot_seq)
            seen_accessions.add(acc)
            seen_nuc_sequences.add(nuc_seq)
            
    # Process BV-BRC
    bvbrc_added = 0
    for seq_id in bvbrc_ids:
        nuc_header, nuc_seq = bvbrc_nuc_records[seq_id]
        prot_header, prot_seq = bvbrc_prot_records[seq_id]
        
        acc = extract_accession(nuc_header)
        if acc not in seen_accessions and nuc_seq not in seen_nuc_sequences:
            final_ids.append(seq_id)
            final_nuc_records[seq_id] = (nuc_header, nuc_seq)
            final_prot_records[seq_id] = (prot_header, prot_seq)
            seen_accessions.add(acc)
            seen_nuc_sequences.add(nuc_seq)
            bvbrc_added += 1
            
    # Process the optional third source last, so the public archives take priority
    extra_ids = sorted(set(extra_nuc_records.keys()) & set(extra_prot_records.keys()))
    extra_added = 0
    for seq_id in extra_ids:
        nuc_header, nuc_seq = extra_nuc_records[seq_id]
        prot_header, prot_seq = extra_prot_records[seq_id]

        acc = extract_accession(nuc_header)
        if acc not in seen_accessions and nuc_seq not in seen_nuc_sequences:
            final_ids.append(seq_id)
            final_nuc_records[seq_id] = (nuc_header, nuc_seq)
            final_prot_records[seq_id] = (prot_header, prot_seq)
            seen_accessions.add(acc)
            seen_nuc_sequences.add(nuc_seq)
            extra_added += 1

    print(f"Merged {len(ncbi_ids)} NCBI and {len(bvbrc_ids)} BV-BRC paired records.")
    print(f"Added {bvbrc_added} unique records from BV-BRC.")
    if extra_ids:
        print(f"Third source supplied {len(extra_ids)} paired records; "
              f"{extra_added} were new after deduplication.")
    print(f"Total unique records before sampling: {len(final_ids)}")
    
    # Cap to max_seq if needed. Take a seeded random sample of the deduplicated
    # set rather than the first N: the leading records are ordered by accession,
    # so slicing would bias the sample toward whichever accessions sort first.
    if len(final_ids) > max_seq:
        rng = random.Random(42)
        final_ids = sorted(rng.sample(final_ids, max_seq))
        print(f"Capped total records to {max_seq} (seeded random sample of the "
              f"deduplicated set).")
        
    os.makedirs(os.path.dirname(out_nuc_file), exist_ok=True)
    os.makedirs(os.path.dirname(out_prot_file), exist_ok=True)
    
    with open(out_nuc_file, 'w') as fn, open(out_prot_file, 'w') as fp:
        for seq_id in final_ids:
            nuc_header, nuc_seq = final_nuc_records[seq_id]
            prot_header, prot_seq = final_prot_records[seq_id]
            fn.write(f"{nuc_header}\n{nuc_seq}\n")
            fp.write(f"{prot_header}\n{prot_seq}\n")
            
    print(f"Successfully wrote paired outputs. Total: {len(final_ids)} sekuens.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge and deduplicate FASTA files from NCBI and BV-BRC while preserving matching nucleotide/protein sequences")
    parser.add_argument("--ncbi-nuc", required=True, help="NCBI nucleotide FASTA path")
    parser.add_argument("--ncbi-prot", required=True, help="NCBI protein FASTA path")
    parser.add_argument("--bvbrc-nuc", required=True, help="BV-BRC nucleotide FASTA path")
    parser.add_argument("--bvbrc-prot", required=True, help="BV-BRC protein FASTA path")
    parser.add_argument("--out-nuc", required=True, help="Output nucleotide FASTA path")
    parser.add_argument("--out-prot", required=True, help="Output protein FASTA path")
    parser.add_argument("--max", type=int, default=100, help="Max sequences to output")
    parser.add_argument("--extra-nuc", default="",
                        help="Optional third-source nucleotide FASTA (e.g. GISAID)")
    parser.add_argument("--extra-prot", default="",
                        help="Optional third-source protein FASTA (e.g. GISAID)")
    args = parser.parse_args()

    merge_and_deduplicate_paired(
        args.ncbi_nuc, args.ncbi_prot,
        args.bvbrc_nuc, args.bvbrc_prot,
        args.out_nuc, args.out_prot,
        args.max,
        args.extra_nuc, args.extra_prot
    )
