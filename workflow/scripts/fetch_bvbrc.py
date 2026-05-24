#!/usr/bin/env python3
import argparse
import json
import urllib.request
import urllib.parse
import time
import re
import random
import os

BVBRC_API = "https://www.bv-brc.org/api"

def translate_dna(dna_seq):
    # Standard codon translation table as fallback (not needed if protein fasta works, but kept just in case)
    codon_table = {
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
    dna_seq = dna_seq.upper().replace('-', '').replace('\n', '').replace('\r', '').strip()
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        protein.append(codon_table.get(codon, 'X'))
    return "".join(protein)

def protein_matches(target, product_val):
    target = target.lower().strip()
    product_val = product_val.lower().strip() if product_val else ""
    
    if target == "gngc":
        return any(x in product_val for x in ["glycoprotein", "gpc", "gn", "gc"]) or "m segment" in product_val
    elif target in ["l_protein", "l"]:
        return any(x in product_val for x in ["polymerase", "large protein", "l protein", "rdrp", "transcriptase"]) or product_val == "l"
    elif target == "g_protein":
        return any(x in product_val for x in ["glycoprotein g", "attachment", "g protein", "g-protein", "receptor-binding"]) or product_val == "glycoprotein"
    elif target == "f_protein":
        return any(x in product_val for x in ["fusion", "f protein", "f-protein"])
    elif target == "gp":
        if "polymerase" in product_val or product_val == "l":
            return False
        return any(x in product_val for x in ["glycoprotein", "gp"])
    else:
        return target in product_val

def bvbrc_get(endpoint, params_str, accept_header="application/json", retries=3):
    url = f"{BVBRC_API}/{endpoint}/?{params_str}"
    headers = {
        "Accept": accept_header,
        "Content-Type": "application/x-www-form-urlencoded",
        "User-Agent": "h2h-genomic-signatures/1.0",
    }
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=90) as resp:
                raw = resp.read()
                if accept_header == "application/json":
                    return json.loads(raw.decode("utf-8"))
                else:
                    return raw.decode("utf-8")
        except urllib.error.HTTPError as e:
            if attempt == retries - 1:
                return None
        except Exception as e:
            if attempt == retries - 1:
                return None
            time.sleep(2 ** attempt)
    return None

def fetch_genomes(taxon_id, geo_filter, host_filter):
    all_records = []
    offset = 0
    page_size = 500
    fields = "genome_id,host_name,host_common_name,geographic_location"
    
    while True:
        rql = f"eq(taxon_id,{taxon_id})&select({fields})&limit({page_size},{offset})"
        data = bvbrc_get("genome", rql)
        
        if data is None:
            break
            
        records = data if isinstance(data, list) else data.get("response", {}).get("docs", [])
        if not records:
            break
            
        all_records.extend(records)
        if len(records) < page_size:
            break
        offset += page_size
        time.sleep(0.5)
        
    valid_ids = []
    for r in all_records:
        host = (r.get("host_name") or r.get("host_common_name") or "").lower()
        geo = (r.get("geographic_location") or "").lower()
        
        if geo_filter and geo_filter.lower() not in geo:
            continue
            
        if host_filter:
            host_keywords = [k.strip().lower() for k in host_filter.split('|')]
            if not any(k in host for k in host_keywords):
                continue
                
        valid_ids.append(r.get("genome_id"))
        
    return valid_ids

def parse_fasta(fasta_text):
    records = {}
    if not fasta_text:
        return records
    current_id = None
    current_product = None
    current_seq = []
    
    for line in fasta_text.split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id:
                records[current_id] = ("".join(current_seq), current_product)
            header = line
            parts = header.split('|')
            if len(parts) >= 3:
                feature_id = parts[1].strip()
                product_part = parts[2]
                product = product_part.split('[')[0].strip()
                last_part = parts[-1].strip().rstrip(']')
                genome_id = last_part.strip()
                
                if feature_id == "undefined":
                    clean_product = re.sub(r'[^a-zA-Z0-9_]', '_', product)
                    feature_id = f"{genome_id}_{clean_product}"
                
                current_id = feature_id
                current_product = product
            else:
                current_id = None
                current_product = None
            current_seq = []
        else:
            if current_id:
                current_seq.append(line)
                
    if current_id:
        records[current_id] = ("".join(current_seq), current_product)
        
    return records

def fetch_features(genome_ids, protein_target, min_len):
    passed_nuc = []
    passed_prot = []
    
    # Process in batches of 100 genomes
    batch_size = 100
    for i in range(0, len(genome_ids), batch_size):
        batch = genome_ids[i:i+batch_size]
        # In RQL, genome_id list should NOT contain quotes, e.g. in(genome_id,(186538.100,186538.1000))
        genomes_str = ",".join(batch)
        
        rql = f"in(genome_id,({genomes_str}))&eq(feature_type,CDS)&limit(25000,0)"
        
        # Fetch DNA and Protein FASTA in parallel/sequence
        dna_fasta = bvbrc_get("genome_feature", rql, accept_header="application/dna+fasta")
        prot_fasta = bvbrc_get("genome_feature", rql, accept_header="application/protein+fasta")
        
        dna_records = parse_fasta(dna_fasta)
        prot_records = parse_fasta(prot_fasta)
        
        # Match by feature ID
        common_ids = set(dna_records.keys()) & set(prot_records.keys())
        
        for fid in common_ids:
            na_seq, product = dna_records[fid]
            aa_seq, _ = prot_records[fid]
            
            if not na_seq:
                continue
                
            if not protein_matches(protein_target, product):
                continue
                
            if min_len and len(na_seq) < min_len:
                continue
                
            if len(na_seq) % 3 != 0:
                continue
                
            # If server returned empty protein seq, translate locally
            if not aa_seq:
                aa_seq = translate_dna(na_seq)
                
            if '*' in aa_seq[:-1]:
                continue
                
            passed_nuc.append((fid, na_seq.upper()))
            passed_prot.append((fid, aa_seq.upper()))
            
        time.sleep(0.5)
            
    return passed_nuc, passed_prot

def write_empty_outputs(out_nuc, out_prot):
    for fpath in [out_nuc, out_prot]:
        os.makedirs(os.path.dirname(fpath), exist_ok=True)
        with open(fpath, 'w') as f:
            pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxon", required=True)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--geo", default="")
    parser.add_argument("--host", default="")
    parser.add_argument("--max", type=int, default=100)
    parser.add_argument("--min-len", type=int, default=1500)
    parser.add_argument("--out-nuc", required=True)
    parser.add_argument("--out-prot", required=True)
    args = parser.parse_args()
    
    valid_genome_ids = fetch_genomes(args.taxon, args.geo, args.host)
    print(f"[{args.protein}] Ditemukan {len(valid_genome_ids)} genome ID yang cocok dari BV-BRC (Host: {args.host}, Geo: {args.geo}).", flush=True)
    
    if not valid_genome_ids:
        write_empty_outputs(args.out_nuc, args.out_prot)
        return
        
    passed_nuc, passed_prot = fetch_features(valid_genome_ids, args.protein, args.min_len)
    print(f"[{args.protein}] Berhasil mengekstrak {len(passed_nuc)} CDS valid dari BV-BRC.", flush=True)
    
    if not passed_nuc:
        write_empty_outputs(args.out_nuc, args.out_prot)
        return
        
    if len(passed_nuc) > args.max:
        random.seed(42)
        indices = random.sample(range(len(passed_nuc)), args.max)
        passed_nuc = [passed_nuc[i] for i in indices]
        passed_prot = [passed_prot[i] for i in indices]
        print(f"[{args.protein}] Downsampled menjadi {args.max} sekuens.", flush=True)
        
    write_empty_outputs(args.out_nuc, args.out_prot) # ensure dirs exist
    
    with open(args.out_nuc, 'w') as f:
        for seq_id, seq in passed_nuc:
            f.write(f">{seq_id}\n{seq}\n")
            
    with open(args.out_prot, 'w') as f:
        for seq_id, seq in passed_prot:
            f.write(f">{seq_id}\n{seq}\n")

if __name__ == "__main__":
    main()
