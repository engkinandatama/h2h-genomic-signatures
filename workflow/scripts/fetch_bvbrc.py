#!/usr/bin/env python3
import argparse
import sys
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


def protein_matches(target, product_val):
    target = target.lower().strip()
    product_val = product_val.lower().strip() if product_val else ""
    
    # Mirrors the exclusion list in extract_cds.py: BV-BRC annotates VP35 as
    # "polymerase complex protein" and the edited filovirus products as
    # "small secreted glycoprotein GP" / "super small secreted glycoprotein GP",
    # all of which match the naive keyword tests below.
    exclusions = {
        "gngc":      ["polymerase", "nucleoprotein", "nucleocapsid", "cofactor"],
        "l_protein": ["cofactor", "complex protein", "vp35", "vp30", "vp24",
                      "vp40", "nucleoprotein", "matrix", "phosphoprotein",
                      "nonstructural", "nss", "ns protein"],
        "l":         ["cofactor", "complex protein", "vp35", "vp30", "vp24",
                      "vp40", "nucleoprotein", "matrix", "phosphoprotein",
                      "nonstructural", "nss", "ns protein"],
        "g_protein": ["polymerase", "fusion", "nucleoprotein", "phosphoprotein"],
        "f_protein": ["polymerase", "nucleoprotein", "phosphoprotein"],
        "gp":        ["secreted", "sgp", "ssgp", "soluble", "delta peptide",
                      "polymerase", "nucleoprotein", "cofactor"],
    }
    for bad in exclusions.get(target, []):
        if bad in product_val:
            return False

    if target == "gngc":
        return (_any_keyword(product_val, ["glycoprotein", "gpc", "gn", "gc"])
                or _has_keyword(product_val, "m segment"))
    elif target in ["l_protein", "l"]:
        return (_any_keyword(product_val, ["polymerase", "large protein",
                                           "l protein", "rdrp", "transcriptase"])
                or product_val == "l")
    elif target == "g_protein":
        return (_any_keyword(product_val, ["glycoprotein g", "attachment",
                                           "g protein", "g-protein",
                                           "receptor-binding"])
                or product_val == "glycoprotein")
    elif target == "f_protein":
        return _any_keyword(product_val, ["fusion", "f protein", "f-protein"])
    elif target == "gp":
        if product_val == "l":
            return False
        return _any_keyword(product_val, ["glycoprotein", "gp"])
    else:
        return _has_keyword(product_val, target)

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
            timeout = 90 if accept_header == "application/json" else 300
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                raw = resp.read()
                if accept_header == "application/json":
                    return json.loads(raw.decode("utf-8"))
                else:
                    return raw.decode("utf-8")
        except urllib.error.HTTPError as e:
            print(f"HTTP Error on attempt {attempt+1} for URL {url}: {e.code} - {e.reason}", flush=True)
            if attempt == retries - 1:
                return None
            time.sleep(2 ** attempt + random.random())
        except Exception as e:
            print(f"Generic error on attempt {attempt+1} for URL {url}: {e}", flush=True)
            if attempt == retries - 1:
                return None
            time.sleep(2 ** attempt + random.random())
    return None

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


def fetch_genomes(taxon_id, geo_filter, host_filter):
    all_records = []
    offset = 0
    page_size = 500
    fields = "genome_id,host_name,host_common_name,geographic_location"
    
    while True:
        rql = f"eq(taxon_id,{taxon_id})&select({fields})&limit({page_size},{offset})"
        data = bvbrc_get("genome", rql)
        
        if data is None:
            print("Error: Gagal menghubungi API BV-BRC saat mengambil metadata genome. Menghentikan pipeline.", flush=True)
            import sys
            sys.exit(1)
            
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
        
        if not matches_any(geo, geo_filter):
            continue

        if not matches_any(host, host_filter):
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
    failed_batches = []
    passed_nuc = []
    passed_prot = []
    
    # Process in batches of genomes. Two BV-BRC behaviours constrain the size:
    #   * without limit() the API silently truncates to 25 records, so a batch of
    #     100 genomes returned at most 25 CDS and the rest were never seen;
    #   * the protein+fasta endpoint hangs beyond roughly 10-15 genomes per
    #     request and returns TRUNCATED output rather than an error. Measured on
    #     taxon 186540: 10 genomes -> 2.8 s and 180 sequences, 20 genomes ->
    #     150 s and only 200 of 360, 39 genomes -> 150 s and 600 of 702.
    # Protein sequences are therefore no longer requested at all; a CDS
    # translation is exact, needs no network call, and cannot fall out of step
    # with its nucleotide record.
    batch_size = 25
    page_limit = 25000
    for i in range(0, len(genome_ids), batch_size):
        batch = genome_ids[i:i+batch_size]
        # In RQL, genome_id list should NOT contain quotes, e.g. in(genome_id,(186538.100,186538.1000))
        genomes_str = ",".join(batch)

        rql = f"in(genome_id,({genomes_str}))&eq(feature_type,CDS)&limit({page_limit})"
        
        # Fetch DNA and Protein FASTA in parallel/sequence
        dna_fasta = bvbrc_get("genome_feature", rql, accept_header="application/dna+fasta")

        # A batch that fails must not discard the dataset. Aborting on the first
        # failure lost all four filovirus datasets even though the batches before
        # it had already succeeded. Skip it, count it, and carry on.
        if dna_fasta is None:
            failed_batches.append(len(batch))
            print(f"WARNING: batch of {len(batch)} genomes could not be retrieved; "
                  f"continuing without it.", file=sys.stderr, flush=True)
            continue
        
        dna_records = parse_fasta(dna_fasta)

        # Guard against silent truncation: if a response comes back exactly at the
        # limit, records were almost certainly dropped and the batch must shrink.
        if len(dna_records) >= page_limit:
            print(f"ERROR: BV-BRC returned {len(dna_records)} CDS for a batch of "
                  f"{len(batch)} genomes, which is the requested limit. Results "
                  f"are truncated; reduce batch_size.", file=sys.stderr)
            sys.exit(1)
        
        for fid in sorted(dna_records):
            na_seq, product = dna_records[fid]
            aa_seq = translate_dna(na_seq)
            
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
            
    if failed_batches:
        print(f"WARNING: {len(failed_batches)} of "
              f"{(len(genome_ids) + batch_size - 1) // batch_size} batches were "
              f"skipped after repeated API failures, covering "
              f"{sum(failed_batches)} genomes. The BV-BRC contribution for this "
              f"dataset is incomplete.", file=sys.stderr, flush=True)
    return passed_nuc, passed_prot

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
    parser.add_argument("--genome-multiplier", type=int, default=5,
                        dest="genome_multiplier",
                        help="Query features for at most this many times --max "
                             "genomes. Ebola matched 3060 genomes for a 100-sequence "
                             "cap, so the pipeline issued ~30x the requests it "
                             "needed and timed out on the FASTA endpoint.")
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
    parser.add_argument("--out-nuc", required=True)
    parser.add_argument("--out-prot", required=True)
    args = parser.parse_args()
    
    valid_genome_ids = fetch_genomes(args.taxon, args.geo, args.host)
    print(f"[{args.protein}] Ditemukan {len(valid_genome_ids)} genome ID yang cocok dari BV-BRC (Host: {args.host}, Geo: {args.geo}).", flush=True)
    
    if not valid_genome_ids:
        write_empty_outputs(args.out_nuc, args.out_prot)
        return
        
    # Cap the genome list before requesting features. Every genome costs two FASTA
    # requests, and beyond a modest multiple of --max the extra genomes only feed
    # sequences that the cap discards later. The sample is seeded, so the subset is
    # reproducible across runs and machines.
    cap = max(args.max * args.genome_multiplier, args.max)
    if len(valid_genome_ids) > cap:
        random.seed(42)
        valid_genome_ids = sorted(random.sample(valid_genome_ids, cap))
        print(f"[{args.protein}] Sampled {cap} of the matching genomes for feature "
              f"retrieval (seeded); the cap is --max x --genome-multiplier.",
              flush=True)

    passed_nuc, passed_prot = fetch_features(valid_genome_ids, args.protein, args.min_len)
    print(f"[{args.protein}] Berhasil mengekstrak {len(passed_nuc)} CDS valid dari BV-BRC.", flush=True)
    
    if args.min_len_fraction > 0:
        reference = load_reference_protein(args.reference)
        keep, notes = reference_gate(
            [sq for _, sq in passed_nuc], [sq for _, sq in passed_prot],
            reference, args.min_len_fraction, args.max_len_fraction,
            args.min_reference_identity, args.protein)
        for line in notes:
            print(line, flush=True)
        passed_nuc = [passed_nuc[i] for i in keep]
        passed_prot = [passed_prot[i] for i in keep]

    if not passed_nuc:
        write_empty_outputs(args.out_nuc, args.out_prot)
        return
        
    if len(passed_nuc) > args.max:
        random.seed(42)
        indices = sorted(random.sample(range(len(passed_nuc)), args.max))
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
