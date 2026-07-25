#!/usr/bin/env python3
"""
Extract a target protein CDS from unannotated whole viral genomes.

GISAID FASTA carries no CDS annotation — headers are only
"virus_name|EPI_ISL_xxxxx|year" and each record is a whole genome. The NCBI
path in extract_cds.py relies on the annotation shipped in the datasets ZIP, so
it cannot consume these records.

This script locates the target gene by open reading frame search across all six
frames, then identifies which ORF is the gene of interest by k-mer similarity to
a reference protein. That is robust for a non-segmented virus with conserved
gene order, tolerates the indels that a strict positional lift-over would break
on, and needs nothing beyond the standard library.

Sequences are only emitted when they pass the same quality gates as the NCBI
path: length within a configurable fraction of the reference, a length that is a
multiple of three, and no internal stop codon.
"""

import argparse
import collections
import csv
import os
import re
import sys

CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def translate(dna):
    return "".join(CODON_TABLE.get(dna[i:i + 3].upper(), 'X')
                   for i in range(0, len(dna) - 2, 3))


def revcomp(dna):
    return dna.translate(COMPLEMENT)[::-1]


def parse_fasta(path):
    """Yield (header, sequence) pairs."""
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def find_orfs(genome, min_aa):
    """
    Yield (nucleotide_cds, protein, strand, start) for every ATG-initiated ORF
    ending in a stop codon, across all six frames.

    Both strands are searched because deposited genomes vary between genome and
    antigenome sense for negative-strand viruses such as Nipah.
    """
    for strand, seq in (("+", genome), ("-", revcomp(genome))):
        for frame in range(3):
            protein = translate(seq[frame:])
            for match in re.finditer(r"M[^*]*\*", protein):
                aa = match.group()[:-1]          # drop the stop
                if len(aa) < min_aa:
                    continue
                start = frame + match.start() * 3
                cds = seq[start:start + (len(aa) + 1) * 3]
                yield cds, aa, strand, start


def kmer_similarity(a, b, k=5):
    """
    Fraction of shared k-mers, normalised by the smaller sequence.

    Used only to decide which ORF corresponds to the target gene, so a cheap
    indel-tolerant measure is sufficient and avoids a Biopython dependency.
    """
    if len(a) < k or len(b) < k:
        return 0.0
    ka = collections.Counter(a[i:i + k] for i in range(len(a) - k + 1))
    kb = collections.Counter(b[i:i + k] for i in range(len(b) - k + 1))
    shared = sum((ka & kb).values())
    return shared / max(1, min(sum(ka.values()), sum(kb.values())))


def load_metadata(path):
    """Map accession -> row dict, keyed on the GISAID 'Accession ID' column."""
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {(r.get("Accession ID") or "").strip(): r for r in reader}


def accession_from_header(header):
    """GISAID headers look like 'virus_name|EPI_ISL_20339844|2004'."""
    for part in header.split("|"):
        part = part.strip()
        if part.startswith("EPI_ISL"):
            return part
    return header.split("|")[0].strip()


def matches(value, pattern):
    """Case-insensitive regex-alternation match, empty pattern matches all."""
    if not pattern:
        return True
    return re.search(pattern, value or "", re.IGNORECASE) is not None


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genomes", required=True, help="Unannotated genome FASTA")
    p.add_argument("--reference-protein", required=True, dest="ref_protein",
                   help="FASTA holding the reference protein for this gene")
    p.add_argument("--metadata", default="", help="GISAID metadata TSV")
    p.add_argument("--out-nuc", required=True, dest="out_nuc")
    p.add_argument("--out-prot", required=True, dest="out_prot")
    p.add_argument("--host-filter", default="", dest="host_filter")
    p.add_argument("--geo-filter", default="", dest="geo_filter")
    p.add_argument("--min-identity", type=float, default=0.35, dest="min_identity",
                   help="Minimum k-mer similarity to the reference protein")
    p.add_argument("--min-length-fraction", type=float, default=0.80,
                   dest="min_len_frac",
                   help="Minimum CDS length as a fraction of the reference")
    p.add_argument("--max-length-fraction", type=float, default=1.20,
                   dest="max_len_frac")
    p.add_argument("--max-ambiguous-fraction", type=float, default=0.02,
                   dest="max_ambig",
                   help="Reject a CDS whose translation is more than this "
                        "fraction 'X', i.e. built from ambiguous codons")
    return p.parse_args()


def main():
    args = parse_args()

    refs = list(parse_fasta(args.ref_protein))
    if not refs:
        sys.exit(f"ERROR: no sequence in {args.ref_protein}")
    ref_aa = refs[0][1].replace("*", "").strip()
    lo = int(len(ref_aa) * args.min_len_frac)
    hi = int(len(ref_aa) * args.max_len_frac)
    print(f"Reference protein: {len(ref_aa)} aa; accepting ORFs of {lo}-{hi} aa")

    meta = load_metadata(args.metadata)
    counts = collections.Counter()
    kept = []

    for header, genome in parse_fasta(args.genomes):
        counts["genomes_read"] += 1
        acc = accession_from_header(header)
        row = meta.get(acc, {})

        if not matches(row.get("Host", ""), args.host_filter):
            counts["dropped_host"] += 1
            continue
        if not matches(row.get("Location", ""), args.geo_filter):
            counts["dropped_geo"] += 1
            continue

        best = None
        for cds, aa, strand, start in find_orfs(genome, lo):
            if not (lo <= len(aa) <= hi):
                continue
            score = kmer_similarity(aa, ref_aa)
            if best is None or score > best[0]:
                best = (score, cds, aa, strand, start)

        if best is None:
            counts["no_orf_in_range"] += 1
            continue
        score, cds, aa, strand, start = best
        if score < args.min_identity:
            counts["below_identity"] += 1
            continue
        if len(cds) % 3 != 0:
            counts["not_multiple_of_three"] += 1
            continue
        if "*" in translate(cds)[:-1]:
            counts["internal_stop"] += 1
            continue
        # Ambiguous nucleotides translate to 'X' and would otherwise pass every
        # check above while contributing no usable codon to the alignment.
        if aa.count("X") > args.max_ambig * len(aa):
            counts["too_ambiguous"] += 1
            continue

        counts["kept"] += 1
        kept.append((f"{acc}:{start}-{start + len(cds)}", cds, aa))

    os.makedirs(os.path.dirname(args.out_nuc) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_prot) or ".", exist_ok=True)
    with open(args.out_nuc, "w") as fn, open(args.out_prot, "w") as fp:
        for name, cds, aa in kept:
            fn.write(f">{name}\n{cds}\n")
            fp.write(f">{name}\n{aa}\n")

    print("Summary:")
    for key in ("genomes_read", "dropped_host", "dropped_geo", "no_orf_in_range",
                "below_identity", "not_multiple_of_three", "internal_stop",
                "too_ambiguous", "kept"):
        print(f"  {key:24s} {counts[key]}")
    if counts["kept"] == 0:
        print("WARNING: no sequences passed. Outputs are empty.", file=sys.stderr)


if __name__ == "__main__":
    main()
