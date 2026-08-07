"""
Build the manuscript tables from the finished pipeline outputs.

Emits a JSON file of table bodies that revise_manuscript_docx.py substitutes
into the submitted document, so no number in a table is ever transcribed by
hand.

Usage:
    python scripts/build_revision_tables.py <results_dir> <branch_dir> <out.json>

`branch_dir` holds the regenerated *_branch_selection.tsv files, which carry the
FAILED / NOT_APPLICABLE distinction for RELAX.
"""
import csv
import json
import os
import sys

import yaml

FAM_ORDER = ["Paramyxoviridae", "Hantaviridae", "Filoviridae"]
SHORT = {"L_protein": "L", "G_protein": "G", "F_protein": "F", "GnGc": "GnGc", "GP": "GP"}


def truthy(row, key):
    return str(row.get(key, "")).strip().lower() == "true"


def num(row, key, fmt="{:.4g}"):
    v = row.get(key)
    try:
        return fmt.format(float(v))
    except (TypeError, ValueError):
        return "n.d."


def main(results, branch_dir, out_path):
    cfg = yaml.safe_load(open("config/config.yaml"))
    groups_cfg = cfg["virus_groups"]

    def klass(g):
        cats = {cfg["viruses"][v]["category"] for v in groups_cfg[g]["viruses"]}
        return "H2H" if "H2H" in cats else "Spillover-only"

    order = sorted(groups_cfg,
                   key=lambda g: (FAM_ORDER.index(groups_cfg[g]["family"]),
                                  klass(g) != "H2H", g))

    def sites(g, p):
        return list(csv.DictReader(
            open(f"{results}/04_selection/{g}/{p}_all_sites.tsv"), delimiter="\t"))

    def relax(g, p):
        f = f"{branch_dir}/{g}__{p}.tsv"
        if not os.path.exists(f):
            return "n.d."
        rows = list(csv.DictReader(open(f), delimiter="\t"))
        if not rows:
            return "n.d."
        k = rows[0].get("relax_k", "NA")
        if k in ("NOT_APPLICABLE", "FAILED"):
            return "not applicable" if k == "NOT_APPLICABLE" else "test failed"
        try:
            return f"{float(k):.3f}"
        except (TypeError, ValueError):
            return "n.d."

    # ---------------------------------------------------------------- Table 1
    t1 = [["Virus group", "Protein", "Consensus sites", "Episodic-only sites",
           "RELAX K", "Purifying (%)"]]
    for g in order:
        label = f"{g.replace('_virus', '').replace('_', ' ')} ({klass(g)})"
        for j, p in enumerate(groups_cfg[g]["proteins"]):
            rs = sites(g, p)
            neg = sum(1 for r in rs if truthy(r, "consensus_neg"))
            t1.append([label if j == 0 else "",
                       SHORT.get(p, p),
                       str(sum(1 for r in rs if truthy(r, "consensus_pos"))),
                       str(sum(1 for r in rs if truthy(r, "episodic_only"))),
                       relax(g, p),
                       f"{100 * neg / len(rs):.1f}"])

    # ---------------------------------------------------------------- Table 2
    # Adaptive sites in the NiV-B clade, the lineage the manuscript follows.
    t2 = [["Protein", "Site", "Ref", "H2H variants", "MEME p", "MEME q",
           "FEL p", "FUBAR pp", "PRIME omnibus p"]]

    def read_fasta(path):
        seqs, head = {}, None
        for line in open(path):
            line = line.rstrip()
            if line.startswith(">"):
                head = line[1:]
                seqs[head] = []
            elif head:
                seqs[head].append(line)
        return {k: "".join(v) for k, v in seqs.items()}

    CODON = {}
    for b1 in "TCAG":
        for b2 in "TCAG":
            for b3 in "TCAG":
                CODON[b1 + b2 + b3] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[
                    ("TCAG".index(b1) * 16 + "TCAG".index(b2) * 4 + "TCAG".index(b3))]

    for p in groups_cfg["Nipah_NiVB"]["proteins"]:
        rs = sites("Nipah_NiVB", p)
        aln = read_fasta(f"{results}/02_aligned/Nipah_NiVB/{p}_codon_aligned.fasta")
        fg = {k: v for k, v in aln.items() if k.startswith("Nipah_NiVB_H2H")}
        ref = {k: v for k, v in aln.items() if k.startswith("Nipah_Reservoir")}
        for r in rs:
            if not (truthy(r, "consensus_pos") or truthy(r, "episodic_only")):
                continue
            i = (int(r["site"]) - 1) * 3
            def aa_counts(d):
                c = {}
                for s in d.values():
                    cod = s[i:i + 3].upper()
                    a = CODON.get(cod)
                    if a and a != "*":
                        c[a] = c.get(a, 0) + 1
                return c
            rc, fc = aa_counts(ref), aa_counts(fg)
            ref_aa = max(rc, key=rc.get) if rc else "?"
            tot = sum(fc.values()) or 1
            var = "; ".join(f"{a} ({100*n/tot:.1f}%)"
                            for a, n in sorted(fc.items(), key=lambda kv: -kv[1])
                            if a != ref_aa)[:38] or "none"
            t2.append([SHORT.get(p, p), r["site"], ref_aa, var,
                       num(r, "meme_p"), num(r, "meme_q"), num(r, "fel_p"),
                       num(r, "fubar_pp"), num(r, "prime_overall_p")])

    # ---------------------------------------------------------------- Table 3
    # Replaces the Andes gap-analysis table, whose premise (50-70% reservoir
    # gaps) no longer holds: every alignment now has 100% median occupancy.
    t3 = [["Virus group", "Protein", "Site", "Role", "MEME q", "FEL q",
           "Nearest breakpoint", "GARD status"]]
    bp = {}
    for g in order:
        for p in groups_cfg[g]["proteins"]:
            f = f"{results}/04_selection/{g}/{p}_gard_summary.txt"
            row = list(csv.DictReader(open(f), delimiter="\t"))[0]
            pos = [] if row["breakpoint_positions"] in ("None", "NA") else \
                  [round(int(x) / 3) for x in row["breakpoint_positions"].split(",")]
            bp[(g, p)] = (row["status"], pos)
    ROLE = {"G_protein": "Entry helper support", "F_protein": "Entry helper",
            "GnGc": "Entry", "GP": "Entry", "L_protein": "Replication"}
    ROLE["G_protein"] = "Entry"
    for g in order:
        for p in groups_cfg[g]["proteins"]:
            st, pos = bp[(g, p)]
            for r in sites(g, p):
                if not truthy(r, "consensus_pos"):
                    continue
                s = int(r["site"])
                d = min((abs(b - s) for b in pos), default=None)
                t3.append([g.replace("_virus", "").replace("_", " "), SHORT.get(p, p),
                           r["site"], ROLE.get(p, "?"),
                           num(r, "meme_q"), num(r, "fel_q"),
                           "none detected" if d is None else str(d), st])

    # ---------------------------------------------------------------- Table 4
    t4 = [["Virus group", "Protein", "Consensus sites", "Sites", "RELAX K"]]
    for g in ("Ebola_Zaire", "Sudan_ebolavirus", "Marburg_virus"):
        for p in groups_cfg[g]["proteins"]:
            rs = sites(g, p)
            sel = [r["site"] for r in rs if truthy(r, "consensus_pos")]
            t4.append([g.replace("_virus", "").replace("_", " "), SHORT.get(p, p),
                       str(len(sel)),
                       ", ".join(sel) if sel else "none", relax(g, p)])

    tables = {"0": t1, "1": t2, "2": t3, "3": t4}
    json.dump(tables, open(out_path, "w"), indent=1)
    for k, v in tables.items():
        print(f"  Tabel #{k}: {len(v)} baris x {len(v[0])} kolom")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
