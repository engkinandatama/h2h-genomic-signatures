"""
Mark up the submitted manuscript with the revised text and artwork.

Non-destructive by design: an original paragraph is struck through in grey and
its replacement inserted immediately below in blue. Nothing is deleted, so every
Zotero citation field survives and the author can see exactly what changed and
re-attach references where a sentence moved.

Artwork is swapped in place and the display box is rescaled to the new file's
aspect ratio, so a replacement figure is never stretched.

Usage:
    python scripts/revise_manuscript_docx.py <in.docx> <out.docx>
"""
import re
import struct
import sys
import zipfile

BLUE = "1F4E79"
GREY = "808080"
FIG = ".dev/manuscript/revisi_figures"

MEDIA = {
    "word/media/image1.png": f"{FIG}/Figure1_H2H_vs_spillover_rate.png",
    "word/media/image2.png": f"{FIG}/Figure2_PRIME_L210.png",
    "word/media/image4.png": f"{FIG}/Figure2_positional_bias.png",
}

REPLACEMENTS = {
35: [
 "The transition of zoonotic viruses from animal reservoirs to sustained human-to-human "
 "(H2H) transmission is a critical evolutionary bottleneck that determines pandemic "
 "potential. Whether that transition leaves a detectable genomic signature is an open "
 "question with direct consequences for genomic surveillance, and it has rarely been "
 "tested against a matched comparator. We addressed it in three negative-sense "
 "single-stranded RNA virus families using a paired design in which every family "
 "contributes both an H2H-capable lineage and a spillover-only relative sharing the same "
 "reservoir contrast: Andes virus against Sin Nombre and Puumala viruses (Hantaviridae), "
 "the Nipah virus NiV-B clade against the Malaysian outbreak lineage (Paramyxoviridae), "
 "and Ebola Zaire and Sudan ebolaviruses against Marburg virus (Filoviridae). Selective "
 "pressure was characterised with Contrast-FEL, MEME, FEL, FUBAR, SLAC, PRIME and RELAX, "
 "with Benjamini-Hochberg correction applied within each alignment, and all codon "
 "alignments were screened for recombination breakpoints with GARD. Across 25,467 codons "
 "in 18 alignments, 26 sites met the consensus criterion for positive selection and one "
 "survived false-discovery correction. H2H-capable lineages did not carry more adaptive "
 "sites than their spillover-only relatives. The rate was lower, not higher: 0.64 against "
 "1.38 sites per 1000 codons (Fisher exact odds ratio 0.46, p = 0.077), and within "
 "Paramyxoviridae the two rates were indistinguishable, 1.179 in both lineages. No bin of "
 "the 100-bin normalised coordinate framework was occupied by more than two of the eight "
 "lineages in any functional role. One positional pattern did hold across all three "
 "families, and it concerns protein architecture rather than transmission mode: adaptive "
 "sites in entry glycoproteins are displaced toward the N-terminus (mean relative position "
 "0.298 against an expected 0.500; q = 0.036 after correction across functional roles), "
 "whereas polymerase sites are not (0.610; q = 0.930). Recombination breakpoints were "
 "detected in four alignments and no differential site fell within ten codons of one. "
 "These results argue that H2H transmission capability, as currently sampled in public "
 "sequence archives, does not leave a selection signature that distinguishes it from "
 "spillover-only transmission in the same viral family. Genomic surveillance strategies "
 "that assume such a signature exists should be reconsidered, or should target far larger "
 "sequence sets than are presently available."],

56: [
 "Maximum-likelihood phylogenetic trees were reconstructed using IQ-TREE v2 with automatic "
 "model selection via ModelFinder and 1000 ultrafast bootstrap replicates (-m TEST -B "
 "1000). Branches with bootstrap support below 70 were collapsed into polytomies. All "
 "codon alignments were screened for recombination breakpoints using GARD implemented in "
 "HyPhy. Screening completed on 15 of the 18 alignments; two terminated on an internal "
 "numerical error after establishing part of the search and are reported as partial, and "
 "one was not screened because of its computational cost. Detected breakpoints are "
 "reported alongside the selection results, and the proximity of each differential site to "
 "the nearest breakpoint is given, rather than sites being removed from the analysis."],

70: [
 "Fisher's exact test was used to assess enrichment of differential selection sites within "
 "predefined protein domains (Gn domain, Gc domain, receptor-binding domain, fusion "
 "peptide, N-terminal, catalytic, and C-terminal regions), with Benjamini-Hochberg false "
 "discovery rate correction across the domains tested."],

74: [
 "All statistical analyses used Python 3.10 with SciPy, NumPy, and Pandas libraries. "
 "Selection rate comparisons between groups used Fisher's exact test on pooled site counts "
 "and the Mann-Whitney U test on per-lineage rates. Benjamini-Hochberg false discovery "
 "rate correction was applied to MEME, FEL, SLAC and PRIME within each alignment, treating "
 "each codon as one test; Contrast-FEL q values are those reported by HyPhy. A site was "
 "called differentially selected only if it also survived HyPhy's permutation test over "
 "branch assignments, which guards against calls driven by a single substitution against a "
 "background rate estimated near zero. Permutation tests reported for functional roles "
 "were corrected across the roles examined."],

304: [
 "In the glycoprotein GP: no site met the consensus criterion for positive selection after "
 "false-discovery correction, and one site was significant under MEME alone without "
 "support from the pervasive models. RELAX was not applied to Ebola Zaire: no natural "
 "reservoir sequences are available for this virus, so the tree has no foreground and "
 "background partition and the branch-aware tests are not applicable. The value previously "
 "reported here was estimated against sequences from experimentally infected macaques, "
 "which do not constitute a natural reservoir."],

306: [
 "In the polymerase L-protein: three sites met the consensus criterion for positive "
 "selection. As for the glycoprotein, RELAX was not applied, because Ebola Zaire has no "
 "natural reservoir comparator. Marburg virus supplies the reservoir contrast for "
 "Filoviridae; across the six lineages that do support a branch contrast, relaxation "
 "parameters did not differ between entry and polymerase proteins (paired Wilcoxon, "
 "p = 1.0; mean K 0.83 against 0.80)."],

334: [
 "No bin was shared by more than two lineages in any functional role. Under permutation, "
 "the best observed co-occurrence was not distinguishable from chance for entry proteins "
 "(2 lineages, p = 0.267), entry helper proteins (1 lineage, p = 1.0) or polymerases "
 "(2 lineages, p = 0.552); all q values exceeded 0.80. This test is underpowered by "
 "construction: with one to nine differential sites per lineage distributed over 100 bins, "
 "two lineages coinciding in the same bin is a rare event even when a shared constraint "
 "exists, so its null result is not evidence of absence.",

 "We therefore asked the same biological question with a test that uses every site rather "
 "than only those that coincide. Adaptive sites in entry glycoproteins are displaced "
 "toward the N-terminus: their mean relative position is 0.298, against 0.500 expected if "
 "they were distributed as the testable positions are (9 sites, 4 lineages, 3 families; "
 "permutation p = 0.018, q = 0.036 after correction across the three functional roles). "
 "Polymerase sites show no such displacement and if anything lean toward the C-terminus "
 "(15 sites, 6 lineages, 3 families; 0.610; p = 0.930, q = 0.930), which serves as an "
 "internal negative control: a displacement present in both roles would indicate the "
 "site-selection rule or the alignment geometry rather than a property of entry proteins. "
 "This analysis was performed post hoc, after the bin co-occurrence test proved "
 "underpowered, and is reported as such."],

337: [
 "No protein domain was enriched for adaptive sites after correction for multiple testing. "
 "The strongest nominal signal was the glycoprotein G stalk domain (4 of 598 residues; "
 "p = 0.060, q = 0.391), and the L-protein N-terminal domain, previously reported as the "
 "strongest enrichment, showed no excess at all (4 of 4792 residues; odds ratio 0.98, "
 "p = 1.0). The earlier estimate was inflated by sequences subsequently found not to "
 "belong to the target gene."],

338: [
 "Recombination breakpoints were detected in four alignments: Andes GnGc (one, at codon "
 "747), Andes L (three, at 284, 629 and 688), Marburg L (two, at 1689 and 1866) and Sin "
 "Nombre GnGc (two, at 291 and 325). None of the 26 differential sites fell within ten "
 "codons of any of them. Three alignments were not fully screened, and their recombination "
 "status is unknown rather than clean."],

341: [
 "The mechanism of human-to-human adaptation is lineage-specific and reflects the receptor "
 "biology of each viral family. In henipaviruses, positive selection at discrete entry "
 "residues is accompanied by a dual-property signal at polymerase site L210, where "
 "hydrophobicity is conserved and isoelectric point diversifies (p = 0.015 and p = 0.017; "
 "omnibus p = 0.031). It is the only site in this alignment with any significant "
 "physicochemical property, but it does not survive correction across the alignment "
 "(omnibus q = 1.0) and is therefore a candidate rather than an established finding. In "
 "hantaviruses, fluid-phase macropinocytic entry combined with the use of conserved "
 "cellular receptors such as PCDH1 may pre-adapt Andes virus to human cells without "
 "requiring host-specific receptor optimisation. We note that the relaxation parameters "
 "previously cited in support of an entry-relaxed and polymerase-intensified dichotomy do "
 "not hold as a systematic pattern across the six lineages that support a branch contrast "
 "(paired Wilcoxon, p = 1.0)."],
}

CAPTIONS = {
163: "Figure 1 Adaptive sites per 1000 codons in H2H-capable and spillover-only lineages, "
     "grouped by virus family. Bar labels give the absolute number of sites meeting the "
     "consensus criterion for positive selection. H2H-capable lineages do not carry more "
     "adaptive sites than their spillover-only relatives; the pooled rate is lower (0.64 "
     "against 1.38 per 1000 codons; Fisher exact odds ratio 0.46, p = 0.077). Within "
     "Paramyxoviridae the two rates are indistinguishable.",
252: "Figure 2 PRIME physicochemical selection profile at Nipah NiV-B L-protein site 210. "
     "Positive lambda: property conserved. Negative lambda: property diversifying. "
     "Hydrophobicity is conserved (lambda = +3.75, p = 0.015) while isoelectric point "
     "diversifies (lambda = -4.90, p = 0.017); volume is unchanged (p = 0.398). The "
     "omnibus test at this site is nominally significant (p = 0.031) but does not survive "
     "correction across the alignment (q = 1.0), so it is reported as a candidate rather "
     "than an established finding.",
332: "Figure 4 Relative position of differential sites along the protein, by functional "
     "role. Each point is one site; the vertical bar marks the mean; the dashed line marks "
     "the position expected if sites were distributed as the testable positions are. Entry "
     "glycoprotein sites are displaced toward the N-terminus (mean 0.298; q = 0.036), "
     "polymerase sites are not (0.610; q = 0.930). q values are permutation-based and "
     "corrected across the three functional roles examined. The entry helper role, "
     "comprising the Nipah F protein, contributed two sites and was not tested.",
}

# Figure 2's caption shares a paragraph with its image, so only text runs there
# may be struck through.
TEXT_ONLY_STRIKE = {252}

INSERTIONS = {
295: ["AUTHOR NOTE, DELETE BEFORE SUBMISSION: Figure 3 is a screenshot of the interactive 3D "
      "viewer and could not be regenerated automatically. Open the viewer from the current "
      "run (05_structure/Andes_virus/GnGc_3d_view.html), which now shows three evidence "
      "tiers in separate colours rather than Contrast-FEL sites alone, and take a fresh "
      "screenshot."],
343: [
 "The number of differential sites per lineage is small, between zero and nine, and the "
 "study is correspondingly underpowered. The absence of a distinguishing signature between "
 "H2H-capable and spillover-only lineages should be read as not detectable at this sample "
 "size, not as established absence. The comparison rests on eight lineages in three "
 "families; a larger panel would be required to exclude an effect of the size that "
 "surveillance applications would need.",

 "Public sequence archives sample H2H-capable and spillover-only lineages unevenly, and "
 "sampling is driven by outbreak investigation rather than by design. Reservoir sequences "
 "in particular are sparse, and for two lineages absent altogether.",

 "Three of eighteen alignments were not fully screened for recombination. Their status is "
 "unknown rather than clean, and results from those alignments should be read with the "
 "possibility of recombination-driven false positives in mind.",

 "The positional analysis is post hoc. It tests the same hypothesis as the pre-specified "
 "bin co-occurrence analysis but uses the data more efficiently; it was selected after the "
 "co-occurrence test proved underpowered, and its p value has not been adjusted for that "
 "choice."],
}


def esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def strike_paragraph(p, text_only=False):
    """
    Grey out and strike through the runs in this paragraph.

    With text_only, runs holding a drawing are left alone: Figure 2's caption
    shares a paragraph with its image, and striking the whole paragraph would
    strike the artwork too.
    """
    if text_only:
        parts = re.split(r"(<w:r>.*?</w:r>|<w:r [^>]*>.*?</w:r>)", p, flags=re.S)
        return "".join(
            x if ("<w:drawing" in x or not x.startswith("<w:r")) else strike_paragraph(x)
            for x in parts)

    mark = f'<w:strike/><w:color w:val="{GREY}"/>'
    p = re.sub(r"(<w:rPr>)", r"\1" + mark, p)
    p = re.sub(r"(<w:r>|<w:r [^>]*>)(?!<w:rPr>)",
               lambda m: m.group(0) + f"<w:rPr>{mark}</w:rPr>", p)
    return p


def blue_paragraph(text):
    return ('<w:p><w:pPr><w:spacing w:after="160" w:line="276" w:lineRule="auto"/></w:pPr>'
            f'<w:r><w:rPr><w:color w:val="{BLUE}"/><w:sz w:val="22"/>'
            f'<w:szCs w:val="22"/></w:rPr>'
            f'<w:t xml:space="preserve">{esc(text)}</w:t></w:r></w:p>')


def rescale_drawings(xml, rid, w, h):
    """Fit the display box of this relationship's drawing to a w:h image."""
    def fix(m):
        block = m.group(0)
        if f'r:embed="{rid}"' not in block:
            return block
        cx = int(re.search(r'<wp:extent cx="(\d+)"', block).group(1))
        cy = int(round(cx * h / w))
        block = re.sub(r'<wp:extent cx="\d+" cy="\d+"/>',
                       f'<wp:extent cx="{cx}" cy="{cy}"/>', block)
        block = re.sub(r'<a:ext cx="\d+" cy="\d+"/>',
                       f'<a:ext cx="{cx}" cy="{cy}"/>', block)
        return block
    return re.sub(r"<w:drawing>.*?</w:drawing>", fix, xml, flags=re.S)


def main(src, dst):
    zin = zipfile.ZipFile(src)
    xml = zin.read("word/document.xml").decode("utf-8")

    spans = [(m.start(), m.end(), m.group(0))
             for m in re.finditer(r"<w:p[ >].*?</w:p>", xml, re.S)]
    out, prev, n_rep, n_cap, n_ins = [], 0, 0, 0, 0
    for i, (a, b, para) in enumerate(spans):
        out.append(xml[prev:a])
        if i in REPLACEMENTS or i in CAPTIONS:
            out.append(strike_paragraph(para, text_only=i in TEXT_ONLY_STRIKE))
            for t in REPLACEMENTS.get(i, []):
                out.append(blue_paragraph(t))
            if i in CAPTIONS:
                out.append(blue_paragraph(CAPTIONS[i])); n_cap += 1
            if i in REPLACEMENTS:
                n_rep += 1
        elif i in INSERTIONS:
            out.append(para)
            out.extend(blue_paragraph(t) for t in INSERTIONS[i])
            n_ins += 1
        else:
            out.append(para)
        prev = b
    out.append(xml[prev:])
    new_xml = "".join(out)

    rels = zin.read("word/_rels/document.xml.rels").decode()
    rid_of = {t: r for r, t in re.findall(r'Id="(rId\d+)"[^>]*Target="(media/[^"]+)"', rels)}
    for part, src_png in MEDIA.items():
        rid = rid_of.get(part.replace("word/", ""))
        if not rid:
            continue
        blob = open(src_png, "rb").read()
        w, h = struct.unpack(">II", blob[16:24])
        new_xml = rescale_drawings(new_xml, rid, w, h)

    with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                data = new_xml.encode("utf-8")
            elif item.filename in MEDIA:
                data = open(MEDIA[item.filename], "rb").read()
            zout.writestr(item, data)

    print(f"paragraf teks diganti : {n_rep}")
    print(f"caption diganti       : {n_cap}")
    print(f"titik sisipan         : {n_ins}")
    print(f"gambar ditukar        : {len(MEDIA)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
