# Convergent Genomic Signatures of Zoonotic Viruses

An end-to-end automated Snakemake pipeline for comparative genomic analysis of Zoonotic RNA viruses (Hantavirus, Nipah, and Ebola). This pipeline aims to identify convergent genomic signatures that distinguish sustained Human-to-Human (H2H) transmission from isolated animal reservoir spillover events.

## 🧬 Pipeline Architecture

The pipeline is fully automated and modular, consisting of 7 stages:

1. **Fetch & QC**: Downloads sequences from three sources — NCBI (`datasets` CLI), BV-BRC (REST API) and, where configured, a manually exported GISAID set whose CDS are located by ORF search against a UniProt reference. Records are filtered by host and geography, then deduplicated across sources by exact nucleotide match. Two length gates apply: an absolute floor and a per-protein relative gate that removes partial fragments.
2. **Alignment**: Performs multiple sequence alignment using MAFFT, followed by codon alignment via `pal2nal` to prepare for evolutionary analysis.
3. **Phylogeny**: Reconstructs Maximum Likelihood phylogenetic trees using IQ-TREE2.
4. **Selection Analysis**: Ten HyPhy analyses. Site-level: MEME (episodic), FEL and FUBAR (pervasive), SLAC, PRIME (physicochemical properties). Branch-contrast: Contrast-FEL and RELAX, run only for groups that contain both a human-derived and a reservoir-derived arm. Gene-wide: BUSTED and aBSREL. Recombination: GARD, whose breakpoints are used to flag nearby sites.

   FEL results are direction-filtered (beta > alpha) before a site counts toward positive selection; its p-value is two-sided, so significance alone also captures purifying sites. Sites are called consensus when at least `min_methods` site-level analyses agree; MEME-only hits are reported separately as episodic candidates because MEME tests a different hypothesis from FEL and FUBAR.
5. **Structural Mapping**: Maps adaptive mutation sites onto 3D protein structures (AlphaFold PDBs) generating interactive HTML viewers via Py3Dmol.
6. **Convergence Statistics**: Maps differential sites onto a 100-bin normalised coordinate and tests bin co-occurrence against a permutation null, so a shared bin carries a p-value rather than only a count. Domain enrichment is conditioned on sites the branch contrast could actually reach. Both uncorrected and FDR-corrected results are written.

7. **Figures**: Per-protein Manhattan plots.

## 📂 Repository Structure

```
.
├── config/
│   └── config.yaml          # Single Source of Truth for pipeline parameters
├── docs/
│   └── usage.md             # Execution instructions (Local & HPC)
├── workflow/
│   ├── envs/                # Conda environments, versions pinned
│   ├── rules/               # Modular Snakemake rules
│   ├── scripts/             # Python/BioPython processing scripts
│   └── Snakefile            # Main pipeline entrypoint
```

## 🚀 Quick Start

Please read the **[Usage Documentation](docs/usage.md)** for detailed instructions on configuring and running the pipeline locally, on a standalone HPC node, or via a SLURM workload manager.
