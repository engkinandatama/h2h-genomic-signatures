# Convergent Genomic Signatures of Zoonotic Viruses

An end-to-end automated Snakemake pipeline for comparative genomic analysis of Zoonotic RNA viruses (Hantavirus, Nipah, and Ebola). This pipeline aims to identify convergent genomic signatures that distinguish sustained Human-to-Human (H2H) transmission from isolated animal reservoir spillover events.

## 🧬 Pipeline Architecture

The pipeline is fully automated and modular, consisting of 6 main stages:

1. **Fetch & QC**: Automatically downloads complete genomes and protein sequences from NCBI using `datasets` CLI, filtering by geographic location and host species (e.g., *Homo sapiens* vs *Reservoir*).
2. **Alignment**: Performs multiple sequence alignment using MAFFT, followed by codon alignment via `pal2nal` to prepare for evolutionary analysis.
3. **Phylogeny**: Reconstructs Maximum Likelihood phylogenetic trees using IQ-TREE2.
4. **Selection Analysis**: Detects episodic and pervasive positive selection (MEME, FEL, FUBAR) on specific lineages using HyPhy.
5. **Structural Mapping**: Maps adaptive mutation sites onto 3D protein structures (AlphaFold PDBs) generating interactive HTML viewers via Py3Dmol.
6. **Convergence Statistics**: Aggregates all viral data to perform statistical tests (Mann-Whitney, Jaccard Index) yielding a final convergent signature matrix.

## 📂 Repository Structure

```
.
├── config/
│   └── config.yaml          # Single Source of Truth for pipeline parameters
├── docs/
│   └── usage.md             # Execution instructions (Local & HPC)
├── workflow/
│   ├── envs/                # Conda environments for reproducibility
│   ├── rules/               # Modular Snakemake rules
│   ├── scripts/             # Python/BioPython processing scripts
│   └── Snakefile            # Main pipeline entrypoint
```

## 🚀 Quick Start

Please read the **[Usage Documentation](docs/usage.md)** for detailed instructions on configuring and running the pipeline locally, on a standalone HPC node, or via a SLURM workload manager.
