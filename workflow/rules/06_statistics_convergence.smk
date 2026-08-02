# Rule 6: Statistical Analysis & Convergence (Agregasi Seluruh Virus Group)

# RETIRED 2026-07-25. Superseded by rule analyze_convergence below, which
# applies one stated site-selection rule, conditions the domain test on
# sites the contrast could reach, and tests bin co-occurrence against a
# permutation null. Kept commented so the provenance of the pre-rebuild
# numbers remains readable.
#
# rule calculate_convergence:
#     input:
#         # Mengambil hasil seleksi dari seluruh virus_group dan protein secara dinamis dari config
#         selection_results=lambda wildcards: [
#             f"{WORKDIR}/04_selection/{vg}/{p}_selection_results.txt"
#             for vg in config["virus_groups"].keys()
#             for p in config["virus_groups"][vg]["proteins"]
#         ]
#     output:
#         matrix=WORKDIR + "/06_statistics/convergent_signature_matrix.csv",
#         plot_jaccard=WORKDIR + "/06_statistics/jaccard_similarity_heatmap.png",
#         plot_correlation=WORKDIR + "/06_statistics/selection_phenotype_correlation.png"
#     log:
#         WORKDIR + "/logs/calculate_convergence/all_virus_groups.log"
#     benchmark:
#         WORKDIR + "/benchmarks/calculate_convergence/all_virus_groups.tsv"
#     conda:
#         "../envs/statistics.yaml"
#     shell:
#         """
#         echo "Starting Level 1-4 Statistical Analysis..." > {log}
#         python workflow/scripts/calculate_convergence.py \
#             --inputs {input.selection_results} \
#             --out_matrix {output.matrix} \
#             --out_jaccard {output.plot_jaccard} \
#             --out_correlation {output.plot_correlation} >> {log} 2>&1
#         echo "Statistical convergence analysis complete." >> {log}
#         """
#
#
rule visualize_selection_manhattan:
    """
    Generate per-protein Manhattan plot showing:
      - MEME + FEL -log10(p) across all codon positions
      - FUBAR posterior probability (secondary panel)
      - SLAC dN-dS (tertiary panel)
      - Contrast-FEL significant sites marked as red diamonds
      - PRIME physicochemical property annotations
    One PNG per virus_group x protein combination.
    """
    input:
        all_sites=WORKDIR + "/04_selection/{virus_group}/{protein}_all_sites.tsv"
    output:
        png=WORKDIR + "/07_figures/{virus_group}/{protein}_manhattan.png"
    log:
        WORKDIR + "/logs/visualize_manhattan/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/visualize_manhattan/{virus_group}_{protein}.tsv"
    conda:
        "../envs/statistics.yaml"
    shell:
        """
        echo "Generating Manhattan plot for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        python workflow/scripts/visualize_selection.py \
            --all-sites   {input.all_sites} \
            --out-png     {output.png} \
            --virus-group {wildcards.virus_group} \
            --protein     {wildcards.protein} >> {log} 2>&1
        echo "Manhattan plot complete." >> {log}
        """


rule analyze_convergence:
    """
    Cross-family convergence over the 100-bin coordinate.

    This was previously a standalone script under .dev/ that produced the
    manuscript's headline numbers and both main figures while sitting outside the
    workflow, so a clean `snakemake` run reproduced none of them. It now runs as a
    rule, applies one stated site-selection rule instead of hardcoded site lists,
    and tests bin co-occurrence against a permutation null.
    """
    input:
        all_sites=expand(
            WORKDIR + "/04_selection/{vg}/{p}_all_sites.tsv",
            zip,
            vg=[vg for vg in VIRUS_GROUPS
                for _ in config["virus_groups"][vg]["proteins"]],
            p=[p for vg in VIRUS_GROUPS
               for p in config["virus_groups"][vg]["proteins"]]),
        gard=expand(
            WORKDIR + "/04_selection/{vg}/{p}_gard_summary.txt",
            zip,
            vg=[vg for vg in VIRUS_GROUPS
                for _ in config["virus_groups"][vg]["proteins"]],
            p=[p for vg in VIRUS_GROUPS
               for p in config["virus_groups"][vg]["proteins"]])
    output:
        sites=WORKDIR + "/06_statistics/differential_sites.tsv",
        counts=WORKDIR + "/06_statistics/differential_site_counts.tsv",
        hotspots=WORKDIR + "/06_statistics/hotspot_bins.tsv",
        permutation=WORKDIR + "/06_statistics/hotspot_permutation_test.tsv",
        domains=WORKDIR + "/06_statistics/domain_enrichment.tsv",
        positional=WORKDIR + "/06_statistics/positional_bias_test.tsv",
        summary=WORKDIR + "/06_statistics/convergence_summary.json"
    log:
        WORKDIR + "/logs/analyze_convergence/all_groups.log"
    benchmark:
        WORKDIR + "/benchmarks/analyze_convergence/all_groups.tsv"
    params:
        outdir=WORKDIR + "/06_statistics",
        pval=config.get("params", {}).get("hyphy", {}).get("pvalue_threshold", 0.05),
        perms=config.get("params", {}).get("convergence_permutations", 20000),
        bp_margin=config.get("params", {}).get("breakpoint_margin_codons", 10)
    conda:
        "../envs/statistics.yaml"
    shell:
        """
        python workflow/scripts/analyze_convergence.py \
            --all-sites {input.all_sites} \
            --config config/config.yaml \
            --outdir {params.outdir} \
            --gard-summaries {input.gard} \
            --pvalue {params.pval} \
            --permutations {params.perms} \
            --breakpoint-margin {params.bp_margin} > {log} 2>&1
        """
