# Rule 6: Statistical Analysis & Convergence (Agregasi Seluruh Virus Group)

rule calculate_convergence:
    input:
        # Mengambil hasil seleksi dari seluruh virus_group dan protein secara dinamis dari config
        selection_results=lambda wildcards: [
            f"{WORKDIR}/04_selection/{vg}/{p}_selection_results.txt"
            for vg in config["virus_groups"].keys()
            for p in config["virus_groups"][vg]["proteins"]
        ]
    output:
        matrix=WORKDIR + "/06_statistics/convergent_signature_matrix.csv",
        plot_jaccard=WORKDIR + "/06_statistics/jaccard_similarity_heatmap.png",
        plot_correlation=WORKDIR + "/06_statistics/selection_phenotype_correlation.png"
    log:
        WORKDIR + "/logs/calculate_convergence/all_virus_groups.log"
    benchmark:
        WORKDIR + "/benchmarks/calculate_convergence/all_virus_groups.tsv"
    conda:
        "../envs/statistics.yaml"
    shell:
        """
        echo "Starting Level 1-4 Statistical Analysis..." > {log}
        python workflow/scripts/calculate_convergence.py \
            --inputs {input.selection_results} \
            --out_matrix {output.matrix} \
            --out_jaccard {output.plot_jaccard} \
            --out_correlation {output.plot_correlation} >> {log} 2>&1
        echo "Statistical convergence analysis complete." >> {log}
        """


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
