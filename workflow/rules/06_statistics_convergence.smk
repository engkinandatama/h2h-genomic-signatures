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
