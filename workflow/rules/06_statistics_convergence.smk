# Rule 6: Statistical Analysis & Convergence (Agregasi Seluruh Virus)
rule calculate_convergence:
    input:
        # Mengambil SEMUA hasil seleksi dari SEMUA virus dan protein
        selection_results=expand(WORKDIR + "/04_selection/{virus}/{protein}_selection_results.txt",
                                 virus=VIRUSES,
                                 protein=[p for v in VIRUSES for p in PROTEINS[v]])
    output:
        matrix=WORKDIR + "/06_statistics/convergent_signature_matrix.csv",
        plot_jaccard=WORKDIR + "/06_statistics/jaccard_similarity_heatmap.png",
        plot_correlation=WORKDIR + "/06_statistics/selection_phenotype_correlation.png"
    log:
        WORKDIR + "/logs/calculate_convergence/all_viruses.log"
    benchmark:
        WORKDIR + "/benchmarks/calculate_convergence/all_viruses.tsv"
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
