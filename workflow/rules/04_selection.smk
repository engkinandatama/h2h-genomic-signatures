# Rule 4: Positive Selection Analysis
rule hyphy_selection:
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus}/{protein}_tree.treefile"
    output:
        results=WORKDIR + "/04_selection/{virus}/{protein}_selection_results.txt"
    log:
        WORKDIR + "/logs/hyphy_selection/{virus}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_selection/{virus}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000)
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy for {wildcards.virus} - {wildcards.protein}..." > {log}
        
        # LABELING TREE (Poin 3: Memisahkan H2H vs Reservoir)
        labeled_tree="results/03_trees/{wildcards.virus}/{wildcards.protein}_labeled.treefile"
        python workflow/scripts/label_tree.py \
            --tree {input.tree} \
            --out $labeled_tree >> {log} 2>&1
            
        echo "Simulated selection results." > {output.results}
        
        # Contoh pemanggilan yang sebenarnya menggunakan tree yang sudah dilabeli:
        # hyphy fubar CPU={threads} --alignment {input.codon_aln} --tree $labeled_tree >> {log} 2>&1
        # hyphy fel CPU={threads} --alignment {input.codon_aln} --tree $labeled_tree >> {log} 2>&1
        # hyphy meme CPU={threads} --alignment {input.codon_aln} --tree $labeled_tree >> {log} 2>&1
        
        echo "Selection analysis complete." >> {log}
        """
