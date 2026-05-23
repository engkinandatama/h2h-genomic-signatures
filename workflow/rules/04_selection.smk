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
    params:
        category=lambda wildcards: config["viruses"][wildcards.virus].get("category", ""),
        pvalue=config.get("params", {}).get("hyphy", {}).get("pvalue_threshold", 0.05)
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy analysis for {wildcards.virus} - {wildcards.protein}..." > {log}
        
        # 1. Melabeli pohon filogeni secara dinamis berdasarkan kategori virus
        labeled_tree="{WORKDIR}/03_trees/{wildcards.virus}/{wildcards.protein}_labeled.treefile"
        python workflow/scripts/label_tree.py \
            --tree {input.tree} \
            --category {params.category} \
            --out $labeled_tree >> {log} 2>&1
            
        # 2. Menjalankan HyPhy MEME yang sesungguhnya
        json_output="{WORKDIR}/04_selection/{wildcards.virus}/{wildcards.protein}_selection.json"
        
        # Panggil hyphy MEME
        hyphy meme CPU={threads} \
            --alignment {input.codon_aln} \
            --tree $labeled_tree \
            --output $json_output >> {log} 2>&1
            
        # 3. Parse output JSON dari HyPhy ke format hasil tabular
        python workflow/scripts/parse_hyphy.py \
            --json $json_output \
            --pvalue {params.pvalue} \
            --out {output.results} >> {log} 2>&1
            
        echo "Selection analysis complete." >> {log}
        """
