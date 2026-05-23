# Rule 4: Positive Selection Analysis (modified for virus_groups)

rule label_tree:
    input:
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_tree.treefile"
    output:
        labeled_tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    log:
        WORKDIR + "/logs/label_tree/{virus_group}_{protein}.log"
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running tree labeling for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ]; then
            echo "SKIP: File treefile kosong. Membuat file labeled_tree kosong." >> {log}
            touch {output.labeled_tree}
            exit 0
        fi
        python workflow/scripts/label_tree.py \
            --tree {input.tree} \
            --out {output.labeled_tree} >> {log} 2>&1
        """

rule hyphy_meme:
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_meme.json"
    log:
        WORKDIR + "/logs/hyphy_meme/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_meme/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000)
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy MEME for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong. Membuat output JSON kosong." >> {log}
            echo '{{}}' > {output.json}
            exit 0
        fi
        
        hyphy meme \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --output {output.json} \
            --cpu {threads} \
            --branches Foreground >> {log} 2>&1 || (echo '{{}}' > {output.json} && echo "Warning: HyPhy MEME failed, created empty JSON" >> {log})
        """

rule hyphy_fel:
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_fel.json"
    log:
        WORKDIR + "/logs/hyphy_fel/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_fel/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000)
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy FEL for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong. Membuat output JSON kosong." >> {log}
            echo '{{}}' > {output.json}
            exit 0
        fi
        
        hyphy fel \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --output {output.json} \
            --cpu {threads} \
            --branches Foreground >> {log} 2>&1 || (echo '{{}}' > {output.json} && echo "Warning: HyPhy FEL failed, created empty JSON" >> {log})
        """

rule hyphy_fubar:
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_tree.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_fubar.json"
    log:
        WORKDIR + "/logs/hyphy_fubar/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_fubar/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000)
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy FUBAR for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong. Membuat output JSON kosong." >> {log}
            echo '{{}}' > {output.json}
            exit 0
        fi
        
        hyphy fubar \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --output {output.json} \
            --cpu {threads} >> {log} 2>&1 || (echo '{{}}' > {output.json} && echo "Warning: HyPhy FUBAR failed, created empty JSON" >> {log})
        """

rule aggregate_selection:
    input:
        meme=WORKDIR + "/04_selection/{virus_group}/{protein}_meme.json",
        fel=WORKDIR + "/04_selection/{virus_group}/{protein}_fel.json",
        fubar=WORKDIR + "/04_selection/{virus_group}/{protein}_fubar.json"
    output:
        results=WORKDIR + "/04_selection/{virus_group}/{protein}_selection_results.txt",
        meme_txt=WORKDIR + "/04_selection/{virus_group}/{protein}_meme_results.txt",
        fel_txt=WORKDIR + "/04_selection/{virus_group}/{protein}_fel_results.txt",
        fubar_txt=WORKDIR + "/04_selection/{virus_group}/{protein}_fubar_results.txt",
        complete=WORKDIR + "/04_selection/{virus_group}/{protein}_selection_complete.tsv"
    log:
        WORKDIR + "/logs/aggregate_selection/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/aggregate_selection/{virus_group}_{protein}.tsv"
    params:
        pval=config.get("params", {}).get("hyphy", {}).get("pvalue_threshold", 0.05),
        fubar_pp=config.get("params", {}).get("hyphy", {}).get("fubar_pp_threshold", 0.90),
        min_methods=config.get("params", {}).get("hyphy", {}).get("min_methods", 2)
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Aggregating HyPhy results for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        python workflow/scripts/parse_hyphy.py \
            --meme {input.meme} \
            --fel {input.fel} \
            --fubar {input.fubar} \
            --pvalue {params.pval} \
            --fubar-pp {params.fubar_pp} \
            --min-methods {params.min_methods} \
            --out-consensus {output.results} \
            --out-meme {output.meme_txt} \
            --out-fel {output.fel_txt} \
            --out-fubar {output.fubar_txt} \
            --out-complete {output.complete} >> {log} 2>&1
        """
