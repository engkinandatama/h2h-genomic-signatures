# Rule 3: Phylogenetics
rule build_tree:
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus}/{protein}_codon_aligned.fasta"
    output:
        tree=WORKDIR + "/03_trees/{virus}/{protein}_tree.treefile"
    log:
        WORKDIR + "/logs/build_tree/{virus}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/build_tree/{virus}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("iqtree", 16)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("iqtree", 16000)
    params:
        iqtree_args=config["params"]["iqtree"],
        prefix=WORKDIR + "/03_trees/{virus}/{protein}_tree"
    conda:
        "../envs/phylogeny.yaml"
    shell:
        """
        echo "Starting phylogenetic tree construction for {wildcards.virus} - {wildcards.protein} with {threads} threads" > {log}
        
        # Check if iqtree2 or iqtree is available in conda environment
        if command -v iqtree2 >/dev/null 2>&1; then
            IQTREE_CMD="iqtree2"
        elif command -v iqtree >/dev/null 2>&1; then
            IQTREE_CMD="iqtree"
        else
            echo "Error: Neither iqtree2 nor iqtree was found in the environment." >> {log}
            exit 127
        fi
        
        echo "Using binary: $IQTREE_CMD" >> {log}
        $IQTREE_CMD -s {input.codon_aln} --prefix {params.prefix} {params.iqtree_args} -T {threads} >> {log} 2>&1
        echo "Tree building complete." >> {log}
        """
