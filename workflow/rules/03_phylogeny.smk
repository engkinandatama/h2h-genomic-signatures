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
        echo "Starting IQ-TREE2 for {wildcards.virus} - {wildcards.protein} with {threads} threads" > {log}
        iqtree2 -s {input.codon_aln} --prefix {params.prefix} {params.iqtree_args} -T {threads} >> {log} 2>&1
        echo "Tree building complete." >> {log}
        """
