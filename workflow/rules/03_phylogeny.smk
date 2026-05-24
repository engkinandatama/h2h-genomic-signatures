# Rule 3: Phylogenetics (modified for virus_groups)
rule build_tree:
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta"
    output:
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_tree.treefile"
    log:
        WORKDIR + "/logs/build_tree/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/build_tree/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("iqtree", 16)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("iqtree", 16000),
        iqtree_jobs=1
    params:
        iqtree_args=config["params"]["iqtree"],
        prefix=WORKDIR + "/03_trees/{virus_group}/{protein}_tree"
    conda:
        "../envs/phylogeny.yaml"
    shell:
        """
        echo "Starting phylogenetic tree construction for {wildcards.virus_group} - {wildcards.protein} with {threads} threads" > {log}
        
        # Check minimum sequences before running IQ-TREE
        if [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input codon alignment kosong. Membuat file treefile kosong." >> {log}
            touch {output.tree}
            exit 0
        fi
        
        n_seqs=$(grep -c "^>" {input.codon_aln} 2>/dev/null || echo "0")
        echo "Jumlah sekuens dalam codon alignment: $n_seqs" >> {log}
        if [ "$n_seqs" -lt 4 ]; then
            echo "SKIP: Hanya $n_seqs sekuens ditemukan (IQ-TREE memerlukan minimal 4). Membuat file treefile kosong." >> {log}
            touch {output.tree}
            exit 0
        fi
        
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
        $IQTREE_CMD -s {input.codon_aln} --prefix {params.prefix} {params.iqtree_args} -T {threads} -redo >> {log} 2>&1
        echo "Tree building complete." >> {log}
        """
