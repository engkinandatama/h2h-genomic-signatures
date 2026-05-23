# Rule 2: Alignment
ruleorder: codon_alignment > align_proteins

rule align_proteins:
    input:
        fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_filtered.faa"
    output:
        msa=WORKDIR + "/02_aligned/{virus}/{protein}_aligned.faa"
    log:
        WORKDIR + "/logs/align_proteins/{virus}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/align_proteins/{virus}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("mafft", 4)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("mafft", 4000)
    params:
        mafft_args=config["params"]["mafft"]
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        echo "Starting MAFFT alignment for {wildcards.virus} - {wildcards.protein}" > {log}
        mafft --thread {threads} {params.mafft_args} {input.fasta} > {output.msa} 2>> {log}
        echo "Alignment complete." >> {log}
        """

rule codon_alignment:
    input:
        msa=WORKDIR + "/02_aligned/{virus}/{protein}_aligned.faa",
        cds=WORKDIR + "/01_raw_fasta/{virus}/{protein}_filtered.fasta"
    output:
        codon_aln=WORKDIR + "/02_aligned/{virus}/{protein}_codon_aligned.fasta"
    log:
        WORKDIR + "/logs/codon_alignment/{virus}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/codon_alignment/{virus}_{protein}.tsv"
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        echo "Starting pal2nal codon alignment for {wildcards.virus} - {wildcards.protein}" > {log}
        pal2nal.pl {input.msa} {input.cds} -output fasta > {output.codon_aln} 2>> {log}
        echo "Codon alignment complete." >> {log}
        """
