# Rule 2: Alignment (modified for virus_groups)
ruleorder: codon_alignment > align_proteins

rule align_proteins:
    input:
        fasta=WORKDIR + "/02_aligned/{virus_group}/{protein}_merged.faa"
    output:
        msa=WORKDIR + "/02_aligned/{virus_group}/{protein}_aligned.faa"
    log:
        WORKDIR + "/logs/align_proteins/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/align_proteins/{virus_group}_{protein}.tsv"
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
        echo "Starting MAFFT alignment for {wildcards.virus_group} - {wildcards.protein}" > {log}
        
        # Check if merged file is empty or has < 4 sequences
        # grep -c prints 0 and ALSO exits 1 when nothing matches, so `|| echo 0`
        # appended a second line and n_seqs became "0\n0". The -lt test then failed
        # with "integer expression expected", the skip branch was never taken, and
        # mafft was handed an empty file.
        n_seqs=$(grep -c "^>" {input.fasta} 2>/dev/null | head -1)
        n_seqs=${{n_seqs:-0}}
        if [ "$n_seqs" -lt 4 ]; then
            echo "SKIP: Hanya $n_seqs sekuens ditemukan. Membuat file output kosong." >> {log}
            touch {output.msa}
            exit 0
        fi
        
        mafft --thread {threads} {params.mafft_args} {input.fasta} > {output.msa} 2>> {log}
        echo "Alignment complete." >> {log}
        """

rule codon_alignment:
    input:
        msa=WORKDIR + "/02_aligned/{virus_group}/{protein}_aligned.faa",
        cds=WORKDIR + "/02_aligned/{virus_group}/{protein}_merged.fasta"
    output:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta"
    log:
        WORKDIR + "/logs/codon_alignment/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/codon_alignment/{virus_group}_{protein}.tsv"
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        echo "Starting pal2nal codon alignment for {wildcards.virus_group} - {wildcards.protein}" > {log}
        
        # Check minimum sequences before running pal2nal
        if [ ! -s "{input.msa}" ]; then
            echo "SKIP: Input MSA kosong. Membuat file output kosong." >> {log}
            touch {output.codon_aln}
            exit 0
        fi
        
        n_seqs=$(grep -c "^>" {input.msa} 2>/dev/null || echo "0")
        echo "Jumlah sekuens dalam MSA: $n_seqs" >> {log}
        if [ "$n_seqs" -lt 4 ]; then
            echo "SKIP: Hanya $n_seqs sekuens ditemukan (minimum 4 diperlukan untuk analisis filogenetik). Membuat file output kosong." >> {log}
            touch {output.codon_aln}
            exit 0
        fi
        
        pal2nal.pl {input.msa} {input.cds} -output fasta > {output.codon_aln} 2>> {log}
        echo "Codon alignment complete." >> {log}
        """
