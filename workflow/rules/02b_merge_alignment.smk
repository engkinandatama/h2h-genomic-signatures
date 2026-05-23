# Rule 2b: Merge sequences per virus group
rule merge_sequences:
    input:
        nuc=lambda wildcards: [f"{WORKDIR}/01_raw_fasta/{v}/{wildcards.protein}_filtered.fasta" for v in config["virus_groups"][wildcards.virus_group]["viruses"]],
        prot=lambda wildcards: [f"{WORKDIR}/01_raw_fasta/{v}/{wildcards.protein}_filtered.faa" for v in config["virus_groups"][wildcards.virus_group]["viruses"]]
    output:
        merged_nuc=WORKDIR + "/02_aligned/{virus_group}/{protein}_merged.fasta",
        merged_prot=WORKDIR + "/02_aligned/{virus_group}/{protein}_merged.faa"
    log:
        WORKDIR + "/logs/merge_sequences/{virus_group}_{protein}.log"
    params:
        viruses=lambda wildcards: " ".join(config["virus_groups"][wildcards.virus_group]["viruses"])
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        echo "Merging sequences for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        python workflow/scripts/merge_sequences.py \
            --inputs {input.nuc} \
            --virus-names {params.viruses} \
            --output {output.merged_nuc} >> {log} 2>&1
            
        python workflow/scripts/merge_sequences.py \
            --inputs {input.prot} \
            --virus-names {params.viruses} \
            --output {output.merged_prot} >> {log} 2>&1
        """
