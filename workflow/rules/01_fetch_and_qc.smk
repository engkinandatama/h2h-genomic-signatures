# Rule 1: Fetch and QC
rule download_genomes:
    output:
        zip=WORKDIR + "/01_raw_fasta/{virus}/{virus}_dataset.zip"
    log:
        WORKDIR + "/logs/download_genomes/{virus}.log"
    benchmark:
        WORKDIR + "/benchmarks/download_genomes/{virus}.tsv"
    params:
        taxon_id=lambda wildcards: config["viruses"][wildcards.virus]["taxon_id"]
    conda:
        "../envs/download.yaml"
    shell:
        """
        echo "Starting download for {wildcards.virus} (Taxon: {params.taxon_id})" > {log}
        datasets download virus genome taxon {params.taxon_id} \
            --include genome,protein,cds,annotation \
            --filename {output.zip} >> {log} 2>&1
        echo "Download complete." >> {log}
        """

rule filter_and_extract_cds:
    input:
        zip=WORKDIR + "/01_raw_fasta/{virus}/{virus}_dataset.zip"
    output:
        fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_filtered.fasta"
    log:
        WORKDIR + "/logs/filter_extract/{virus}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/filter_extract/{virus}_{protein}.tsv"
    params:
        geo_filter=lambda wildcards: config["viruses"][wildcards.virus].get("geo_filter", ""),
        host_filter=lambda wildcards: config["viruses"][wildcards.virus].get("host_filter", ""),
        max_seq=config["max_sequences_per_group"],
        min_len=config.get("min_length_cds", 1500)
    conda:
        "../envs/download.yaml"
    shell:
        """
        echo "Extracting CDS for {wildcards.virus} - {wildcards.protein}" > {log}
        python workflow/scripts/extract_cds.py \
            --zip {input.zip} \
            --protein {wildcards.protein} \
            --geo "{params.geo_filter}" \
            --host "{params.host_filter}" \
            --max {params.max_seq} \
            --min-len {params.min_len} \
            --out {output.fasta} >> {log} 2>&1
        echo "Extraction and filtering complete." >> {log}
        """
