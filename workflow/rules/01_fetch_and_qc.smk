# Rule 1: Fetch and QC
rule download_genomes:
    output:
        zip=WORKDIR + "/01_raw_fasta/{virus}/{virus}_dataset.zip"
    log:
        WORKDIR + "/logs/download_genomes/{virus}.log"
    benchmark:
        WORKDIR + "/benchmarks/download_genomes/{virus}.tsv"
    params:
        taxon_id=lambda wildcards: config["viruses"][wildcards.virus]["taxon_id"],
        ncbi_api_key=os.environ.get("NCBI_API_KEY", "")
    resources:
        ncbi_api=1
    conda:
        "../envs/download.yaml"
    shell:
        """
        if [ -n "{params.ncbi_api_key}" ]; then
            export NCBI_API_KEY="{params.ncbi_api_key}"
            echo "Using NCBI_API_KEY from environment" > {log}
        else
            echo "Warning: NCBI_API_KEY is not set" > {log}
        fi
        
        echo "Starting download for {wildcards.virus} (Taxon: {params.taxon_id})" >> {log}
        datasets download virus genome taxon {params.taxon_id} \
            --include genome,protein,cds \
            --filename {output.zip} >> {log} 2>&1
        
        if [ -f "{output.zip}" ]; then
            echo "Download complete. File size: $(du -sh {output.zip} | cut -f1)" >> {log}
            echo "Contents of downloaded ZIP:" >> {log}
            unzip -l {output.zip} >> {log} 2>&1 || echo "Warning: unzip -l failed" >> {log}
        else
            echo "Error: Downloaded ZIP file {output.zip} was not created!" >> {log}
        fi
        """

rule filter_and_extract_cds_ncbi:
    input:
        zip=WORKDIR + "/01_raw_fasta/{virus}/{virus}_dataset.zip"
    output:
        fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_ncbi.fasta",
        faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_ncbi.faa"
    log:
        WORKDIR + "/logs/filter_extract_ncbi/{virus}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/filter_extract_ncbi/{virus}_{protein}.tsv"
    params:
        geo_filter=lambda wildcards: config["viruses"][wildcards.virus].get("geo_filter", ""),
        host_filter=lambda wildcards: config["viruses"][wildcards.virus].get("host_filter", ""),
        max_seq=config["max_sequences_per_group"],
        min_len=config.get("min_length_cds", 1500)
    conda:
        "../envs/download.yaml"
    shell:
        """
        echo "Extracting CDS and protein translation for {wildcards.virus} - {wildcards.protein} (NCBI)" > {log}
        python workflow/scripts/extract_cds.py \
            --zip {input.zip} \
            --protein {wildcards.protein} \
            --geo "{params.geo_filter}" \
            --host "{params.host_filter}" \
            --max {params.max_seq} \
            --min-len {params.min_len} \
            --out-nuc {output.fasta} \
            --out-prot {output.faa} >> {log} 2>&1
        echo "Extraction and filtering complete." >> {log}
        """

rule fetch_bvbrc:
    output:
        fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_bvbrc.fasta",
        faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_bvbrc.faa"
    log:
        WORKDIR + "/logs/fetch_bvbrc/{virus}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/fetch_bvbrc/{virus}_{protein}.tsv"
    params:
        taxon_id=lambda wildcards: config["viruses"][wildcards.virus]["taxon_id"],
        geo_filter=lambda wildcards: config["viruses"][wildcards.virus].get("geo_filter", ""),
        host_filter=lambda wildcards: config["viruses"][wildcards.virus].get("host_filter", ""),
        max_seq=config["max_sequences_per_group"],
        min_len=config.get("min_length_cds", 1500)
    conda:
        "../envs/download.yaml"
    shell:
        """
        echo "Fetching from BV-BRC for {wildcards.virus} - {wildcards.protein}" > {log}
        python workflow/scripts/fetch_bvbrc.py \
            --taxon "{params.taxon_id}" \
            --protein "{wildcards.protein}" \
            --geo "{params.geo_filter}" \
            --host "{params.host_filter}" \
            --max {params.max_seq} \
            --min-len {params.min_len} \
            --out-nuc {output.fasta} \
            --out-prot {output.faa} >> {log} 2>&1
        echo "BV-BRC fetch complete." >> {log}
        """

rule merge_and_deduplicate:
    input:
        ncbi_fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_ncbi.fasta",
        ncbi_faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_ncbi.faa",
        bvbrc_fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_bvbrc.fasta",
        bvbrc_faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_bvbrc.faa"
    output:
        fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_filtered.fasta",
        faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_filtered.faa"
    log:
        WORKDIR + "/logs/merge_and_dedup/{virus}_{protein}.log"
    params:
        max_seq=config["max_sequences_per_group"]
    conda:
        "../envs/download.yaml"
    shell:
        """
        echo "Merging and deduplicating nucleotides and proteins together..." > {log}
        python workflow/scripts/merge_fasta.py \
            --ncbi-nuc {input.ncbi_fasta} \
            --ncbi-prot {input.ncbi_faa} \
            --bvbrc-nuc {input.bvbrc_fasta} \
            --bvbrc-prot {input.bvbrc_faa} \
            --out-nuc {output.fasta} \
            --out-prot {output.faa} \
            --max {params.max_seq} >> {log} 2>&1
        echo "Merge and paired deduplication complete." >> {log}
        """
