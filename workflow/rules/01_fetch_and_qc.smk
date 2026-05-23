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

rule filter_and_extract_cds:
    input:
        zip=WORKDIR + "/01_raw_fasta/{virus}/{virus}_dataset.zip"
    output:
        fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_filtered.fasta",
        faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_filtered.faa"
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
        echo "Extracting CDS and protein translation for {wildcards.virus} - {wildcards.protein}" > {log}
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
