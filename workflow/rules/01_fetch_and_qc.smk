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
        min_len=config.get("min_length_cds", 1500),
        min_len_fraction=config.get("min_cds_length_fraction", 0.70)
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
            --min-len-fraction {params.min_len_fraction} \
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
        min_len=config.get("min_length_cds", 1500),
        min_len_fraction=config.get("min_cds_length_fraction", 0.70)
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
            --min-len-fraction {params.min_len_fraction} \
            --out-nuc {output.fasta} \
            --out-prot {output.faa} >> {log} 2>&1
        echo "BV-BRC fetch complete." >> {log}
        """

rule merge_and_deduplicate:
    input:
        ncbi_fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_ncbi.fasta",
        ncbi_faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_ncbi.faa",
        bvbrc_fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_bvbrc.fasta",
        bvbrc_faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_bvbrc.faa",
        gisaid_fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_gisaid.fasta",
        gisaid_faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_gisaid.faa"
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
            --extra-nuc {input.gisaid_fasta} \
            --extra-prot {input.gisaid_faa} \
            --out-nuc {output.fasta} \
            --out-prot {output.faa} \
            --max {params.max_seq} >> {log} 2>&1
        echo "Merge and paired deduplication complete." >> {log}
        """


# =============================================================================
# Optional third sequence source (GISAID)
# =============================================================================
# GISAID distributes whole genomes with no CDS annotation, so the NCBI path in
# extract_cds.py cannot read them. These rules fetch a UniProt reference protein
# for the target gene and use it to locate the corresponding ORF in each genome.
# Viruses whose taxon has no configured source produce empty files, which the
# merge step ignores.

rule fetch_reference_protein:
    output:
        faa=WORKDIR + "/00_reference/{virus}/{protein}_reference.faa"
    log:
        WORKDIR + "/logs/fetch_reference/{virus}_{protein}.log"
    params:
        uniprot_id=lambda w: reference_uniprot(w.virus, w.protein)
    conda:
        "../envs/download.yaml"
    shell:
        """
        mkdir -p $(dirname {output.faa})
        if [ -z "{params.uniprot_id}" ]; then
            echo "ERROR: no uniprot_id configured for {wildcards.virus}/{wildcards.protein}" > {log}
            exit 1
        fi
        echo "Fetching UniProt reference {params.uniprot_id}..." > {log}
        curl -fsSL --retry 3 --max-time 60 \
            "https://rest.uniprot.org/uniprotkb/{params.uniprot_id}.fasta" \
            -o {output.faa} 2>> {log}
        if [ ! -s "{output.faa}" ]; then
            echo "ERROR: UniProt returned nothing for {params.uniprot_id}" >> {log}
            exit 1
        fi
        """


rule fetch_gisaid:
    input:
        reference=WORKDIR + "/00_reference/{virus}/{protein}_reference.faa"
    output:
        fasta=WORKDIR + "/01_raw_fasta/{virus}/{protein}_gisaid.fasta",
        faa=WORKDIR + "/01_raw_fasta/{virus}/{protein}_gisaid.faa"
    log:
        WORKDIR + "/logs/fetch_gisaid/{virus}_{protein}.log"
    params:
        gisaid_fasta=lambda w: config.get("gisaid_sources", {}).get(
            str(config["viruses"][w.virus]["taxon_id"]), {}).get("fasta", ""),
        gisaid_meta=lambda w: config.get("gisaid_sources", {}).get(
            str(config["viruses"][w.virus]["taxon_id"]), {}).get("metadata", ""),
        host=lambda w: config["viruses"][w.virus].get("host_filter", ""),
        geo=lambda w: config["viruses"][w.virus].get("geo_filter", ""),
        gis=config.get("params", {}).get("gisaid", {})
    conda:
        "../envs/download.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        GENOMES="{params.gisaid_fasta}"
        META="{params.gisaid_meta}"

        if [ -z "$GENOMES" ] || [ ! -s "$GENOMES" ]; then
            echo "No GISAID source configured for {wildcards.virus}; writing empty outputs." > {log}
            : > {output.fasta}
            : > {output.faa}
            exit 0
        fi

        python workflow/scripts/extract_cds_from_genome.py \
            --genomes "$GENOMES" \
            --metadata "$META" \
            --reference-protein {input.reference} \
            --host-filter "{params.host}" \
            --geo-filter "{params.geo}" \
            --min-identity {params.gis[min_identity]} \
            --min-length-fraction {params.gis[min_length_fraction]} \
            --max-length-fraction {params.gis[max_length_fraction]} \
            --max-ambiguous-fraction {params.gis[max_ambiguous_fraction]} \
            --out-nuc {output.fasta} \
            --out-prot {output.faa} > {log} 2>&1
        """
