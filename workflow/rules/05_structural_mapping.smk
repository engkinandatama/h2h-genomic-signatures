# Rule 5: Structural Mapping 3D (AlphaFold + Py3Dmol) - modified for virus_groups
rule map_structure:
    input:
        selection_results=WORKDIR + "/04_selection/{virus_group}/{protein}_all_sites.tsv",
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta"
    output:
        mapped_pdb=WORKDIR + "/05_structure/{virus_group}/{protein}_mapped.pdb",
        html_view=WORKDIR + "/05_structure/{virus_group}/{protein}_3d_view.html",
        mapping=WORKDIR + "/05_structure/{virus_group}/{protein}_coordinate_map.tsv"
    log:
        WORKDIR + "/logs/map_structure/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/map_structure/{virus_group}_{protein}.tsv"
    params:
        uniprot_id=lambda wildcards: (
            config["virus_groups"][wildcards.virus_group]
            .get("uniprot_ids", {})
            .get(wildcards.protein,
                 config["virus_groups"][wildcards.virus_group].get("uniprot_id", "UNKNOWN"))
        )
    conda:
        "../envs/structure.yaml"
    shell:
        """
        echo "Starting Structural Mapping for UniProt {params.uniprot_id}..." > {log}
        
        # Skip jika tidak ada UniProt ID yang dikonfigurasi
        if [ "{params.uniprot_id}" = "UNKNOWN" ]; then
            echo "SKIP: UniProt ID tidak dikonfigurasi untuk {wildcards.virus_group}/{wildcards.protein}." >> {log}
            mkdir -p $(dirname {output.mapped_pdb})
            touch {output.mapped_pdb} {output.html_view} {output.mapping}
            exit 0
        fi
        
        # Skip jika file selection_results tidak ada sama sekali
        if [ ! -f "{input.selection_results}" ]; then
            echo "SKIP: File selection results tidak ditemukan. Membuat output kosong." >> {log}
            mkdir -p $(dirname {output.mapped_pdb})
            touch {output.mapped_pdb} {output.html_view} {output.mapping}
            exit 0
        fi
        
        python workflow/scripts/map_mutations_3d.py \
            --uniprot {params.uniprot_id} \
            --selection {input.selection_results} \
            --alignment {input.codon_aln} \
            --virus-group {wildcards.virus_group} \
            --protein {wildcards.protein} \
            --out_mapping {output.mapping} \
            --out_pdb {output.mapped_pdb} \
            --out_html {output.html_view} >> {log} 2>&1
        echo "Structural mapping complete." >> {log}
        """
