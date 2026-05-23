# Rule 5: Structural Mapping 3D (AlphaFold + Py3Dmol)
rule map_structure:
    input:
        selection_results=WORKDIR + "/04_selection/{virus}/{protein}_selection_results.txt"
    output:
        mapped_pdb=WORKDIR + "/05_structure/{virus}/{protein}_mapped.pdb",
        html_view=WORKDIR + "/05_structure/{virus}/{protein}_3d_view.html"
    log:
        WORKDIR + "/logs/map_structure/{virus}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/map_structure/{virus}_{protein}.tsv"
    params:
        uniprot_id=lambda wildcards: config["viruses"][wildcards.virus].get("uniprot_id", "UNKNOWN")
    conda:
        "../envs/structure.yaml"
    shell:
        """
        echo "Starting Structural Mapping for UniProt {params.uniprot_id}..." > {log}
        python workflow/scripts/map_mutations_3d.py \
            --uniprot {params.uniprot_id} \
            --selection {input.selection_results} \
            --out_pdb {output.mapped_pdb} \
            --out_html {output.html_view} >> {log} 2>&1
        echo "Structural mapping complete." >> {log}
        """
