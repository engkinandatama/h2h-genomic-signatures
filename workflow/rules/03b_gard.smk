# Rule 03b: GARD — Recombination Detection
# MUST run AFTER codon alignment, BEFORE IQ-TREE and all selection analyses.
# If GARD detects significant breakpoints, downstream analyses are flagged.

rule hyphy_gard:
    """
    GARD: Genetic Algorithm Recombination Detection.
    Tests whether the codon alignment shows evidence of recombination.
    Recombination violates the phylogenetic assumptions of FEL/MEME/FUBAR/aBSREL/RELAX,
    so this step is a prerequisite validation for all selection analyses.

    Key outputs:
      - gard.json : full GARD output with breakpoint positions and improvement scores
      - gard_summary.txt : human-readable summary (recombination detected? Y/N, breakpoints)
    """
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_gard.json",
        summary=WORKDIR + "/04_selection/{virus_group}/{protein}_gard_summary.txt"
    log:
        WORKDIR + "/logs/hyphy_gard/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_gard/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    params:
        rate_classes=config.get("params", {}).get("gard", {}).get("rate_classes", 4),
        skip=lambda w: f"{w.virus_group}/{w.protein}" in (
            config.get("params", {}).get("gard", {}).get("skip_datasets", []) or [])
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy GARD (recombination detection) for {wildcards.virus_group} - {wildcards.protein}..." > {log}

        if [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Codon alignment kosong. Melewati GARD." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            python workflow/scripts/parse_gard.py \
                --json    {output.json} \
                --out     {output.summary} \
                --virus-group {wildcards.virus_group} \
                --protein {wildcards.protein} \
            --log {log} >> {log} 2>&1
            exit 0
        fi

        # Count sequences in alignment; GARD needs at least 4
        N_SEQ=$(grep -c '^>' {input.codon_aln} || echo 0)
        if [ "$N_SEQ" -lt 4 ]; then
            echo "SKIP: Hanya $N_SEQ sekuens, minimum 4 diperlukan untuk GARD." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            python workflow/scripts/parse_gard.py \
                --json    {output.json} \
                --out     {output.summary} \
                --virus-group {wildcards.virus_group} \
                --protein {wildcards.protein} \
            --log {log} >> {log} 2>&1
            exit 0
        fi

        if [ "{params.skip}" = "True" ]; then
            echo "SKIPPED by config (params.gard.skip_datasets): recombination" >> {log}
            echo "status is UNKNOWN for this alignment, not clean." >> {log}
            echo '{{"status": "SKIPPED"}}' > {output.json}
            # The summary is a declared output too. Writing only the JSON here made
            # Snakemake fail the rule for a missing file, so turning a dataset off
            # broke the run instead of shortening it.
            python workflow/scripts/parse_gard.py \
                --json    {output.json} \
                --out     {output.summary} \
                --virus-group {wildcards.virus_group} \
                --protein {wildcards.protein} \
            --log {log} >> {log} 2>&1
            exit 0
        fi

        # ENV must be a command-line argument; HyPhy does not read it from the process
        # environment, and its value is HBL source, so the trailing semicolon belongs to
        # HyPhy and must be quoted. Unquoted, bash split the line there: hyphy ran as its
        # own command with no redirect, and '>> log || fallback' became a separate no-op
        # that always succeeded, so GARD output went to the console and a crash never
        # wrote the FAILED sentinel. Writing --output-lf to /dev/null makes HyPhy try to
        # create '/dev/null.fit.bf', which fails, so use a real path.
        hyphy gard \
            --alignment {input.codon_aln} \
            --output {output.json} \
            --output-lf {output.json}.lf \
            --rate-classes {params.rate_classes} \
            CPU={threads} \
            ENV='TOLERATE_NUMERICAL_ERRORS=1;' >> {log} 2>&1 \
            || (echo '{{"status": "FAILED"}}' > {output.json} \
                && echo "ERROR: GARD did not complete; see above." >> {log})

        # Parse JSON to write human-readable summary
        python workflow/scripts/parse_gard.py \
            --json    {output.json} \
            --out     {output.summary} \
            --virus-group {wildcards.virus_group} \
            --protein {wildcards.protein} \
            --log {log} >> {log} 2>&1

        echo "GARD complete for {wildcards.virus_group} - {wildcards.protein}." >> {log}
        """
