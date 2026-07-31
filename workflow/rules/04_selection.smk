# Rule 4: Positive Selection Analysis (modified for virus_groups)

rule label_tree:
    input:
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_tree.treefile"
    output:
        labeled_tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    log:
        WORKDIR + "/logs/label_tree/{virus_group}_{protein}.log"
    params:
        bootstrap_threshold=config.get("params", {}).get("bootstrap_threshold", 0)
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running tree labeling for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ]; then
            echo "SKIP: File treefile kosong. Membuat file labeled_tree kosong." >> {log}
            touch {output.labeled_tree}
            exit 0
        fi
        python workflow/scripts/label_tree.py \
            --tree {input.tree} \
            --out  {output.labeled_tree} \
            --bootstrap-threshold {params.bootstrap_threshold} >> {log} 2>&1
        """

rule hyphy_meme:
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_meme.json"
    log:
        WORKDIR + "/logs/hyphy_meme/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_meme/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy MEME for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong. Membuat output JSON kosong." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            exit 0
        fi
        
        hyphy meme \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --output {output.json} \
            CPU={threads} \
            --branches Foreground >> {log} 2>&1 || (echo '{{"status": "FAILED"}}' > {output.json} && echo "ERROR: HyPhy MEME failed, created empty JSON" >> {log})
        """

rule hyphy_fel:
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_fel.json"
    log:
        WORKDIR + "/logs/hyphy_fel/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_fel/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy FEL for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong. Membuat output JSON kosong." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            exit 0
        fi
        
        hyphy fel \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --output {output.json} \
            CPU={threads} \
            --branches Foreground >> {log} 2>&1 || (echo '{{"status": "FAILED"}}' > {output.json} && echo "ERROR: HyPhy FEL failed, created empty JSON" >> {log})
        """

rule hyphy_fubar:
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_tree.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_fubar.json"
    log:
        WORKDIR + "/logs/hyphy_fubar/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_fubar/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy FUBAR for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong. Membuat output JSON kosong." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            exit 0
        fi
        
        # Remove any stale grid cache: HyPhy reuses *.FUBAR.cache without checking
        # that it matches the current alignment, which both skips the computation
        # and crashes when alignment dimensions have changed.
        rm -f {input.codon_aln}.FUBAR.cache

        hyphy fubar \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --output {output.json} \
            CPU={threads} >> {log} 2>&1 || (echo '{{"status": "FAILED"}}' > {output.json} && echo "ERROR: HyPhy FUBAR failed, created empty JSON" >> {log})
        """

rule aggregate_selection:
    input:
        meme=WORKDIR + "/04_selection/{virus_group}/{protein}_meme.json",
        fel=WORKDIR + "/04_selection/{virus_group}/{protein}_fel.json",
        fubar=WORKDIR + "/04_selection/{virus_group}/{protein}_fubar.json"
    output:
        results=WORKDIR + "/04_selection/{virus_group}/{protein}_selection_results.txt",
        meme_txt=WORKDIR + "/04_selection/{virus_group}/{protein}_meme_results.txt",
        fel_txt=WORKDIR + "/04_selection/{virus_group}/{protein}_fel_results.txt",
        fubar_txt=WORKDIR + "/04_selection/{virus_group}/{protein}_fubar_results.txt",
        complete=WORKDIR + "/04_selection/{virus_group}/{protein}_selection_complete.tsv"
    log:
        WORKDIR + "/logs/aggregate_selection/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/aggregate_selection/{virus_group}_{protein}.tsv"
    params:
        pval=config.get("params", {}).get("hyphy", {}).get("pvalue_threshold", 0.05),
        fubar_pp=config.get("params", {}).get("hyphy", {}).get("fubar_pp_threshold", 0.90),
        min_methods=config.get("params", {}).get("hyphy", {}).get("min_methods", 2)
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Aggregating HyPhy results for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        python workflow/scripts/parse_hyphy.py \
            --meme {input.meme} \
            --fel {input.fel} \
            --fubar {input.fubar} \
            --pvalue {params.pval} \
            --fubar-pp {params.fubar_pp} \
            --min-methods {params.min_methods} \
            --out-consensus {output.results} \
            --out-meme {output.meme_txt} \
            --out-fel {output.fel_txt} \
            --out-fubar {output.fubar_txt} \
            --out-complete {output.complete} >> {log} 2>&1
        """


# =============================================================================
# Branch-Specific Selection Analysis (Host-Range Change Signatures)
# =============================================================================

rule hyphy_busted:
    """
    BUSTED: Branch-Site Unrestricted Test for Episodic Diversification.
    Gene-level gate test — has this gene experienced positive selection
    anywhere on the Foreground (H2H) branches?
    Faster to run than aBSREL; use as a pre-filter.
    """
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_busted.json"
    log:
        WORKDIR + "/logs/hyphy_busted/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_busted/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy BUSTED for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong. Membuat output JSON kosong." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            exit 0
        fi

        hyphy busted \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --output {output.json} \
            --branches Foreground \
            CPU={threads} >> {log} 2>&1 \
            || (echo '{{"status": "FAILED"}}' > {output.json} && echo "ERROR: HyPhy BUSTED failed, created empty JSON" >> {log})
        """


rule hyphy_absrel:
    """
    aBSREL: Adaptive Branch-Site Random Effects Likelihood.
    Branch-level resolution — identifies WHICH specific branches
    have experienced episodic diversifying selection.
    More granular than BUSTED, heavier computation.
    """
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_absrel.json"
    log:
        WORKDIR + "/logs/hyphy_absrel/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_absrel/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy aBSREL for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong. Membuat output JSON kosong." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            exit 0
        fi

        hyphy absrel \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --output {output.json} \
            --branches Foreground \
            CPU={threads} >> {log} 2>&1 \
            || (echo '{{"status": "FAILED"}}' > {output.json} && echo "ERROR: HyPhy aBSREL failed, created empty JSON" >> {log})
        """


rule hyphy_relax:
    """
    RELAX: Tests whether selection pressure is relaxed (k < 1)
    or intensified (k > 1) on Foreground (H2H) branches compared
    to Reference (Reservoir/Spillover) branches.
    Key question: Did host-range change accompany intensified
    or relaxed purifying/positive selection?
    """
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_relax.json"
    params:
        branch_contrast=lambda w: has_branch_contrast(w.virus_group)
    log:
        WORKDIR + "/logs/hyphy_relax/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_relax/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy RELAX for {wildcards.virus_group} - {wildcards.protein}..." > {log}

        if [ "{params.branch_contrast}" != "True" ]; then
            echo "NOT APPLICABLE: {wildcards.virus_group} has no reservoir group, so there" >> {log}
            echo "is no reference branch set. RELAX is not run." >> {log}
            echo '{{"status": "NOT_APPLICABLE"}}' > {output.json}
            exit 0
        fi

        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "ERROR: Input tree atau alignment kosong." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            exit 0
        fi

        hyphy relax \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --output {output.json} \
            --test Foreground \
            --reference Reference \
            CPU={threads} >> {log} 2>&1 \
            || (echo '{{"status": "FAILED"}}' > {output.json} && echo "ERROR: HyPhy RELAX failed, created empty JSON" >> {log})
        """


rule aggregate_branch_selection:
    """
    Parse BUSTED + aBSREL + RELAX results into a single per-protein TSV summary.
    This feeds into the final convergence statistics rule.
    """
    input:
        busted=WORKDIR + "/04_selection/{virus_group}/{protein}_busted.json",
        absrel=WORKDIR + "/04_selection/{virus_group}/{protein}_absrel.json",
        relax=WORKDIR + "/04_selection/{virus_group}/{protein}_relax.json"
    output:
        summary=WORKDIR + "/04_selection/{virus_group}/{protein}_branch_selection.tsv"
    log:
        WORKDIR + "/logs/aggregate_branch_selection/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/aggregate_branch_selection/{virus_group}_{protein}.tsv"
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Aggregating branch-specific selection for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        python workflow/scripts/parse_branch_selection.py \
            --busted {input.busted} \
            --absrel {input.absrel} \
            --relax  {input.relax} \
            --out    {output.summary} \
            --virus-group {wildcards.virus_group} \
            --protein {wildcards.protein} >> {log} 2>&1
        """


# =============================================================================
# Additional Site-Level Analyses
# =============================================================================

rule hyphy_slac:
    """
    SLAC: Single Likelihood Ancestor Counting.
    Fast, conservative dN/dS estimator using parsimony-based counting.
    Serves as a methodological sanity check alongside FEL/FUBAR.
    Reports per-site dS, dN, dN-dS, and significance of deviation from neutrality.
    """
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_slac.json"
    log:
        WORKDIR + "/logs/hyphy_slac/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_slac/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy SLAC for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            exit 0
        fi

        hyphy slac \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --branches All \
            --output {output.json} \
            CPU={threads} >> {log} 2>&1 \
            || (echo '{{"status": "FAILED"}}' > {output.json} && echo "ERROR: SLAC did not complete" >> {log})
        """


rule hyphy_contrast_fel:
    """
    Contrast-FEL: Directly compares dN/dS between H2H (Foreground) branches
    and Reservoir/Spillover (Background) branches at each individual codon site.

    This is the most direct test for host-range-specific selection:
      - Significant site + beta_H2H > beta_Reservoir => site under STRONGER
        positive selection in H2H than in Reservoir
      - This identifies the specific amino acid positions that may drive
        human-to-human transmissibility

    Requires labeled tree (Foreground = H2H branches, Background = Reservoir).
    """
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_contrast_fel.json"
    params:
        branch_contrast=lambda w: has_branch_contrast(w.virus_group)
    log:
        WORKDIR + "/logs/hyphy_contrast_fel/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_contrast_fel/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy Contrast-FEL for {wildcards.virus_group} - {wildcards.protein}..." > {log}

        if [ "{params.branch_contrast}" != "True" ]; then
            echo "NOT APPLICABLE: {wildcards.virus_group} has no reservoir group, so no" >> {log}
            echo "foreground/background partition exists. Contrast-FEL is not run." >> {log}
            echo '{{"status": "NOT_APPLICABLE"}}' > {output.json}
            exit 0
        fi

        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "ERROR: Input tree atau alignment kosong." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            exit 0
        fi

        hyphy contrast-fel \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --branch-set Foreground \
            --output {output.json} \
            CPU={threads} \
            ENV='TOLERATE_NUMERICAL_ERRORS=1;' >> {log} 2>&1 \
            || (echo '{{"status": "FAILED"}}' > {output.json} \
                && echo "ERROR: Contrast-FEL did not complete." >> {log})
        """


rule hyphy_prime:
    """
    PRIME: PRoperty Informed Models of Evolution.
    Tests whether amino acid substitutions at positively selected sites
    are biased toward specific physicochemical properties:
      - Volume (size of the amino acid)
      - Polarity
      - Charge
      - Hydrophobicity
      - Composition

    Answers: "Are the selected mutations changing the CHARGE or HYDROPHOBICITY
    of the protein surface?" — critical for receptor binding interpretation.
    """
    input:
        codon_aln=WORKDIR + "/02_aligned/{virus_group}/{protein}_codon_aligned.fasta",
        tree=WORKDIR + "/03_trees/{virus_group}/{protein}_labeled.treefile"
    output:
        json=WORKDIR + "/04_selection/{virus_group}/{protein}_prime.json"
    log:
        WORKDIR + "/logs/hyphy_prime/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/hyphy_prime/{virus_group}_{protein}.tsv"
    threads:
        config.get("resources", {}).get("threads", {}).get("hyphy", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", {}).get("hyphy", 8000),
        hyphy_jobs=1
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Running HyPhy PRIME for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        if [ ! -s "{input.tree}" ] || [ ! -s "{input.codon_aln}" ]; then
            echo "SKIP: Input tree atau alignment kosong." >> {log}
            echo '{{"status": "FAILED"}}' > {output.json}
            exit 0
        fi

        hyphy prime \
            --alignment {input.codon_aln} \
            --tree {input.tree} \
            --branches Foreground \
            --output {output.json} \
            CPU={threads} >> {log} 2>&1 \
            || (echo '{{"status": "FAILED"}}' > {output.json} && echo "ERROR: PRIME did not complete" >> {log})
        """


rule aggregate_all_selection:
    """
    Merge ALL site-level selection results into one comprehensive TSV:
      - MEME episodic p-value
      - FEL pervasive p-value
      - FUBAR posterior probability
      - SLAC dN, dS, dN-dS, p-value
      - Contrast-FEL: beta_H2H, beta_Reservoir, differential p-value
      - PRIME: property biases (volume, polarity, charge, hydrophobicity, composition)
      - Consensus call (how many methods agree)
    This is the master per-site table used for 3D visualization and convergence.
    """
    input:
        meme=WORKDIR + "/04_selection/{virus_group}/{protein}_meme.json",
        fel=WORKDIR + "/04_selection/{virus_group}/{protein}_fel.json",
        fubar=WORKDIR + "/04_selection/{virus_group}/{protein}_fubar.json",
        slac=WORKDIR + "/04_selection/{virus_group}/{protein}_slac.json",
        contrast_fel=WORKDIR + "/04_selection/{virus_group}/{protein}_contrast_fel.json",
        prime=WORKDIR + "/04_selection/{virus_group}/{protein}_prime.json",
        gard=WORKDIR + "/04_selection/{virus_group}/{protein}_gard_summary.txt"
    output:
        full_table=WORKDIR + "/04_selection/{virus_group}/{protein}_all_sites.tsv"
    log:
        WORKDIR + "/logs/aggregate_all_selection/{virus_group}_{protein}.log"
    benchmark:
        WORKDIR + "/benchmarks/aggregate_all_selection/{virus_group}_{protein}.tsv"
    params:
        pval=config.get("params", {}).get("hyphy", {}).get("pvalue_threshold", 0.05),
        min_methods=config.get("params", {}).get("hyphy", {}).get("min_methods", 2),
        fubar_pp=config.get("params", {}).get("hyphy", {}).get("fubar_pp_threshold", 0.90),
        cfdr=config.get("params", {}).get("hyphy", {}).get("contrast_fel_fdr", 0.20)
    conda:
        "../envs/selection.yaml"
    shell:
        """
        echo "Aggregating ALL site-level selection results for {wildcards.virus_group} - {wildcards.protein}..." > {log}
        python workflow/scripts/parse_additional_selection.py \
            --meme         {input.meme} \
            --fel          {input.fel} \
            --fubar        {input.fubar} \
            --slac         {input.slac} \
            --contrast-fel {input.contrast_fel} \
            --prime        {input.prime} \
            --gard         {input.gard} \
            --pvalue       {params.pval} \
            --fubar-pp     {params.fubar_pp} \
            --contrast-fdr {params.cfdr} \
            --min-methods  {params.min_methods} \
            --virus-group  {wildcards.virus_group} \
            --protein      {wildcards.protein} \
            --out          {output.full_table} >> {log} 2>&1
        """
