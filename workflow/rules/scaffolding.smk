# Scaffolding rules: HapHiC + RagTag for species with HiC data.
#
# Strategy (per Aaryan's supervisor):
#   Pass 1: HapHiC on p_utg → manual Juicebox curation → curated p_utg backbone
#   Pass 2: RagTag (curated p_utg as ref, hap1+hap2 as query)
#           → HapHiC on RagTag output → manual curation → final hap1+hap2
#
# This file currently implements pass 1 only. Pass 2 rules (ragtag + haphic
# pass 2) will be added after manual curation reveals the pass 1 output works.
#
# Pass 1 input: initial p_utg (NOT decontaminated). Contaminant contigs are
# small and end up as unplaced scaffolds, dropped during curation. The pass 1
# output is used purely as a scaffolding reference for pass 2, where the
# DECONTAMINATED hap1/hap2 are the actual query.


HAPHIC_BIN = "methods/HapHiC/haphic"
HAPHIC_FILTER_BAM = "methods/HapHiC/utils/filter_bam.py"
HAPHIC_JUICER_POST = "methods/HapHiC/utils/juicer"


def hic_r1_for_species(wildcards):
    return hic_r1_files(wildcards.species)


def hic_r2_for_species(wildcards):
    return hic_r2_files(wildcards.species)

def pass1_putg_dir(species: str) -> str:
    """Directory holding the p_utg FASTA+GFA that pass-1 scaffolds.
    onepass_putg species (binata) scaffold the DECONTAMINATED p_utg so the
    contamination report covers the actual scaffolding input; all others use
    the initial p_utg (contaminants dropped during curation)."""
    stage = "decontaminated" if scaffolding_strategy(species) == "onepass_putg" else "initial"
    return f"results/{species}/assembly/{stage}/p_utg"


rule scaffold_haphic_pass1:
    """HapHiC pass 1: scaffold the initial p_utg into chromosome-scale scaffolds.
    Output is ready for manual Juicebox curation.

    User workflow after this rule completes:
      1. cd results/{species}/scaffolding/pass1/04.build/
      2. bash juicebox.sh  (already invoked by this rule, produces .hic + .assembly)
      3. Open out_JBAT.hic + out_JBAT.assembly in Juicebox GUI
      4. Edit, save as out_JBAT.review.assembly in the same dir
      5. Re-run snakemake; the scaffold_haphic_pass1_post_juicebox rule will pick up
    """
    input:
        ref=lambda wc: f"{pass1_putg_dir(wc.species)}/{wc.species}.fa",
        gfa=lambda wc: f"{pass1_putg_dir(wc.species)}/{wc.species}.gfa",
        hic_r1=hic_r1_for_species,
        hic_r2=hic_r2_for_species,
    output:
        scaffolds_fa="results/{species}/scaffolding/pass1/04.build/scaffolds.fa",
        raw_agp="results/{species}/scaffolding/pass1/04.build/scaffolds.raw.agp",
        liftover_agp="results/{species}/scaffolding/pass1/04.build/out_JBAT.liftover.agp",
        jbat_hic="results/{species}/scaffolding/pass1/04.build/out_JBAT.hic",
        jbat_assembly="results/{species}/scaffolding/pass1/04.build/out_JBAT.assembly",
        done=touch("results/{species}/scaffolding/pass1/pass1.done"),
    params:
        chr_num=lambda wc: chr_num_for_species(wc.species),
        workdir="results/{species}/scaffolding/pass1/work",
        outdir="results/{species}/scaffolding/pass1",
        ref_local=lambda wc: f"results/{wc.species}/scaffolding/pass1/work/{wc.species}.p_utg.fa",
        gfa_local=lambda wc: f"results/{wc.species}/scaffolding/pass1/work/{wc.species}.p_utg.gfa",
        haphic_bin=HAPHIC_BIN,
        haphic_filter_bam=HAPHIC_FILTER_BAM,
    threads: 40
    resources:
        mem_mb=120000,
    conda:
        "/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/haphic"
    log:
        "logs/scaffold_haphic_pass1/{species}.log"
    benchmark:
        "benchmarks/scaffold_haphic_pass1/{species}.tsv"
    shell:
        r"""
        set -euo pipefail

        # Absolute paths — we cd around a lot below, relative paths break.
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        WORKDIR="$REPO_ROOT/{params.workdir}"
        OUTDIR="$REPO_ROOT/{params.outdir}"

        mkdir -p "$WORKDIR" "$OUTDIR" "$(dirname $LOG)"

        # Copy ref + gfa to workdir (bwa index writes alongside the ref;
        # results/{wildcards.species}/assembly/initial/p_utg/ is a managed
        # dir we don't want to pollute with .bwt/.pac/etc index files)
        cp -L "$REPO_ROOT/{input.ref}" "$WORKDIR/{wildcards.species}.p_utg.fa"
        cp -L "$REPO_ROOT/{input.gfa}" "$WORKDIR/{wildcards.species}.p_utg.gfa"

        cd "$WORKDIR"
        REF="{wildcards.species}.p_utg.fa"
        GFA="{wildcards.species}.p_utg.gfa"

        # 1. Index reference
        echo "=== bwa index ===" >> "$LOG"
        bwa index "$REF" 2>> "$LOG"

        # 2. Align each HiC pair → per-library BAM
        R1_FILES=({input.hic_r1})
        R2_FILES=({input.hic_r2})

        if [ "${{#R1_FILES[@]}}" -ne "${{#R2_FILES[@]}}" ]; then
            echo "ERROR: R1/R2 count mismatch" >&2
            exit 1
        fi

        BAM_FILES=()
        for i in "${{!R1_FILES[@]}}"; do
            R1="${{R1_FILES[$i]}}"
            R2="${{R2_FILES[$i]}}"
            BAM="HIC_${{i}}.bam"
            BAM_FILES+=("$BAM")
            echo "=== bwa mem pair $((i+1))/${{#R1_FILES[@]}}: $R1 + $R2 ===" >> "$LOG"
            bwa mem -t {threads} -5SP "$REF" "$R1" "$R2" 2>> "$LOG" \
                | samblaster 2>> "$LOG" \
                | samtools view -@ {threads} -S -h -b -F 3340 -o "$BAM" 2>> "$LOG"
        done

        # 3. Merge per-library BAMs
        echo "=== samtools merge ===" >> "$LOG"
        samtools merge -@ {threads} -f HIC.bam "${{BAM_FILES[@]}}" 2>> "$LOG"

        # 4. Sort by read name
        echo "=== samtools sort -n ===" >> "$LOG"
        samtools sort -n -@ {threads} -o HiC.sorted.bam HIC.bam 2>> "$LOG"

        # 5. HapHiC filter_bam
        echo "=== HapHiC filter_bam ===" >> "$LOG"
        "$REPO_ROOT/{params.haphic_filter_bam}" HiC.sorted.bam 2 --NM 3 --threads {threads} 2>> "$LOG" \
            | samtools view -b -@ {threads} -o HiC.filtered.bam 2>> "$LOG"

        # 6. Run HapHiC pipeline (in outdir, not workdir, so 01.cluster etc end up clean)
        cd "$OUTDIR"
        ln -srf "$WORKDIR/HiC.filtered.bam" ./HiC.filtered.bam
        ln -srf "$WORKDIR/$REF" ./"$REF"
        ln -srf "$WORKDIR/$GFA" ./"$GFA"
        


        # HapHiC's pipeline mkdir's 01.cluster..04.build and ABORTS if any exist
        # (happens on any rerun). Clear the step dirs so reruns start clean.
        rm -rf 01.cluster 02.reassign 03.sort 04.build



        echo "=== haphic pipeline ===" >> "$LOG"
        "$REPO_ROOT/{params.haphic_bin}" pipeline \
            "$REF" HiC.filtered.bam {params.chr_num} \
            --threads {threads} --processes {threads} \
            --gfa "$GFA" \
            2>> "$LOG"

        # 7. juicebox.sh
        echo "=== juicebox.sh ===" >> "$LOG"
        cd "$OUTDIR/04.build/"
        rm -f "{wildcards.species}.p_utg.fa" #symlinked by juicebox.sh

        source /opt/share/software/scs/appStore/modules/init/profile.sh
        module load java/jdk-17.0.10
        export PATH=/usr/bin:$PATH #in case module load doesnt point to anything
        bash juicebox.sh 2>> "$LOG"

        echo "=== pass 1 algorithmic step complete ===" >> "$LOG"
        echo "Next: open $OUTDIR/04.build/out_JBAT.hic in Juicebox," >> "$LOG"
        echo "      edit, save as out_JBAT.review.assembly," >> "$LOG"
        echo "      re-run snakemake → post_juicebox rule will run." >> "$LOG"
        """       

rule scaffold_haphic_pass1_post_juicebox:
    """Apply user's manual Juicebox curation to produce FINAL.fa.

    Input includes out_JBAT.review.assembly, which the USER manually creates
    after curating in Juicebox. If this file is missing, snakemake fails with
    a clear 'missing input' error. Re-run snakemake after creating the file.
    """
    input:
        review_assembly="results/{species}/scaffolding/pass1/04.build/out_JBAT.review.assembly",
        liftover_agp="results/{species}/scaffolding/pass1/04.build/out_JBAT.liftover.agp",
        ref_fa="results/{species}/scaffolding/pass1/work/{species}.p_utg.fa",
    output:
        final_fa="results/{species}/scaffolding/pass1/05.post_juicebox/out_JBAT.FINAL.fa",
    params:
        outdir="results/{species}/scaffolding/pass1/05.post_juicebox",
        haphic_juicer=lambda wc: f"{os.getcwd()}/{HAPHIC_JUICER_POST}",
    threads: 2
    resources:
        mem_mb=8000,
    conda:
        "/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/haphic"
    log:
        "logs/scaffold_haphic_pass1_post_juicebox/{species}.log"
    shell:
        r"""
        set -euo pipefail
        # Absolute paths — we cd around a lot below, relative paths break.
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        OUTDIR="$REPO_ROOT/{params.outdir}"

        mkdir -p {params.outdir} logs/scaffold_haphic_pass1_post_juicebox

        # Sanity: review.assembly should differ from the auto-generated assembly.
        AUTO_ASSEMBLY="results/{wildcards.species}/scaffolding/pass1/04.build/out_JBAT.assembly"
        if [ -f "$AUTO_ASSEMBLY" ] && cmp -s "$AUTO_ASSEMBLY" "{input.review_assembly}"; then
            echo "WARNING: review.assembly is byte-identical to out_JBAT.assembly" >> "$LOG"
            echo "         did you actually curate in Juicebox? Continuing anyway." >> "$LOG"
        fi

        cd {params.outdir}
        {params.haphic_juicer} post \
            -o out_JBAT \
            "$REPO_ROOT/{input.review_assembly}" \
            "$REPO_ROOT/{input.liftover_agp}" \
            "$REPO_ROOT/{input.ref_fa}" \
            2>> "$LOG"
        echo "=== pass 1 post-juicebox complete ===" >> "$LOG"
        echo "FINAL.fa is at: {output.final_fa}" >> "$LOG"
        """


# ════════════════════════════════════════════════════════════════════════
# PASS 2: RagTag (curated pass-1 p_utg as ref; decontaminated hap1..hapN
# concatenated as query) → HapHiC pass 2 → manual curation → final scaffolds.
#
# Query hap count is species-dependent: diploids contribute hap1+hap2,
# autotetraploid scorpioides contributes hap1..hap4. Driven by
# decontamination_targets(species) so no hap count is hardcoded here.
#
# Pass 2 ref is the RagTag scaffold (not a gfa-backed assembly → no --gfa),
# and may take per-species HapHiC tuning via config['haphic_pass2_extra'].
# ════════════════════════════════════════════════════════════════════════

RAGTAG_ENV = "/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/ragtag"
HAPHIC_ENV = "/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/haphic"


def ragtag_query_haps(wildcards):
    """Decontaminated hap FASTAs to concatenate as the RagTag query.
    hap1..hapN for the species (p_utg targets are excluded — those species
    don't reach pass 2 since they have no HiC)."""
    return [
        f"results/{wildcards.species}/assembly/decontaminated/{tgt}/{wildcards.species}.fa"
        for tgt in decontamination_targets(wildcards.species)
        if tgt.startswith("hap")
    ]


def haphic_pass2_extra(species: str) -> str:
    """Extra HapHiC pipeline flags for pass 2, per species (default empty).
    Populate config['haphic_pass2_extra'][species] when plain output needs
    tuning, e.g. scorpioides:
    '--correct_nrounds 5 --region_len_ratio 0.5 --remove_allelic_links 4 --max_inflation 10'
    """
    return config.get("haphic_pass2_extra", {}).get(species, "")


rule scaffold_ragtag:
    """RagTag scaffold: order+orient decontaminated hap1..hapN against the
    curated pass-1 p_utg (FINAL.fa). Rescues small contigs that pure Hi-C
    scaffolding leaves unplaced."""
    input:
        ref="results/{species}/scaffolding/pass1/05.post_juicebox/out_JBAT.FINAL.fa",
        haps=ragtag_query_haps,
    output:
        scaffold_fa="results/{species}/scaffolding/pass2/ragtag/ragtag.scaffold.fasta",
    params:
        workdir="results/{species}/scaffolding/pass2/ragtag",
    threads: 40
    resources:
        mem_mb=120000,
    conda:
        RAGTAG_ENV
    log:
        "logs/scaffold_ragtag/{species}.log"
    benchmark:
        "benchmarks/scaffold_ragtag/{species}.tsv"
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        WORKDIR="$REPO_ROOT/{params.workdir}"
        mkdir -p "$WORKDIR" "$(dirname $LOG)"

        # Concatenate decontaminated haplotypes as the query. hifiasm names
        # contigs per-hap (h1tg/h2tg/...), so no header collisions across haps.
        cat {input.haps} > "$WORKDIR/all_haps.fa"

        cd "$WORKDIR"
        ragtag.py scaffold \
            "$REPO_ROOT/{input.ref}" \
            all_haps.fa \
            -t {threads} \
            -o ragtag_output \
            2>> "$LOG"

        cp ragtag_output/ragtag.scaffold.fasta ./ragtag.scaffold.fasta
        echo "=== ragtag scaffold complete ===" >> "$LOG"
        """


rule scaffold_haphic_pass2:
    """HapHiC pass 2: Hi-C scaffold the RagTag output. Same alignment/filtering
    as pass 1; ref is ragtag.scaffold.fasta, no --gfa, optional tuning params."""
    input:
        ref="results/{species}/scaffolding/pass2/ragtag/ragtag.scaffold.fasta",
        hic_r1=hic_r1_for_species,
        hic_r2=hic_r2_for_species,
    output:
        scaffolds_fa="results/{species}/scaffolding/pass2/04.build/scaffolds.fa",
        raw_agp="results/{species}/scaffolding/pass2/04.build/scaffolds.raw.agp",
        liftover_agp="results/{species}/scaffolding/pass2/04.build/out_JBAT.liftover.agp",
        jbat_hic="results/{species}/scaffolding/pass2/04.build/out_JBAT.hic",
        jbat_assembly="results/{species}/scaffolding/pass2/04.build/out_JBAT.assembly",
        done=touch("results/{species}/scaffolding/pass2/pass2.done"),
    params:
        chr_num=lambda wc: chr_num_for_species(wc.species),
        extra=lambda wc: haphic_pass2_extra(wc.species),
        workdir="results/{species}/scaffolding/pass2/work",
        outdir="results/{species}/scaffolding/pass2",
        haphic_bin=HAPHIC_BIN,
        haphic_filter_bam=HAPHIC_FILTER_BAM,
    threads: 40
    resources:
        mem_mb=120000,
    conda:
        HAPHIC_ENV
    log:
        "logs/scaffold_haphic_pass2/{species}.log"
    benchmark:
        "benchmarks/scaffold_haphic_pass2/{species}.tsv"
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        WORKDIR="$REPO_ROOT/{params.workdir}"
        OUTDIR="$REPO_ROOT/{params.outdir}"
        mkdir -p "$WORKDIR" "$OUTDIR" "$(dirname $LOG)"

        cp -L "$REPO_ROOT/{input.ref}" "$WORKDIR/{wildcards.species}.ragtag.fa"
        cd "$WORKDIR"
        REF="{wildcards.species}.ragtag.fa"

        echo "=== bwa index ===" >> "$LOG"
        bwa index "$REF" 2>> "$LOG"

        R1_FILES=({input.hic_r1})
        R2_FILES=({input.hic_r2})
        if [ "${{#R1_FILES[@]}}" -ne "${{#R2_FILES[@]}}" ]; then
            echo "ERROR: R1/R2 count mismatch" >&2
            exit 1
        fi

        BAM_FILES=()
        for i in "${{!R1_FILES[@]}}"; do
            R1="${{R1_FILES[$i]}}"; R2="${{R2_FILES[$i]}}"
            BAM="HIC_${{i}}.bam"; BAM_FILES+=("$BAM")
            echo "=== bwa mem pair $((i+1))/${{#R1_FILES[@]}} ===" >> "$LOG"
            bwa mem -t {threads} -5SP "$REF" "$R1" "$R2" 2>> "$LOG" \
                | samblaster 2>> "$LOG" \
                | samtools view -@ {threads} -S -h -b -F 3340 -o "$BAM" 2>> "$LOG"
        done

        echo "=== samtools merge ===" >> "$LOG"
        samtools merge -@ {threads} -f HIC.bam "${{BAM_FILES[@]}}" 2>> "$LOG"
        echo "=== samtools sort -n ===" >> "$LOG"
        samtools sort -n -@ {threads} -o HiC.sorted.bam HIC.bam 2>> "$LOG"
        echo "=== HapHiC filter_bam ===" >> "$LOG"
        "$REPO_ROOT/{params.haphic_filter_bam}" HiC.sorted.bam 2 --NM 3 --threads {threads} 2>> "$LOG" \
            | samtools view -b -@ {threads} -o HiC.filtered.bam 2>> "$LOG"

        cd "$OUTDIR"
        ln -srf "$WORKDIR/HiC.filtered.bam" ./HiC.filtered.bam
        ln -srf "$WORKDIR/$REF" ./"$REF"


        # HapHiC's pipeline mkdir's 01.cluster..04.build and ABORTS if any exist
        # (happens on any rerun). Clear the step dirs so reruns start clean.
        rm -rf 01.cluster 02.reassign 03.sort 04.build

        



        echo "=== haphic pipeline (pass2, extra: {params.extra}) ===" >> "$LOG"
        "$REPO_ROOT/{params.haphic_bin}" pipeline \
            "$REF" HiC.filtered.bam {params.chr_num} \
            --threads {threads} --processes {threads} \
            {params.extra} \
            2>> "$LOG"

        echo "=== juicebox.sh ===" >> "$LOG"
        cd "$OUTDIR/04.build/"
        rm -f "{wildcards.species}.ragtag.fa"  # symlinked by juicebox.sh
        export PATH=/usr/bin:$PATH  # system java 17 for juicer_tools
        bash juicebox.sh 2>> "$LOG"

        echo "=== pass 2 algorithmic step complete ===" >> "$LOG"
        echo "Next: curate $OUTDIR/04.build/out_JBAT.hic in Juicebox," >> "$LOG"
        echo "      save out_JBAT.review.assembly, rerun for post_juicebox." >> "$LOG"
        """


rule scaffold_haphic_pass2_post_juicebox:
    """Apply manual Juicebox curation of pass 2 → final scaffolds FASTA.
    Requires user-created out_JBAT.review.assembly (snakemake fails clean if absent)."""
    input:
        review_assembly="results/{species}/scaffolding/pass2/04.build/out_JBAT.review.assembly",
        liftover_agp="results/{species}/scaffolding/pass2/04.build/out_JBAT.liftover.agp",
        ref_fa="results/{species}/scaffolding/pass2/work/{species}.ragtag.fa",
    output:
        final_fa="results/{species}/scaffolding/pass2/05.post_juicebox/out_JBAT.FINAL.fa",
    params:
        outdir="results/{species}/scaffolding/pass2/05.post_juicebox",
        haphic_juicer=HAPHIC_JUICER_POST,
    threads: 2
    resources:
        mem_mb=8000,
    conda:
        HAPHIC_ENV
    log:
        "logs/scaffold_haphic_pass2_post_juicebox/{species}.log"
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        OUTDIR="$REPO_ROOT/{params.outdir}"
        mkdir -p "$OUTDIR" "$(dirname $LOG)"

        AUTO_ASSEMBLY="results/{wildcards.species}/scaffolding/pass2/04.build/out_JBAT.assembly"
        if [ -f "$AUTO_ASSEMBLY" ] && cmp -s "$AUTO_ASSEMBLY" "{input.review_assembly}"; then
            echo "WARNING: review.assembly is byte-identical to out_JBAT.assembly" >> "$LOG"
            echo "         did you actually curate in Juicebox? Continuing anyway." >> "$LOG"
        fi

        cd "$OUTDIR"
        "$REPO_ROOT/{params.haphic_juicer}" post \
            -o out_JBAT \
            "$REPO_ROOT/{input.review_assembly}" \
            "$REPO_ROOT/{input.liftover_agp}" \
            "$REPO_ROOT/{input.ref_fa}" \
            2>> "$LOG"
        echo "=== pass 2 post-juicebox complete ===" >> "$LOG"
        echo "FINAL.fa is at: {output.final_fa}" >> "$LOG"
        """


# ════════════════════════════════════════════════════════════════════════
# ONE-PASS (gfa-phased): HapHiC on cat(decon hap1..hapN) using the per-hap
# decontaminated GFAs as a comma-separated --gfa list. HapHiC uses the graphs
# to keep haplotypes in separate groups (removes inter-haplotype Hi-C links),
# resolving all 2n scaffolds. For well-phased genomes where Hi-C alone cannot
# split low-divergence haplotypes but the assembly graph can (binata, paradoxa).
#
# Selected per species via scaffolding_strategy(species) == "onepass_gfa".
# ════════════════════════════════════════════════════════════════════════


def decon_hap_fastas_input(wildcards):
    return decon_hap_fastas(wildcards.species)


def decon_hap_gfas_input(wildcards):
    return decon_hap_gfas(wildcards.species)


def onepass_gfa_arg(wildcards):
    """Comma-separated decontaminated hap GFA paths for HapHiC --gfa.
    HapHiC keeps each gfa's contigs in separate haplotype groups."""
    return ",".join(decon_hap_gfas(wildcards.species))


def haphic_onepass_extra(species: str) -> str:
    """Extra HapHiC flags for one-pass, per species (default empty).
    Populate config['haphic_onepass_extra'][species] if tuning is needed."""
    return config.get("haphic_onepass_extra", {}).get(species, "")


rule scaffold_haphic_onepass:
    """One-pass HapHiC: scaffold cat(decon hap1..hapN) with per-hap decon GFAs
    (comma-separated --gfa) for haplotype phasing. Output ready for Juicebox
    curation (same downstream workflow as pass1)."""
    input:
        haps=decon_hap_fastas_input,
        gfas=decon_hap_gfas_input,
        hic_r1=hic_r1_for_species,
        hic_r2=hic_r2_for_species,
    output:
        scaffolds_fa="results/{species}/scaffolding/onepass/04.build/scaffolds.fa",
        raw_agp="results/{species}/scaffolding/onepass/04.build/scaffolds.raw.agp",
        liftover_agp="results/{species}/scaffolding/onepass/04.build/out_JBAT.liftover.agp",
        jbat_hic="results/{species}/scaffolding/onepass/04.build/out_JBAT.hic",
        jbat_assembly="results/{species}/scaffolding/onepass/04.build/out_JBAT.assembly",
        done=touch("results/{species}/scaffolding/onepass/onepass.done"),
    params:
        chr_num=lambda wc: chr_num_for_species(wc.species),
        gfa_arg=onepass_gfa_arg,
        extra=lambda wc: haphic_onepass_extra(wc.species),
        workdir="results/{species}/scaffolding/onepass/work",
        outdir="results/{species}/scaffolding/onepass",
        haphic_bin=HAPHIC_BIN,
        haphic_filter_bam=HAPHIC_FILTER_BAM,
    threads: 40
    resources:
        mem_mb=120000,
    conda:
        HAPHIC_ENV
    log:
        "logs/scaffold_haphic_onepass/{species}.log"
    benchmark:
        "benchmarks/scaffold_haphic_onepass/{species}.tsv"
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        WORKDIR="$REPO_ROOT/{params.workdir}"
        OUTDIR="$REPO_ROOT/{params.outdir}"
        mkdir -p "$WORKDIR" "$OUTDIR" "$(dirname $LOG)"

        # Concatenate decontaminated hap FASTAs → allhaps reference. hifiasm
        # names contigs per-hap (h1tg/h2tg/...), so no header collisions.
        cat {input.haps} > "$WORKDIR/{wildcards.species}.allhaps.fa"

        # GFAs stay SEPARATE — passed as comma-separated list to --gfa so HapHiC
        # treats each as its own haplotype. Build absolute paths.
        GFA_ARG=""
        for g in {input.gfas}; do
            abs="$REPO_ROOT/$g"
            if [ -z "$GFA_ARG" ]; then GFA_ARG="$abs"; else GFA_ARG="$GFA_ARG,$abs"; fi
        done
        echo "=== gfa arg: $GFA_ARG ===" >> "$LOG"

        cd "$WORKDIR"
        REF="{wildcards.species}.allhaps.fa"

        echo "=== bwa index ===" >> "$LOG"
        bwa index "$REF" 2>> "$LOG"

        R1_FILES=({input.hic_r1})
        R2_FILES=({input.hic_r2})
        if [ "${{#R1_FILES[@]}}" -ne "${{#R2_FILES[@]}}" ]; then
            echo "ERROR: R1/R2 count mismatch" >&2
            exit 1
        fi

        BAM_FILES=()
        for i in "${{!R1_FILES[@]}}"; do
            R1="${{R1_FILES[$i]}}"; R2="${{R2_FILES[$i]}}"
            BAM="HIC_${{i}}.bam"; BAM_FILES+=("$BAM")
            echo "=== bwa mem pair $((i+1))/${{#R1_FILES[@]}} ===" >> "$LOG"
            bwa mem -t {threads} -5SP "$REF" "$R1" "$R2" 2>> "$LOG" \
                | samblaster 2>> "$LOG" \
                | samtools view -@ {threads} -S -h -b -F 3340 -o "$BAM" 2>> "$LOG"
        done

        echo "=== samtools merge ===" >> "$LOG"
        samtools merge -@ {threads} -f HIC.bam "${{BAM_FILES[@]}}" 2>> "$LOG"
        echo "=== samtools sort -n ===" >> "$LOG"
        samtools sort -n -@ {threads} -o HiC.sorted.bam HIC.bam 2>> "$LOG"
        echo "=== HapHiC filter_bam ===" >> "$LOG"
        "$REPO_ROOT/{params.haphic_filter_bam}" HiC.sorted.bam 2 --NM 3 --threads {threads} 2>> "$LOG" \
            | samtools view -b -@ {threads} -o HiC.filtered.bam 2>> "$LOG"

        cd "$OUTDIR"
        ln -srf "$WORKDIR/HiC.filtered.bam" ./HiC.filtered.bam
        ln -srf "$WORKDIR/$REF" ./"$REF"


        # HapHiC's pipeline mkdir's 01.cluster..04.build and ABORTS if any exist
        # (happens on any rerun). Clear the step dirs so reruns start clean.
        rm -rf 01.cluster 02.reassign 03.sort 04.build

        echo "=== haphic pipeline (onepass, gfa-phased, extra: {params.extra}) ===" >> "$LOG"
        "$REPO_ROOT/{params.haphic_bin}" pipeline \
            "$REF" HiC.filtered.bam {params.chr_num} \
            --threads {threads} --processes {threads} \
            --gfa "$GFA_ARG" \
            {params.extra} \
            2>> "$LOG"

        echo "=== juicebox.sh ===" >> "$LOG"
        cd "$OUTDIR/04.build/"
        rm -f "{wildcards.species}.allhaps.fa"  # symlinked by juicebox.sh
        export PATH=/usr/bin:$PATH  # system java 17 for juicer_tools
        bash juicebox.sh 2>> "$LOG"

        echo "=== onepass algorithmic step complete ===" >> "$LOG"
        echo "Next: curate $OUTDIR/04.build/out_JBAT.hic in Juicebox," >> "$LOG"
        echo "      save out_JBAT.review.assembly, rerun for post_juicebox." >> "$LOG"
        """


rule scaffold_haphic_onepass_post_juicebox:
    """Apply manual Juicebox curation of one-pass → final scaffolds FASTA.
    Requires user-created out_JBAT.review.assembly (snakemake fails clean if absent)."""
    input:
        review_assembly="results/{species}/scaffolding/onepass/04.build/out_JBAT.review.assembly",
        liftover_agp="results/{species}/scaffolding/onepass/04.build/out_JBAT.liftover.agp",
        ref_fa="results/{species}/scaffolding/onepass/work/{species}.allhaps.fa",
    output:
        final_fa="results/{species}/scaffolding/onepass/05.post_juicebox/out_JBAT.FINAL.fa",
    params:
        outdir="results/{species}/scaffolding/onepass/05.post_juicebox",
        haphic_juicer=HAPHIC_JUICER_POST,
    threads: 2
    resources:
        mem_mb=8000,
    conda:
        HAPHIC_ENV
    log:
        "logs/scaffold_haphic_onepass_post_juicebox/{species}.log"
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        OUTDIR="$REPO_ROOT/{params.outdir}"
        mkdir -p "$OUTDIR" "$(dirname $LOG)"

        AUTO_ASSEMBLY="results/{wildcards.species}/scaffolding/onepass/04.build/out_JBAT.assembly"
        if [ -f "$AUTO_ASSEMBLY" ] && cmp -s "$AUTO_ASSEMBLY" "{input.review_assembly}"; then
            echo "WARNING: review.assembly is byte-identical to out_JBAT.assembly" >> "$LOG"
        fi

        cd "$OUTDIR"
        "$REPO_ROOT/{params.haphic_juicer}" post \
            -o out_JBAT \
            "$REPO_ROOT/{input.review_assembly}" \
            "$REPO_ROOT/{input.liftover_agp}" \
            "$REPO_ROOT/{input.ref_fa}" \
            2>> "$LOG"
        echo "=== onepass post-juicebox complete ===" >> "$LOG"
        echo "FINAL.fa is at: {output.final_fa}" >> "$LOG"
        """
