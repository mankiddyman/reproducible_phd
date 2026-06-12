# Structural annotation: per-tool gene prediction on frozen chromosome-scale
# assemblies, for the methods comparison (Helixer / Tiberius / ANNEVO / BRAKER).
#
# GPU tools run via Singularity on the GPU node. Because gpu01 is currently
# DOWN in SLURM, these rules are run by a LOCAL snakemake invocation ON gpu01
# (snakemake --cores N --resources gpu=2), not via the slurm profile. They
# declare `resources: gpu=1` so that if gpu01 is restored to SLURM the same
# rules submit cleanly with a gres mapping.
#
# All tools annotate the frozen {species}_chr.fa (chromosomes only, no debris).
# Genome split per-chromosome for parallel GPU annotation; per-chr GFF3s merged
# + standardized with gffread.

annotation_df = pd.read_csv(config["annotation_table"]).fillna("").set_index("species_id", drop=False)

# External (collaborator-assembled) genomes — enter the pipeline at the frozen
# stage; no hifiasm/decon/scaffold. Driven by config/external_genomes.csv, kept
# OUT of species.csv so the assembly DAG never sees them.
external_genomes_df = pd.read_csv(config["external_genomes_table"]).fillna("").set_index("species_id", drop=False)

def keep_n_scaffolds(species: str) -> int:
    """Number of chromosome-scale scaffolds to retain from a collapsed external
    genome = chr_number // ploidy (collapsed → one scaffold per chromosome)."""
    row = external_genomes_df.loc[species]
    return int(row["chr_number"]) // int(row["ploidy"])


ANNOT_TIBERIUS_SIF = config["annotation"]["tiberius_sif"]
ANNOT_HELIXER_SIF  = config["annotation"]["helixer_sif"]
GFFREAD            = config["annotation"]["gffread"]
AGAT_ENV           = config["annotation"]["agat_env"]
ANNEVO_ENV         = config["annotation"]["annevo_env"]
ANNEVO_REPO        = config["annotation"]["annevo_repo"]
ANNEVO_MODEL       = config["annotation"]["annevo_model"]


def frozen_chr_fasta(species: str) -> str:
    """Frozen chromosome-scale FASTA ({species}_chr.fa), stage from annotation.csv."""
    stage = str(annotation_df.loc[species, "frozen_stage"]).strip()
    if not stage:
        raise ValueError(f"frozen_stage not set for {species} in annotation.csv.")
    return f"results/{species}/assembly_final/{stage}/{species}_chr.fa"


# Pick least-used GPU (by used memory) and export CUDA_VISIBLE_DEVICES.
PICK_GPU = r'''
        GPU_ID=$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits \
                 | sort -t, -k2 -n | head -1 | cut -d, -f1 | tr -d ' ')
        export CUDA_VISIBLE_DEVICES="$GPU_ID"
        echo "Using GPU $GPU_ID" >> "$LOG"
'''


rule annotate_split:
    """Split frozen {species}_chr.fa into one FASTA per chromosome."""
    input:
        chr_fa=lambda wc: frozen_chr_fasta(wc.species),
    output:
        split_done=touch("results/{species}/annotation/_split/.split.done"),
    params:
        split_dir="results/{species}/annotation/_split",
    log:
        "logs/annotate_split/{species}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.split_dir}" "$(dirname {log})"
        awk '/^>/ {{ name=substr($1,2); file="{params.split_dir}/" name ".fa" }}
                  {{ print > file }}' "{input.chr_fa}" 2> {log}
        echo "split $(ls {params.split_dir}/*.fa | wc -l) chromosomes" >> {log}
        """


rule annotate_tiberius:
    """Tiberius WHOLE-GENOME gene prediction on GPU (angiosperms model, unmasked
    input), AGAT-standardized. Tiberius restarts gene IDs (g1,g2..) per scaffold even in whole-genome mode,
    so IDs collide across scaffolds; we prefix the scaffold name onto gene_id/
    transcript_id before AGAT-standardizing to make them globally unique."""
    input:
        genome=lambda wc: frozen_chr_fasta(wc.species),
    output:
        gff="results/{species}/annotation/tiberius/{species}.tiberius.gff3",
    params:
        outdir="results/{species}/annotation/tiberius",
        sif=ANNOT_TIBERIUS_SIF,
        model=lambda wc: annotation_df.loc[wc.species, "tiberius_model"],
        agat_env=AGAT_ENV,
    resources:
        gpu=1,
        mem_mb=64000,
    threads: 8
    log:
        "logs/annotate_tiberius/{species}.log",
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        mkdir -p "{params.outdir}" "$(dirname $LOG)"
        module load singularity-binds 2>/dev/null || module load singularity 2>/dev/null || true
        {PICK_GPU}

        echo "=== Tiberius whole-genome: {wildcards.species} ===" >> "$LOG"
        singularity exec --nv \
            -B "$REPO_ROOT:$REPO_ROOT" \
            "{params.sif}" \
            tiberius.py \
            --genome "$REPO_ROOT/{input.genome}" \
            --model_cfg {params.model} \
            --out "{params.outdir}/tiberius_raw.gtf" \
            2>> "$LOG"

        # Tiberius numbers genes per-scaffold (g1,g2.. restart each sequence),
        # even in whole-genome mode. Naive use collides IDs across scaffolds
        # (e.g. g1104.t1 on 480 scaffolds) -> AGAT merges unrelated genes into
        # giant chimeras. Prefix scaffold (col1) onto gene_id/transcript_id to
        # make them globally unique before standardizing.
        echo "=== namespace IDs by scaffold ===" >> "$LOG"
        awk -F'\t' 'BEGIN{{OFS="\t"}} /^#/{{print; next}} {{
            scaf=$1
            gsub(/gene_id "/, "gene_id \"" scaf "_", $9)
            gsub(/transcript_id "/, "transcript_id \"" scaf "_", $9)
            print
        }}' "{params.outdir}/tiberius_raw.gtf" > "{params.outdir}/tiberius_raw.namespaced.gtf"

        echo "=== AGAT standardize ===" >> "$LOG"
        micromamba run -p "{params.agat_env}" agat_convert_sp_gxf2gxf.pl \
            -g "{params.outdir}/tiberius_raw.namespaced.gtf" \
            -o "$REPO_ROOT/{output.gff}" >> "$LOG" 2>&1

        echo "=== tiberius complete: {output.gff} ===" >> "$LOG"
        grep -c -P "\tgene\t" "$REPO_ROOT/{output.gff}" >> "$LOG" 2>&1 || true
        """


# ---------------------------------------------------------------------------
# BRAKER4 (RNA + protein evidence, ETP mode) — the RNA-informed annotation arm.
#
# BRAKER4 is itself a Snakemake pipeline that orchestrates ~13 Singularity
# containers. We do NOT integrate its rules into this workflow; instead this
# rule is a thin wrapper that (1) generates BRAKER4's samples.csv for the
# species, (2) runs BRAKER4's own snakemake in LOCAL mode inside one big SLURM
# allocation (no nested SLURM submission), (3) ingests the resulting GFF3.
#
# BRAKER4 runs RepeatModeler2+RepeatMasker automatically (genome_masked empty),
# aligns RNA with HISAT2 internally (ETP mode), and emits an already-clean,
# AGAT-normalized braker.gff3.gz. We standardize once more via AGAT for ID-
# prefix consistency with the other tools' outputs.
# ---------------------------------------------------------------------------

BRAKER4_REPO   = config["annotation"]["braker4_repo"]
BRAKER4_ENV    = config["annotation"]["braker4_env"]
BRAKER4_CACHE  = config["annotation"]["braker4_cache"]
ORTHODB_VIRID  = config["annotation"]["orthodb_viridiplantae"]
BRAKER4_BUSCO  = config["annotation"]["braker4_busco_lineage"]


def _rna_fastq_pairs(species: str):
    """Return (r1_colon_joined, r2_colon_joined) for a species' RNA dir.
    Globs *_R1_*.fastq.gz (sorted) and derives matching R2 by name."""
    import glob, os
    rna_dir = str(annotation_df.loc[species, "rna_dir"]).strip()
    if not rna_dir:
        raise ValueError(f"no rna_dir for {species} (BRAKER needs RNA for ETP mode)")
    r1s = sorted(glob.glob(os.path.join(rna_dir, "*_R1_*.fastq.gz")))
    if not r1s:
        raise ValueError(f"no *_R1_*.fastq.gz found in {rna_dir} for {species}")
    r2s = [r1.replace("_R1_", "_R2_") for r1 in r1s]
    for r2 in r2s:
        if not os.path.exists(r2):
            raise ValueError(f"missing R2 mate: {r2}")
    return ":".join(r1s), ":".join(r2s)
rule annotate_helixer:
    """Helixer per-chromosome gene prediction on GPU (land_plant model, unmasked
    input), merged + AGAT-standardized. Per-chr split for memory (Helixer OOMs
    on large genomes); Helixer embeds the seqname in gene IDs (_Chr01_000001),
    so per-chr outputs merge with a plain cat WITHOUT the ID collision that
    forced Tiberius to run whole-genome."""
    input:
        genome=lambda wc: frozen_chr_fasta(wc.species),
    output:
        gff="results/{species}/annotation/helixer/{species}.helixer.gff3",
    params:
        outdir="results/{species}/annotation/helixer",
        sif=ANNOT_HELIXER_SIF,
        lineage=lambda wc: annotation_df.loc[wc.species, "helixer_lineage"],
        agat_env=AGAT_ENV,
    resources:
        gpu=1,
        mem_mb=64000,
    threads: 8
    log:
        "logs/annotate_helixer/{species}.log",
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        OUT="$REPO_ROOT/{params.outdir}"
        SPLIT="$OUT/per_chr"
        ANNO="$OUT/per_chr_gff"
        mkdir -p "$SPLIT" "$ANNO" "$(dirname $LOG)"
        module load singularity-binds 2>/dev/null || module load singularity 2>/dev/null || true
        {PICK_GPU}

        SCRATCH=$(mktemp -d /scratch/$USER/helixer_{wildcards.species}_XXXXXX 2>/dev/null \
                  || mktemp -d /tmp/helixer_{wildcards.species}_XXXXXX)
        trap "rm -rf $SCRATCH" EXIT

        echo "=== split {wildcards.species} into per-chr ===" >> "$LOG"
        awk -v d="$SPLIT" '/^>/ {{ n=substr($1,2); f=d"/"n".fa" }} {{ print > f }}' \
            "$REPO_ROOT/{input.genome}"

        for chr_fa in "$SPLIT"/*.fa; do
            base=$(basename "$chr_fa" .fa)
            out_gff="$ANNO/${{base}}.helixer.gff3"
            if [ -s "$out_gff" ] && grep -q "^#" "$out_gff"; then
                echo "=== skip $base (already done) ===" >> "$LOG"; continue
            fi
            echo "=== Helixer: $base ===" >> "$LOG"
            run_tmp="$SCRATCH/${{base}}_tmp"; mkdir -p "$run_tmp"
            singularity run --nv \
                -B "$REPO_ROOT:$REPO_ROOT" \
                -B "$SCRATCH:$SCRATCH" -B "$SCRATCH:/tmp" -B "$SCRATCH:/var/tmp" \
                --env TMPDIR="$run_tmp" \
                "{params.sif}" Helixer.py \
                --fasta-path "$chr_fa" \
                --lineage {params.lineage} \
                --gff-output-path "$out_gff" \
                --temporary-dir "$run_tmp" \
                --write-by 5000000 \
                --subsequence-length 29997 \
                --batch-size 4 \
                --no-overlap \
                >> "$LOG" 2>&1
            rm -rf "$run_tmp"
        done

        echo "=== merge + AGAT standardize ===" >> "$LOG"
        cat "$ANNO"/*.helixer.gff3 > "$OUT/helixer_raw.gff3"
        micromamba run -p "{params.agat_env}" agat_convert_sp_gxf2gxf.pl \
            -g "$OUT/helixer_raw.gff3" -o "$REPO_ROOT/{output.gff}" >> "$LOG" 2>&1

        echo "=== helixer complete: {output.gff} ===" >> "$LOG"
        grep -c -P "\tgene\t" "$REPO_ROOT/{output.gff}" >> "$LOG" 2>&1 || true
        """




rule annotate_braker:
    """BRAKER4 ETP-mode annotation (RNA+protein) for one species. Runs BRAKER4's
    own snakemake locally inside this allocation; ingests + AGAT-standardizes."""
    input:
        genome=lambda wc: frozen_chr_fasta(wc.species),
        proteins=ORTHODB_VIRID,
    output:
        gff="results/{species}/annotation/braker/{species}.braker.gff3",
    params:
        workdir="results/{species}/annotation/braker",
        repo=BRAKER4_REPO,
        env=BRAKER4_ENV,
        cache=BRAKER4_CACHE,
        busco=BRAKER4_BUSCO,
        agat_env=AGAT_ENV,
        rna=lambda wc: _rna_fastq_pairs(wc.species),
    resources:
        mem_mb=400000,
        runtime=28800,      # 20 days (minutes)
    threads: 64
    log:
        "logs/annotate_braker/{species}.log",
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        WORK="$REPO_ROOT/{params.workdir}"
        mkdir -p "$WORK" "$(dirname $LOG)"

        module load singularity-binds 2>/dev/null || module load singularity 2>/dev/null || true

        R1="{params.rna[0]}"
        R2="{params.rna[1]}"

        # --- generate samples.csv (ETP: genome + proteins + RNA fastqs) ---
        # genome_masked left empty -> BRAKER4 runs RepeatModeler2+RepeatMasker
        echo "sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf" > "$WORK/samples.csv"
        echo "{wildcards.species},$REPO_ROOT/{input.genome},,{input.proteins},,$R1,$R2,,,,,,{params.busco}," >> "$WORK/samples.csv"

        # --- generate config.ini (BRAKER4 requires this file; pointed to via
        #     BRAKER4_CONFIG). Production: optimize AUGUSTUS on, OMArk+ncRNA off,
        #     gm_max_intergenic UNSET (docs: never set on real data). ---
        cat > "$WORK/config.ini" << 'INI'
[paths]
samples_file = samples.csv
augustus_config_path = augustus_config

[PARAMS]
run_best_by_compleasm = 1
run_omark = 0
skip_optimize_augustus = 0
skip_single_exon_downsampling = 0
downsampling_lambda = 2
use_varus = 0
use_compleasm_hints = 1
skip_busco = 1
use_dev_shm = 0
min_contig = 10000
no_cleanup = 0
run_ncrna = 0
INI
        export BRAKER4_CONFIG="$WORK/config.ini"

        # --- production config via env-var overrides (skip_optimize_augustus=0,
        #     run_omark=0, run_ncrna=0) ---
        export USE_SLURM=false
        export BRAKER4_SKIP_OPTIMIZE_AUGUSTUS=0
        export BRAKER4_RUN_OMARK=0
        export BRAKER4_RUN_NCRNA=0
        export BRAKER4_MEM_OF_NODE={resources.mem_mb}
        export BRAKER4_MAX_RUNTIME={resources.runtime}
        export SINGULARITYENV_PREPEND_PATH=/opt/conda/bin

        echo "=== launching BRAKER4 (ETP) for {wildcards.species} ===" >> "$LOG"
        echo "R1=$R1" >> "$LOG"; echo "R2=$R2" >> "$LOG"

        cd "$WORK"
        "{params.env}/bin/snakemake" \
            --snakefile "{params.repo}/Snakefile" \
            --cores {threads} --jobs {threads} \
            --printshellcmds --rerun-incomplete --latency-wait 120 \
            --restart-times 3 --nolock \
            --use-singularity \
            --singularity-prefix "{params.cache}" \
            --singularity-args "-B /home -B /netscratch -B /biodata --env PREPEND_PATH=/opt/conda/bin" \
            >> "$LOG" 2>&1

        # --- ingest: braker.gff3.gz -> gunzip -> AGAT standardize ---
        echo "=== ingest + AGAT standardize ===" >> "$LOG"
        RESULT="$WORK/output/{wildcards.species}/results/braker.gff3.gz"
        zcat "$RESULT" > "$WORK/braker_raw.gff3"
        micromamba run -p "{params.agat_env}" agat_convert_sp_gxf2gxf.pl \
            -g "$WORK/braker_raw.gff3" -o "$REPO_ROOT/{output.gff}" >> "$LOG" 2>&1

        echo "=== braker complete: {output.gff} ===" >> "$LOG"
        grep -c -P "\tgene\t" "$REPO_ROOT/{output.gff}" >> "$LOG" 2>&1 || true
        """


rule annotate_annevo:
    """ANNEVO v2.3.2 whole-genome ab-initio annotation on GPU (Magnoliopsida
    model), AGAT-standardized. Whole-genome is safe: ANNEVO embeds the seqname
    in gene IDs (NC_xxx-g1) so IDs are globally unique across sequences (like
    Helixer, unlike Tiberius). Uses node-scratch --tmp_path for the per-run
    model_prediction.h5 (GitHub #13: fills disk on big genomes otherwise)."""
    input:
        genome=lambda wc: frozen_chr_fasta(wc.species),
    output:
        gff="results/{species}/annotation/annevo/{species}.annevo.gff3",
    params:
        outdir="results/{species}/annotation/annevo",
        env=ANNEVO_ENV,
        repo=ANNEVO_REPO,
        model=ANNEVO_MODEL,
        lineage=lambda wc: annotation_df.loc[wc.species, "annevo_lineage"],
        agat_env=AGAT_ENV,
    resources:
        gpu=1,
        mem_mb=64000,
    threads: 16
    log:
        "logs/annotate_annevo/{species}.log",
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        OUT="$REPO_ROOT/{params.outdir}"
        mkdir -p "$OUT" "$(dirname $LOG)"
        {PICK_GPU}

        TMPBASE="/netscratch/dep_mercier/grp_marques/Aaryan/tmp/annevo"
        mkdir -p "$TMPBASE"
        SCRATCH=$(mktemp -d "$TMPBASE/{wildcards.species}_XXXXXX")
        trap "rm -rf $SCRATCH" EXIT

        echo "=== ANNEVO {wildcards.species} (GPU $CUDA_VISIBLE_DEVICES) ===" >> "$LOG"
        PYTHONNOUSERSITE=1 "{params.env}/bin/python" "{params.repo}/annotation.py" \
            -g "$REPO_ROOT/{input.genome}" \
            -m "{params.repo}/{params.model}" \
            -l {params.lineage} \
            -o "$OUT/annevo_raw.gff3" \
            --tmp_path "$SCRATCH" \
            -s 100 \
            -t {threads} \
            --batch_size 8 \
            --num_workers 1 \
            --show_log >> "$LOG" 2>&1

        echo "=== AGAT standardize ===" >> "$LOG"
        micromamba run -p "{params.agat_env}" agat_convert_sp_gxf2gxf.pl \
            -g "$OUT/annevo_raw.gff3" -o "$REPO_ROOT/{output.gff}" >> "$LOG" 2>&1

        echo "=== annevo complete: {output.gff} ===" >> "$LOG"
        grep -c -P "\tgene\t" "$REPO_ROOT/{output.gff}" >> "$LOG" 2>&1 || true
        """


rule import_external_genome:
    """Import a collaborator-assembled COLLAPSED genome into the pipeline at the
    frozen stage. Keeps the top-N chromosome-scale scaffolds (N = chr_number //
    ploidy), renames scaffold{k} -> chr{k}_collapsed (preserving the
    collaborators' scaffold numbering, since their FASTA is length-sorted and
    scaffold1..N == the N chromosomes), sorts output by chromosome number, and
    retains the dropped small-contig tail as {species}_debris.fa.

    These genomes have no HiFi/HiC accessible, so they skip
    hifiasm/standardize/decon/scaffold and feed annotation directly."""
    input:
        src=lambda wc: external_genomes_df.loc[wc.species, "source_fasta"],
    output:
        chr="results/{species}/assembly_final/external_collapsed/{species}_chr.fa",
        debris="results/{species}/assembly_final/external_collapsed/{species}_debris.fa",
    params:
        keep_n=lambda wc: keep_n_scaffolds(wc.species),
    log:
        "logs/import_external_genome/{species}.log",
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        OUT=$(dirname "$REPO_ROOT/{output.chr}")
        mkdir -p "$OUT" "$(dirname $LOG)"
        module load seqkit 2>/dev/null || true

        N={params.keep_n}
        SRC="{input.src}"
        echo "=== import_external_genome {wildcards.species}: keep scaffold1..scaffold$N ===" > "$LOG"

        # 1. kept set: scaffold1..scaffoldN -> renamed chr{{k}}_collapsed, sorted by k
        : > "$OUT/_kept_unsorted.fa"
        for k in $(seq 1 $N); do
            samtools faidx "$SRC" "scaffold$k" \
              | sed "1s/^>scaffold$k.*/>chr${{k}}_collapsed/" >> "$OUT/_kept_unsorted.fa"
        done

        # 2. sort kept records by chromosome number (chr1_collapsed..chrN_collapsed)
        #    seqkit sort by natural name order -> human-readable frozen genome
        seqkit sort --natural-order --by-name "$OUT/_kept_unsorted.fa" > "$REPO_ROOT/{output.chr}" 2>> "$LOG"
        rm -f "$OUT/_kept_unsorted.fa"

        # 3. debris = every scaffold NOT in the kept set (retain for the record)
        KEEPLIST="$OUT/_keep.list"
        for k in $(seq 1 $N); do echo "scaffold$k"; done > "$KEEPLIST"
        ALL=$(cut -f1 "$SRC.fai" 2>/dev/null || (samtools faidx "$SRC" && cut -f1 "$SRC.fai"))
        : > "$REPO_ROOT/{output.debris}"
        echo "$ALL" | grep -vxF -f "$KEEPLIST" | while read s; do
            [ -n "$s" ] && samtools faidx "$SRC" "$s" >> "$REPO_ROOT/{output.debris}"
        done
        rm -f "$KEEPLIST"

        # 4. report
        echo "kept scaffolds:" >> "$LOG"
        grep -c "^>" "$REPO_ROOT/{output.chr}" >> "$LOG"
        echo "debris scaffolds:" >> "$LOG"
        grep -c "^>" "$REPO_ROOT/{output.debris}" >> "$LOG" 2>&1 || echo 0 >> "$LOG"
        echo "=== chr names (should be chr1_collapsed..chr${{N}}_collapsed in order) ===" >> "$LOG"
        grep "^>" "$REPO_ROOT/{output.chr}" >> "$LOG"
        """


rule filter_braker:
    """Filter raw BRAKER GFF3 to evidence-supported transcripts using BRAKER's
    own gene_support.tsv. Keep a transcript if:
      - multi-exon  & >=1 intron supported by RNA/protein, OR
      - single-exon & its exon supported by RNA/protein.
    Drops the unsupported single-exon over-prediction reservoir (TE ORFs,
    fragments). 'any' = RNA or protein evidence (rescues real-but-unexpressed
    genes); a parallel rnaseq-only count is logged for sensitivity.
    Produces the production GFF3 for RNA species (Part B)."""
    input:
        gff="results/{species}/annotation/braker/{species}.braker.gff3",
        support="results/{species}/annotation/braker/output/{species}/results/gene_support.tsv",
        awk="workflow/scripts/filter_braker_by_support.awk",
    output:
        gff="results/{species}/annotation/braker/{species}.braker_filtered.gff3",
        stats="results/{species}/annotation/braker/{species}.braker_filter_stats.tsv",
    params:
        agat_env=AGAT_ENV,
        mode="any",
    log:
        "logs/filter_braker/{species}.log",
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"; mkdir -p "$(dirname $LOG)"
        TMP=$(mktemp -d)
        trap "rm -rf $TMP" EXIT

        echo "=== filter_braker {wildcards.species} (mode={params.mode}) ===" > "$LOG"

        # 1. filter to supported transcripts
        awk -f {input.awk} -v mode={params.mode} "{input.support}" "{input.gff}" > "$TMP/filtered.gff3" 2>> "$LOG"

        # 2. AGAT re-standardize
        micromamba run -p "{params.agat_env}" agat_convert_sp_gxf2gxf.pl \
            -g "$TMP/filtered.gff3" -o "$REPO_ROOT/{output.gff}" >> "$LOG" 2>&1

        # 3. stats: raw vs filtered(any) vs filtered(rnaseq)
        cntg () {{ awk -F'\t' -v f="$1" '$3==f' "$2" | wc -l; }}
        RAW_G=$(cntg gene "{input.gff}"); RAW_T=$(cntg transcript "{input.gff}")
        FA_G=$(cntg gene "$REPO_ROOT/{output.gff}"); FA_T=$(cntg transcript "$REPO_ROOT/{output.gff}")
        FR_T=$(awk -f {input.awk} -v mode=rnaseq "{input.support}" "{input.gff}" | awk -F'\t' '$3=="transcript"' | wc -l)
        FR_G=$(awk -f {input.awk} -v mode=rnaseq "{input.support}" "{input.gff}" | awk -F'\t' '$3=="gene"' | wc -l)
        {{
          echo -e "metric\traw\tfiltered_any\tfiltered_rnaseq"
          echo -e "genes\t$RAW_G\t$FA_G\t$FR_G"
          echo -e "transcripts\t$RAW_T\t$FA_T\t$FR_T"
        }} > "$REPO_ROOT/{output.stats}"
        cat "$REPO_ROOT/{output.stats}" >> "$LOG"
        echo "=== filter_braker complete ===" >> "$LOG"
        """


# tool -> standardized GFF3 path. ab-initio: {sp}.{tool}.gff3; filtered braker special-cased.
def annotation_tool_gff3(species, tool):
    if tool == "braker_filtered":
        return f"results/{species}/annotation/braker/{species}.braker_filtered.gff3"
    return f"results/{species}/annotation/{tool}/{species}.{tool}.gff3"


rule extract_proteins:
    """Extract protein FASTA from a tool's GFF3 + frozen genome (AGAT -p).
    Feeds OMArk + compleasm. tool in {annevo, helixer, tiberius, braker_filtered}."""
    input:
        gff=lambda wc: annotation_tool_gff3(wc.species, wc.tool),
        genome=lambda wc: frozen_chr_fasta(wc.species),
    output:
        faa="results/{species}/annotation/{tool}/{species}.{tool}.proteins.fa",
    params:
        agat_env=AGAT_ENV,
    wildcard_constraints:
        tool="annevo|helixer|tiberius|braker_filtered",
    log:
        "logs/extract_proteins/{species}.{tool}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        micromamba run -p "{params.agat_env}" agat_sp_extract_sequences.pl \
            -g "{input.gff}" -f "{input.genome}" -p \
            -o "{output.faa}" > "{log}" 2>&1
        echo "proteins: $(grep -c '^>' {output.faa})" >> "{log}"
        """


# ---- OMArk QC (Tier 1): per (species, tool, db) ----
OMAMER_DB = {
    "luca": config["annotation"]["omamer_db_luca"],
    "viridiplantae": config["annotation"]["omamer_db_viridiplantae"],
}
DROSERA_TAXID = 4363
OMARK_TEST_ENV = "/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/omark_test"

rule qc_omark:
    """OMArk completeness/consistency for one proteome vs one OMAmer DB.
    omamer search (multithreaded) -> omark (single-threaded, -t Drosera 4363).
    db: viridiplantae = sharp plant scoring (bake-off); luca = tree-of-life
    context for taxon attribution of unknowns (the 2-DB merge)."""
    input:
        faa="results/{species}/annotation/{tool}/{species}.{tool}.proteins.fa",
    output:
        omamer="results/{species}/annotation/{tool}/omark/{db}/{species}.{tool}.{db}.omamer",
        summ="results/{species}/annotation/{tool}/omark/{db}/{species}.{tool}.{db}.sum",
        ump="results/{species}/annotation/{tool}/omark/{db}/{species}.{tool}.{db}.ump",
        done="results/{species}/annotation/{tool}/omark/{db}/{species}.{tool}.{db}.done",
    params:
        db=lambda wc: OMAMER_DB[wc.db],
        taxid=DROSERA_TAXID,
        outdir="results/{species}/annotation/{tool}/omark/{db}",
        prefix="{species}.{tool}.{db}",
    threads: 16
    resources:
        mem_mb=32000,
        runtime=180,
    wildcard_constraints:
        tool="annevo|helixer|tiberius|braker_filtered",
        db="luca|viridiplantae",
    conda:
        "../../envs/omark.yaml"
    log:
        "logs/qc_omark/{species}.{tool}.{db}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname {log})"
        omamer search --db "{params.db}" --query "{input.faa}" \
            --out "{output.omamer}" --nthreads {threads} > "{log}" 2>&1
        omark -f "{output.omamer}" -d "{params.db}" -t {params.taxid} \
            -o "{params.outdir}" >> "{log}" 2>&1
        # omark names outputs by the omamer basename ({params.prefix}); guard it landed
        test -s "{output.summ}" || {{ echo "ERR: .sum missing"; ls -la "{params.outdir}" >> "{log}"; exit 1; }}
        test -s "{output.ump}"  || {{ echo "ERR: .ump missing"; ls -la "{params.outdir}" >> "{log}"; exit 1; }}
        touch "{output.done}"
        """


rule omark_merge:
    """Cross-reference Viridiplantae-unknown proteins against the LUCA run to
    attribute taxon (plant_orphan KEEP / nonplant DROP / unplaceable KEEP).
    Emits the drop-list consumed by refine_proteome (decontaminated ab-initio set)."""
    input:
        viridi_ump="results/{species}/annotation/{tool}/omark/viridiplantae/{species}.{tool}.viridiplantae.ump",
        luca_omamer="results/{species}/annotation/{tool}/omark/luca/{species}.{tool}.luca.omamer",
        script="workflow/scripts/omark_merge_unknowns.py",
    output:
        detail="results/{species}/annotation/{tool}/omark/{species}.{tool}.attribution.tsv",
        rollup="results/{species}/annotation/{tool}/omark/{species}.{tool}.rollup.tsv",
        droplist="results/{species}/annotation/{tool}/omark/{species}.{tool}.drop.txt",
        orphans="results/{species}/annotation/{tool}/omark/{species}.{tool}.orphans.tsv",
    params:
        env=OMARK_TEST_ENV,
    wildcard_constraints:
        tool="annevo|helixer|tiberius|braker_filtered",
    log:
        "logs/omark_merge/{species}.{tool}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        micromamba run -p "{params.env}" python3 {input.script} \
            --viridi-ump "{input.viridi_ump}" \
            --luca-omamer "{input.luca_omamer}" \
            --out-detail "{output.detail}" \
            --out-rollup "{output.rollup}" \
            --out-droplist "{output.droplist}" \
            --out-orphans "{output.orphans}" \
            --species "{wildcards.species}" --tool "{wildcards.tool}" > "{log}" 2>&1
        """


# ---- RNA-support (STAR independent splice junctions) ----
def _star_reads(species):
    """(comma-joined R1, comma-joined R2) for STAR, from the species rna_dir.
    Same fastqs BRAKER used; STAR wants comma-separated (BRAKER used colon)."""
    import glob, os
    rna_dir = str(annotation_df.loc[species, "rna_dir"]).strip()
    r1s = sorted(glob.glob(os.path.join(rna_dir, "*_R1_*.fastq.gz")))
    if not r1s:
        raise ValueError(f"no *_R1_*.fastq.gz in {rna_dir} for {species}")
    r2s = [r.replace("_R1_", "_R2_") for r in r1s]
    for r2 in r2s:
        if not os.path.exists(r2):
            raise ValueError(f"missing R2 mate: {r2}")
    return ",".join(r1s), ",".join(r2s)


rule star_index:
    """STAR genome index for RNA-support eval. NO annotation (sjdbGTFfile) ->
    junctions discovered de novo, keeping RNA-support independent of any tool's
    annotation (avoids circularity)."""
    input:
        genome=lambda wc: frozen_chr_fasta(wc.species),
    output:
        sa="results/{species}/rna_support/star_index/SAindex",
    params:
        idxdir="results/{species}/rna_support/star_index",
        sa_nbases=14,
    threads: 16
    resources:
        mem_mb=64000,
        runtime=240,
    log: "logs/star_index/{species}.log",
    shell:
        r"""
        set -euo pipefail
        module load star/v2.7.11b 2>/dev/null || true
        mkdir -p "{params.idxdir}" "$(dirname {log})"
        STAR --runMode genomeGenerate \
            --genomeDir "{params.idxdir}" \
            --genomeFastaFiles "{input.genome}" \
            --genomeSAindexNbases {params.sa_nbases} \
            --runThreadN {threads} > "{log}" 2>&1
        """


rule star_align:
    """2-pass STAR align of all RNA -> SJ.out.tab (independent splice junctions).
    The RNA-support truth set: each tool's introns scored against these later."""
    input:
        sa="results/{species}/rna_support/star_index/SAindex",
    output:
        sj="results/{species}/rna_support/star_align/SJ.out.tab",
        bam="results/{species}/rna_support/star_align/Aligned.sortedByCoord.out.bam",
    params:
        idxdir="results/{species}/rna_support/star_index",
        outprefix="results/{species}/rna_support/star_align/",
        r1=lambda wc: _star_reads(wc.species)[0],
        r2=lambda wc: _star_reads(wc.species)[1],
    threads: 16
    resources:
        mem_mb=64000,
        runtime=480,
    log: "logs/star_align/{species}.log",
    shell:
        r"""
        set -euo pipefail
        module load star/v2.7.11b 2>/dev/null || true
        mkdir -p "{params.outprefix}" "$(dirname {log})"
        STAR --runMode alignReads \
            --genomeDir "{params.idxdir}" \
            --readFilesIn "{params.r1}" "{params.r2}" \
            --readFilesCommand zcat \
            --twopassMode Basic \
            --outSAMtype BAM SortedByCoordinate \
            --runThreadN {threads} \
            --outFileNamePrefix "{params.outprefix}" > "{log}" 2>&1
        """


# ---- RNA intron-support (bake-off PRIMARY metric) ----
RNA_SUPPORT_SPECIES = ["Drosera_paradoxa", "Drosera_binata"]
RNA_SUPPORT_TOOLS = ["annevo", "helixer", "tiberius", "braker_filtered"]

rule rna_intron_support:
    """Fraction of a tool's predicted introns confirmed by STAR junctions
    (>=3 uniq reads). Independent RNA-based structural-accuracy metric for the
    ab-initio bake-off. Precision = of predicted introns, % STAR-confirmed;
    recall = of STAR junctions, % predicted. braker_filtered is partly circular
    (used RNA) -> reported as the RNA-informed reference, not an ab-initio peer."""
    input:
        gff=lambda wc: annotation_tool_gff3(wc.species, wc.tool),
        sj="results/{species}/rna_support/star_align/SJ.out.tab",
        script="workflow/scripts/rna_intron_support.py",
    output:
        tsv="results/{species}/rna_support/intron_support/{species}.{tool}.intron_support.tsv",
    params:
        min_uniq=3,
    wildcard_constraints:
        tool="annevo|helixer|tiberius|braker_filtered",
    log:
        "logs/rna_intron_support/{species}.{tool}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv})" "$(dirname {log})"
        printf 'species\ttool\tn_pred\tn_conf\tprecision\trecall\tmin_uniq\tn_truth\n' > "{output.tsv}"
        python3 {input.script} --gff "{input.gff}" --sj "{input.sj}" \
            --species "{wildcards.species}" --tool "{wildcards.tool}" \
            --min-uniq {params.min_uniq} >> "{output.tsv}" 2> "{log}"
        """

rule rna_support_aggregate:
    """Concatenate all per-(species,tool) intron-support rows into one scoreboard."""
    input:
        tsvs=expand(
            "results/{species}/rna_support/intron_support/{species}.{tool}.intron_support.tsv",
            species=RNA_SUPPORT_SPECIES, tool=RNA_SUPPORT_TOOLS,
        ),
    output:
        scoreboard="results/rna_support_scoreboard.tsv",
    shell:
        r"""
        set -euo pipefail
        first=1
        for f in {input.tsvs}; do
            if [ "$first" = 1 ]; then cat "$f" > "{output.scoreboard}"; first=0;
            else tail -n +2 "$f" >> "{output.scoreboard}"; fi
        done
        """


# ---- Final annotation: OMArk-filtered, track-aware ----
def _is_rna(species):
    return str(annotation_df.loc[species, "has_rna"]).strip().lower() in ("yes", "true", "1")

def _final_source_tool(species):
    """RNA species -> braker_filtered; no-RNA -> annevo (chosen ab-initio winner)."""
    return "braker_filtered" if _is_rna(species) else "annevo"

rule omark_classify_genes:
    """Per-gene keep/drop table from OMArk consistency + LUCA-merge attribution.
    Keep: Consistent + plant_orphan + unplaceable (+ unscored if RNA-evidenced).
    Drop: Inconsistent + Contamination + nonplant (+ unscored if ab-initio)."""
    input:
        ump="results/{species}/annotation/{tool}/omark/viridiplantae/{species}.{tool}.viridiplantae.ump",
        attribution="results/{species}/annotation/{tool}/omark/{species}.{tool}.attribution.tsv",
        gff=lambda wc: annotation_tool_gff3(wc.species, wc.tool),
        script="workflow/scripts/omark_classify_genes.py",
    output:
        tsv="results/{species}/annotation/{tool}/{species}.{tool}.gene_class.tsv",
    params:
        rna_flag=lambda wc: "--rna-evidenced" if wc.tool == "braker_filtered" else "",
    wildcard_constraints:
        tool="annevo|braker_filtered",
    log:
        "logs/omark_classify_genes/{species}.{tool}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        python3 {input.script} {params.rna_flag} \
            --ump "{input.ump}" --attribution "{input.attribution}" \
            --gff "{input.gff}" --out "{output.tsv}" \
            --species "{wildcards.species}" --tool "{wildcards.tool}" 2> "{log}"
        """

rule filter_gff3_by_class:
    """Apply the keep/drop table -> cleaned, omark_class-tagged GFF3."""
    input:
        gff=lambda wc: annotation_tool_gff3(wc.species, wc.tool),
        classes="results/{species}/annotation/{tool}/{species}.{tool}.gene_class.tsv",
        script="workflow/scripts/filter_gff3_by_class.py",
    output:
        gff="results/{species}/annotation/{tool}/{species}.{tool}.omark_clean.gff3",
    wildcard_constraints:
        tool="annevo|braker_filtered",
    log:
        "logs/filter_gff3_by_class/{species}.{tool}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        python3 {input.script} --gff "{input.gff}" \
            --classes "{input.classes}" --out "{output.gff}" 2> "{log}"
        """

rule final_annotation:
    """Authoritative annotation for downstream (OrthoFinder, FANTASIA):
    OMArk-cleaned GFF3 (track tool: RNA->braker_filtered, no-RNA->annevo) copied
    to a stable path, plus a peptide FASTA extracted from that cleaned GFF3
    (kept genes only; one protein per kept transcript)."""
    input:
        gff=lambda wc: "results/{sp}/annotation/{t}/{sp}.{t}.omark_clean.gff3".format(
            sp=wc.species, t=_final_source_tool(wc.species)),
        genome=lambda wc: frozen_chr_fasta(wc.species),
    output:
        gff="results/{species}/annotation/final/{species}.final.gff3",
        faa="results/{species}/annotation/final/{species}.final.proteins.fa",
    params:
        agat_env=AGAT_ENV,
    log:
        "logs/final_annotation/{species}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gff})" "$(dirname {log})"
        cp "{input.gff}" "{output.gff}"
        micromamba run -p "{params.agat_env}" agat_sp_extract_sequences.pl \
            -g "{output.gff}" -f "{input.genome}" -p \
            -o "{output.faa}" > "{log}" 2>&1
        echo "final proteins: $(grep -c '^>' {output.faa})" >> "{log}"
        """
