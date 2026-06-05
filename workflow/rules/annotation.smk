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

ANNOT_TIBERIUS_SIF = config["annotation"]["tiberius_sif"]
ANNOT_HELIXER_SIF  = config["annotation"]["helixer_sif"]
GFFREAD            = config["annotation"]["gffread"]


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
    """Tiberius per chromosome on GPU (angiosperms model, unmasked input),
    merged + gffread-standardized."""
    input:
        split_done="results/{species}/annotation/_split/.split.done",
    output:
        gff="results/{species}/annotation/tiberius/{species}.tiberius.gff3",
    params:
        split_dir="results/{species}/annotation/_split",
        outdir="results/{species}/annotation/tiberius",
        sif=ANNOT_TIBERIUS_SIF,
        model=lambda wc: annotation_df.loc[wc.species, "tiberius_model"],
        gffread=GFFREAD,
    resources:
        gpu=1,
        mem_mb=32000,
    threads: 8
    log:
        "logs/annotate_tiberius/{species}.log",
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        mkdir -p "{params.outdir}/per_chr" "$(dirname $LOG)"

        module load singularity-binds 2>/dev/null || module load singularity 2>/dev/null || true
        {PICK_GPU}

        for chr_fa in {params.split_dir}/*.fa; do
            base=$(basename "$chr_fa" .fa)
            echo "=== Tiberius: $base ===" >> "$LOG"
            singularity exec --nv \
                -B "$REPO_ROOT:$REPO_ROOT" \
                "{params.sif}" \
                tiberius.py \
                --genome "$REPO_ROOT/$chr_fa" \
                --model_cfg {params.model} \
                --out "{params.outdir}/per_chr/${{base}}.tiberius.gff3" \
                2>> "$LOG"
        done

        echo "=== merge + gffread standardize ===" >> "$LOG"
        cat {params.outdir}/per_chr/*.tiberius.gff3 > "{params.outdir}/merged_raw.gff3"
        {params.gffread} "{params.outdir}/merged_raw.gff3" \
            -o "{output.gff}" --keep-genes 2>> "$LOG"

        echo "=== tiberius complete: {output.gff} ===" >> "$LOG"
        grep -c -P "\tgene\t" "{output.gff}" >> "$LOG" 2>&1 || true
        """
