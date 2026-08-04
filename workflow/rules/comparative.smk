# Comparative genomics: genus-wide orthology across all final-annotated species.

import os as _os

# All species in annotation.csv. NOT a disk check: ORTHO_SPECIES is evaluated at
# DAG-construction time, so globbing for existing final.proteins.fa silently
# excluded any species whose annotation is produced later in the SAME run
# (scorpioides, 2026-08). Snakemake resolves the dependency itself.
ORTHO_SPECIES = list(annotation_df.index)

rule orthofinder_prep:
    """Longest-isoform-per-gene proteomes for all final-annotated species into
    one OrthoFinder input dir. Collapses isoforms (gene= field); keeps distinct
    gene ids incl. haplotype copies (ploidy-corrected downstream, not here)."""
    input:
        faas=expand("results/{sp}/annotation/final/{sp}.final.proteins.fa",
                    sp=ORTHO_SPECIES),
        script="workflow/scripts/primary_transcript.py",
    output:
        prepdir=directory("results/comparative/orthofinder/input"),
    log:
        "logs/orthofinder/prep.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{output.prepdir}" "$(dirname {log})"
        : > "{log}"
        for f in {input.faas}; do
            sp=$(basename "$f" .final.proteins.fa)
            python3 {input.script} "$f" "{output.prepdir}/${{sp}}.fa" 2>> "{log}"
        done
        echo "=== prep done: $(ls {output.prepdir}/*.fa | wc -l) species ===" >> "{log}"
        """

rule orthofinder:
    """OrthoFinder v3 all-vs-all orthology across all species (primary
    transcripts). Comprehensive (no orthogroup filtering — filter downstream).
    Fresh conda env. Thread/mem guardrailed against the all-vs-all OOM."""
    input:
        prepdir="results/comparative/orthofinder/input",
    output:
        orthogroups="results/comparative/orthofinder/out/Results_drosera/Orthogroups/Orthogroups.tsv",
    params:
        outdir="results/comparative/orthofinder/out",
        name="drosera",
    conda:
        "../envs/orthofinder.yaml"
    threads: 48
    resources:
        mem_mb=350000,
        runtime=2880,
    log:
        "logs/orthofinder/run.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        # OrthoFinder v3 refuses an -o dir that already exists; it creates it itself.
        rm -rf "{params.outdir}"
        orthofinder -f "{input.prepdir}" \
            -t {threads} -a 8 \
            -o "{params.outdir}" -n "{params.name}" \
            > "{log}" 2>&1
        echo "=== orthofinder done ===" >> "{log}"
        """
