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

# OrthoFinder v3 (all-vs-all orthology on primary transcripts, comprehensive: no
# orthogroup filtering, filter downstream) runs as two jobs since the 19-species
# run (2026-11: 11 Droseraceae/Nepenthes + holocentric/monocentric pairs from 4
# other clades). All four 19-species attempts (hpc16, 2026-11-03 and -04) died
# with "[Errno 12] Cannot allocate memory" at "Starting STRIDE", the rooting of
# the species tree: STRIDE runs a pool of -t worker processes, each a fork of
# OrthoFinder's main process (~28 GB with 19 species; a fresh fork for every 10
# gene trees), and under the cluster's strict memory overcommit 120 such forks
# cannot be committed. Every other step forks at most -a copies of the main
# process, and runs DIAMOND, FAMSA and FastTree through a helper process forked
# at startup while OrthoFinder is still small. Hence:
#   orthofinder_trees  -t 120: search, orthogroups, alignments, gene trees and
#                      the unrooted species tree (-ogt stops there);
#   orthofinder        -t 4, from those trees (-fgt): STRIDE rooting,
#                      orthologues, HOGs, and the results (Orthogroups.tsv
#                      included) in out/Results_drosera/ as before.
# Tested with OrthoFinder 3.1.5 (bioconda): the two jobs give the same
# orthogroups, gene counts, HOGs, orthologues and rooted species tree as one
# uninterrupted run. OrthoFinder exits 0 even when it fails (its finally: block
# calls sys.exit()), so each rule checks that OrthoFinder got to its end.
# Slurm values are in profiles/slurm/config.yaml (set-resources); the node's
# memory and commit headroom at the start go to logs/orthofinder/node_memory*.txt.

rule orthofinder_trees:
    """OrthoFinder up to the gene trees and the unrooted species tree."""
    input:
        prepdir="results/comparative/orthofinder/input",
    output:
        trees=directory("results/comparative/orthofinder/out/Results_drosera_trees"),
    params:
        outdir="results/comparative/orthofinder/out",
        name="drosera_trees",
    conda:
        "../envs/orthofinder.yaml"
    threads: 120
    resources:
        mem_mb=900000,
        runtime=5760,
    log:
        "logs/orthofinder/trees.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        (hostname; date; grep -E '^(MemTotal|MemAvailable|CommitLimit|Committed_AS)' /proc/meminfo) \
            > "$(dirname {log})/node_memory_trees.txt"
        # OrthoFinder v3 refuses an -o dir that already exists; it creates it itself.
        rm -rf "{params.outdir}"
        orthofinder -f "{input.prepdir}" \
            -t {threads} -a 4 -ogt \
            -o "{params.outdir}" -n "{params.name}" \
            > "{log}" 2>&1
        grep -q "stopped after gene tree" "{log}" \
            && grep -q "^WorkingDirectory_Trees: " "{output.trees}/Log.txt" \
            && [ -s "{output.trees}/WorkingDirectory/SpeciesTree_unrooted_ids.txt" ] \
            || {{ echo "OrthoFinder did not get to the end of the gene trees: see {log}" >&2; exit 1; }}
        echo "=== orthofinder trees done ===" >> "{log}"
        """

rule orthofinder:
    """OrthoFinder from the gene trees on: rooted species tree, orthologues,
    HOGs and all result files, in out/Results_drosera/."""
    input:
        trees="results/comparative/orthofinder/out/Results_drosera_trees",
    output:
        orthogroups="results/comparative/orthofinder/out/Results_drosera/Orthogroups/Orthogroups.tsv",
    params:
        resdir="results/comparative/orthofinder/out/Results_drosera",
        name="drosera",
    conda:
        "../envs/orthofinder.yaml"
    threads: 4
    resources:
        mem_mb=900000,
        runtime=2880,
    log:
        "logs/orthofinder/run.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        (hostname; date; grep -E '^(MemTotal|MemAvailable|CommitLimit|Committed_AS)' /proc/meminfo) \
            > "$(dirname {log})/node_memory.txt"
        # An earlier attempt's results must go: OrthoFinder would write Results_drosera_1 next to them.
        rm -rf "{params.resdir}"
        orthofinder -fgt "{input.trees}" \
            -t {threads} -a {threads} \
            -n "{params.name}" \
            > "{log}" 2>&1
        grep -q "OrthoFinder run completed" "{params.resdir}/Log.txt" \
            || {{ echo "OrthoFinder did not complete: see {log}" >&2; exit 1; }}
        echo "=== orthofinder done ===" >> "{log}"
        """
