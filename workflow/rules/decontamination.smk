# Stage 4: decontamination of initial assemblies based on blobtoolkit output.
#
# Per (species, assembly), filters contigs using the buscoregions_class
# taxonomy assignments from blobtoolkit. Default keep set is defined in
# config/decontamination.yaml. Per-species overrides supported there.
#
# Outputs:
#   results/{species}/assembly/decontaminated/{assembly}/
#     ├── {species}.fa                          decontaminated FASTA
#     ├── keep_contigs.tsv / drop_contigs.tsv   per-contig decisions
#     ├── keep_contigs.list                     for seqkit grep
#     ├── decontamination_audit.yaml            full provenance
#     └── stats.tsv                             summary
#
# Per-species assembly targets (hap1+hap2 for diploids, p_utg for polyploids)
# come from decontamination_targets(species) in common.smk.

wildcard_constraints:
    assembly=r"hap\d+|p_utg",

rule decontaminate_initial:
    input:
        fasta="results/{species}/assembly/initial/{assembly}/{species}.fa",
        blobdir_done="results/{species}/blobtoolkit/initial/{assembly}/.done",
        config="config/decontamination.yaml",
    output:
        fasta="results/{species}/assembly/decontaminated/{assembly}/{species}.fa",
        keep_tsv="results/{species}/assembly/decontaminated/{assembly}/keep_contigs.tsv",
        drop_tsv="results/{species}/assembly/decontaminated/{assembly}/drop_contigs.tsv",
        keep_list="results/{species}/assembly/decontaminated/{assembly}/keep_contigs.list",
        audit="results/{species}/assembly/decontaminated/{assembly}/decontamination_audit.yaml",
        stats="results/{species}/assembly/decontaminated/{assembly}/stats.tsv",
    params:
        blobdir=lambda wc: f"results/{wc.species}/blobtoolkit/initial/{wc.assembly}/output/blobtoolkit/{wc.species}",
    log:
        "logs/decontamination/{species}_{assembly}.log",
    threads: 4
    resources:
        mem_mb=8000,
    conda:
        "../../envs/seqkit.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fasta})" logs/decontamination

        python3 workflow/scripts/decontaminate_assembly.py \
            --blobdir "{params.blobdir}" \
            --config "{input.config}" \
            --species "{wildcards.species}" \
            --assembly "{wildcards.assembly}" \
            --out-keep-tsv "{output.keep_tsv}" \
            --out-drop-tsv "{output.drop_tsv}" \
            --out-keep-list "{output.keep_list}" \
            --out-audit-yaml "{output.audit}" \
            --out-stats-tsv "{output.stats}" \
            >> {log} 2>&1

        # Extract kept contigs from the input FASTA
        seqkit grep -n -f "{output.keep_list}" "{input.fasta}" -o "{output.fasta}" \
            2>> {log}

        echo "Decontamination complete for {wildcards.species}/{wildcards.assembly}" >> {log}
        """


rule decontaminate_gfa:
    """Filter the initial hap GFA to the contigs kept after decontamination,
    producing a GFA matched to the decontaminated FASTA. Used as HapHiC --gfa
    input for the onepass_gfa scaffolding strategy (the assembly graph carries
    haplotype phasing).

    Keeps: S-lines whose segment ∈ keep_contigs.list; L-lines whose BOTH
    endpoints survive; A-lines whose contig survives; header (H) lines.
    Drops links/records referencing removed (contaminant) contigs.
    """
    input:
        gfa="results/{species}/assembly/initial/{assembly}/{species}.gfa",
        keep="results/{species}/assembly/decontaminated/{assembly}/keep_contigs.list",
    output:
        gfa="results/{species}/assembly/decontaminated/{assembly}/{species}.gfa",
    wildcard_constraints:
        assembly=r"hap\d+",
    log:
        "logs/decontaminate_gfa/{species}_{assembly}.log",
    threads: 1
    resources:
        mem_mb=8000,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gfa})" logs/decontaminate_gfa

        awk '
          NR==FNR {{ keep[$1]=1; next }}
          /^S/    {{ if ($2 in keep) print; next }}
          /^L/    {{ if (($2 in keep) && ($4 in keep)) print; next }}
          /^A/    {{ if ($2 in keep) print; next }}
          /^H/    {{ print; next }}
                  {{ print }}
        ' {input.keep} {input.gfa} > {output.gfa} 2> {log}

        echo "decontaminate_gfa {wildcards.species}/{wildcards.assembly}: kept $(grep -c '^S' {output.gfa}) S, $(grep -c '^L' {output.gfa}) L, $(grep -c '^A' {output.gfa}) A lines" >> {log}
        """


rule aggregate_decontamination_report:
    input:
        # Depend on all per-species audit yamls (uses common.smk helper)
        audits=lambda wc: [
            f"results/{sp}/assembly/decontaminated/{tgt}/decontamination_audit.yaml"
            for sp in SPECIES
            for tgt in decontamination_targets(sp)
        ],
        stats=lambda wc: [
            f"results/{sp}/assembly/decontaminated/{tgt}/stats.tsv"
            for sp in SPECIES
            for tgt in decontamination_targets(sp)
        ],
        # Also include the per-assembly drop_contigs for taxonomy aggregation
        drop_tsvs=lambda wc: [
            f"results/{sp}/assembly/decontaminated/{tgt}/drop_contigs.tsv"
            for sp in SPECIES
            for tgt in decontamination_targets(sp)
        ],
        # Compleasm before/after for the QC comparison
        compleasm_initial=lambda wc: [
            f"results/{sp}/qc/compleasm/initial/{tgt}/summary.txt"
            for sp in SPECIES
            for tgt in decontamination_targets(sp)
            if (Path(f"results/{sp}/qc/compleasm/initial/{tgt}/summary.txt")).is_file() or True
        ],
        compleasm_decontaminated=lambda wc: [
            f"results/{sp}/qc/compleasm/decontaminated/{tgt}/summary.txt"
            for sp in SPECIES
            for tgt in decontamination_targets(sp)
        ],
    output:
        html="results/aggregated/decontamination_report.html",
        tsv="results/aggregated/decontamination_summary.tsv",
    log:
        "logs/decontamination/aggregate_report.log",
    threads: 1
    resources:
        mem_mb=4000,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.html})" logs/decontamination

        python3 workflow/scripts/aggregate_decontamination_report.py \
            --html "{output.html}" \
            --tsv "{output.tsv}" \
            > {log} 2>&1
        """
