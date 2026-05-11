# Stage 3: per-species blobtoolkit on the initial collapsed assembly.
#
# Wraps sanger-tol/blobtoolkit (Nextflow) as a Snakemake rule.
# Parameters are passed via a YAML params file (Nextflow CLI cannot pass
# booleans cleanly when nf-schema strict validation is enabled — see
# https://github.com/nextflow-io/nextflow/issues/6760).

rule make_blobtoolkit_samplesheet:
    input:
        species_table=config["species_table"],
    output:
        samplesheet="results/{species}/blobtoolkit/initial/samplesheet.csv",
    log:
        "logs/blobtoolkit/{species}_samplesheet.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.samplesheet})" logs/blobtoolkit
        python3 scripts/make_blobtoolkit_samplesheet.py \
            --species "{wildcards.species}" \
            --species-table "{input.species_table}" \
            --out "{output.samplesheet}" \
            > {log} 2>&1
        """


rule make_blobtoolkit_params:
    input:
        fasta="results/{species}/assembly/initial/collapsed/{species}.fa",
        samplesheet="results/{species}/blobtoolkit/initial/samplesheet.csv",
    output:
        params="results/{species}/blobtoolkit/initial/params.yaml",
    params:
        outdir=lambda wc: f"results/{wc.species}/blobtoolkit/initial/output",
        taxdump=config["blobtoolkit"]["taxdump"],
        blastn=config["blobtoolkit"]["blastn_nal"],
        blastp=config["blobtoolkit"]["blastp_dmnd"],
        blastx=config["blobtoolkit"]["blastx_dmnd"],
        taxon=4363,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.params})"

        cat > {output.params} <<YAML
input: $(realpath {input.samplesheet})
outdir: $(realpath -m {params.outdir})
fasta: $(realpath {input.fasta})
align: true
taxon: {params.taxon}
taxdump: {params.taxdump}
blastn: {params.blastn}
blastp: {params.blastp}
blastx: {params.blastx}
skip_taxon_filtering: true
use_work_dir_as_temp: true
YAML
        """


rule run_blobtoolkit_initial:
    input:
        fasta="results/{species}/assembly/initial/collapsed/{species}.fa",
        samplesheet="results/{species}/blobtoolkit/initial/samplesheet.csv",
        params_file="results/{species}/blobtoolkit/initial/params.yaml",
    output:
        sentinel=touch("results/{species}/blobtoolkit/initial/.done"),
    params:
        outdir="results/{species}/blobtoolkit/initial/output",
        workdir=lambda wc: f"{config['blobtoolkit']['nextflow_work_base']}/{wc.species}",
        revision=config["blobtoolkit"]["nextflow_revision"],
        nf_config="workflow/nextflow.config",
    log:
        "logs/blobtoolkit/{species}_run.log",
    threads: 1
    resources:
        nextflow_slot=1,
    shell:
        r"""
        set -euo pipefail
        export NXF_OPTS='-Xmx8g -Xms2g'

        mkdir -p "{params.outdir}" "{params.workdir}" logs/blobtoolkit

        LOG_ABS="$(realpath {log})"

        nextflow run sanger-tol/blobtoolkit \
            -r {params.revision} \
            -profile singularity \
            -c {params.nf_config} \
            -work-dir "{params.workdir}" \
            -resume \
            -params-file "{input.params_file}" \
            > "$LOG_ABS" 2>&1
        """


rule blobtoolkit_all:
    input:
        expand("results/{species}/blobtoolkit/initial/.done", species=SPECIES),
