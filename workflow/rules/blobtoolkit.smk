# Stage 3: per-species, per-haplotype blobtoolkit on the initial assembly.
#
# Wraps sanger-tol/blobtoolkit (Nextflow) as a Snakemake rule, parameterized
# by haplotype ({hap} ∈ {hap1, hap2}). Each (species, hap) gets isolated
# nextflow work and launch directories so all combinations can run in parallel
# without session lock collisions.
#
# Parameters are passed via a YAML params file because Nextflow CLI cannot pass
# booleans cleanly when nf-schema strict validation is enabled — see
# https://github.com/nextflow-io/nextflow/issues/6760.

# Constrain {hap} so it cannot match anything other than hap1 or hap2.
# Without this, snakemake could try to match e.g. "collapsed" and explode the DAG.
wildcard_constraints:
    hap=r"hap[12]",


HAPS = ["hap1", "hap2"]


rule make_blobtoolkit_samplesheet:
    input:
        species_table=config["species_table"],
    output:
        samplesheet="results/{species}/blobtoolkit/initial/{hap}/samplesheet.csv",
    log:
        "logs/blobtoolkit/{species}_{hap}_samplesheet.log",
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
        fasta="results/{species}/assembly/initial/{hap}/{species}.fa",
        samplesheet="results/{species}/blobtoolkit/initial/{hap}/samplesheet.csv",
    output:
        params="results/{species}/blobtoolkit/initial/{hap}/params.yaml",
    params:
        outdir=lambda wc: f"results/{wc.species}/blobtoolkit/initial/{wc.hap}/output",
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
        fasta="results/{species}/assembly/initial/{hap}/{species}.fa",
        samplesheet="results/{species}/blobtoolkit/initial/{hap}/samplesheet.csv",
        params_file="results/{species}/blobtoolkit/initial/{hap}/params.yaml",
    output:
        sentinel=touch("results/{species}/blobtoolkit/initial/{hap}/.done"),
    params:
        outdir=lambda wc: f"results/{wc.species}/blobtoolkit/initial/{wc.hap}/output",
        workdir=lambda wc: f"{config['blobtoolkit']['nextflow_work_base']}/{wc.species}_{wc.hap}",
        launchdir=lambda wc: f"{config['blobtoolkit']['nextflow_work_base']}/{wc.species}_{wc.hap}_launch",
        revision=config["blobtoolkit"]["nextflow_revision"],
        nf_config="workflow/nextflow.config",
    log:
        "logs/blobtoolkit/{species}_{hap}_run.log",
    threads: 1
    resources:
        nextflow_slot=1,
    shell:
        r"""
        set -euo pipefail
        export NXF_OPTS='-Xmx8g -Xms2g'

        OUTDIR_ABS="$(realpath -m {params.outdir})"
        WORKDIR_ABS="{params.workdir}"
        LAUNCH_ABS="{params.launchdir}"
        PARAMS_ABS="$(realpath {input.params_file})"
        NF_CONFIG_ABS="$(realpath {params.nf_config})"
        LOG_ABS="$(realpath -m {log})"

        mkdir -p "$OUTDIR_ABS" "$WORKDIR_ABS" "$LAUNCH_ABS" "$(dirname $LOG_ABS)"

        # Per-(species, hap) launch dir gives each nextflow invocation its own
        # .nextflow/ state, history, and log — prevents session lock collisions
        # when 14 combinations run in parallel.
        cd "$LAUNCH_ABS"

        nextflow run sanger-tol/blobtoolkit \
            -r {params.revision} \
            -profile singularity \
            -c "$NF_CONFIG_ABS" \
            -work-dir "$WORKDIR_ABS" \
            -resume \
            -params-file "$PARAMS_ABS" \
            > "$LOG_ABS" 2>&1
        """


rule blobtoolkit_all:
    input:
        expand(
            "results/{species}/blobtoolkit/initial/{hap}/.done",
            species=SPECIES, hap=HAPS,
        ),
