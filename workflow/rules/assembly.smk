rule run_hifiasm:
    input:
        hifi_reads=get_fastqs,
    output:
        done="results/{species}/hifiasm/{species}.done",
        collapsed_gfa="results/{species}/hifiasm/{species}.p_ctg.gfa",
        hap1_gfa="results/{species}/hifiasm/{species}.hap1.p_ctg.gfa",
        hap2_gfa="results/{species}/hifiasm/{species}.hap2.p_ctg.gfa",
        rutg_gfa="results/{species}/hifiasm/{species}.r_utg.gfa",
    params:
        out_prefix="results/{species}/hifiasm/{species}",
        hic_flags=lambda wc: hifiasm_hic_flags(wc.species),
        extra=config["hifiasm"].get("extra", ""),
    threads: 32
    resources:
        mem_mb=lambda wc, attempt: 400000 * attempt,
        runtime=lambda wc, attempt: 1440 * attempt,
    retries: 3
    log:
        "logs/hifiasm/{species}.log"
    benchmark:
        "benchmarks/hifiasm/{species}.tsv"
    conda:
        "../../envs/hifiasm.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/{wildcards.species}/hifiasm logs/hifiasm benchmarks/hifiasm

        hifiasm \
            -o {params.out_prefix} \
            -t {threads} \
            {params.extra} \
            {params.hic_flags} \
            {input.hifi_reads} \
            > {log} 2>&1

        # Hifiasm uses different output prefixes depending on mode (hic vs bp).
        # Symlink whichever set exists to canonical names so downstream rules
        # don't need to know which mode produced them.
        prefix="{params.out_prefix}"
        if [ -f "${{prefix}}.hic.p_ctg.gfa" ]; then
            mode_infix="hic"
        elif [ -f "${{prefix}}.bp.p_ctg.gfa" ]; then
            mode_infix="bp"
        else
            echo "ERROR: hifiasm produced no recognised output" >&2
            ls -la "${{prefix}}".*.gfa 2>&1 || true
            exit 1
        fi
        echo "Detected hifiasm output mode: ${{mode_infix}}" >> {log}

        ln -srf "${{prefix}}.${{mode_infix}}.p_ctg.gfa"       "${{prefix}}.p_ctg.gfa"
        ln -srf "${{prefix}}.${{mode_infix}}.hap1.p_ctg.gfa"  "${{prefix}}.hap1.p_ctg.gfa"
        ln -srf "${{prefix}}.${{mode_infix}}.hap2.p_ctg.gfa"  "${{prefix}}.hap2.p_ctg.gfa"
        ln -srf "${{prefix}}.${{mode_infix}}.r_utg.gfa"       "${{prefix}}.r_utg.gfa"

        touch {output.done}
        """


rule standardize_initial_assembly:
    input:
        done="results/{species}/hifiasm/{species}.done",
        collapsed_gfa_src="results/{species}/hifiasm/{species}.p_ctg.gfa",
        hap1_gfa_src="results/{species}/hifiasm/{species}.hap1.p_ctg.gfa",
        hap2_gfa_src="results/{species}/hifiasm/{species}.hap2.p_ctg.gfa",
        rutg_gfa_src="results/{species}/hifiasm/{species}.r_utg.gfa",
    output:
        collapsed_gfa="results/{species}/assembly/initial/collapsed/{species}.gfa",
        hap1_gfa="results/{species}/assembly/initial/hap1/{species}.gfa",
        hap2_gfa="results/{species}/assembly/initial/hap2/{species}.gfa",
        rutg_gfa="results/{species}/assembly/initial/r_utg/{species}.gfa",
        manifest="results/{species}/assembly/initial/manifest.tsv",
    log:
        "logs/assembly_standardization/{species}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
          results/{wildcards.species}/assembly/initial/collapsed \
          results/{wildcards.species}/assembly/initial/hap1 \
          results/{wildcards.species}/assembly/initial/hap2 \
          results/{wildcards.species}/assembly/initial/r_utg \
          logs/assembly_standardization

        ln -sf "$(realpath {input.collapsed_gfa_src})" {output.collapsed_gfa}
        ln -sf "$(realpath {input.hap1_gfa_src})" {output.hap1_gfa}
        ln -sf "$(realpath {input.hap2_gfa_src})" {output.hap2_gfa}
        ln -sf "$(realpath {input.rutg_gfa_src})" {output.rutg_gfa}

        cat > {output.manifest} << EOF
role	gfa
collapsed	{output.collapsed_gfa}
hap1	{output.hap1_gfa}
hap2	{output.hap2_gfa}
r_utg	{output.rutg_gfa}
EOF

        echo "Standardized initial GFA outputs for {wildcards.species}" > {log}
        """


rule gfa_to_fasta_initial:
    input:
        collapsed_gfa="results/{species}/assembly/initial/collapsed/{species}.gfa",
        hap1_gfa="results/{species}/assembly/initial/hap1/{species}.gfa",
        hap2_gfa="results/{species}/assembly/initial/hap2/{species}.gfa",
        rutg_gfa="results/{species}/assembly/initial/r_utg/{species}.gfa",
    output:
        collapsed_fa="results/{species}/assembly/initial/collapsed/{species}.fa",
        hap1_fa="results/{species}/assembly/initial/hap1/{species}.fa",
        hap2_fa="results/{species}/assembly/initial/hap2/{species}.fa",
        rutg_fa="results/{species}/assembly/initial/r_utg/{species}.fa",
    log:
        "logs/gfa_to_fasta/{species}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p logs/gfa_to_fasta

        awk '/^S/{{print ">"$2"\n"$3}}' {input.collapsed_gfa} > {output.collapsed_fa}
        awk '/^S/{{print ">"$2"\n"$3}}' {input.hap1_gfa} > {output.hap1_fa}
        awk '/^S/{{print ">"$2"\n"$3}}' {input.hap2_gfa} > {output.hap2_fa}
        awk '/^S/{{print ">"$2"\n"$3}}' {input.rutg_gfa} > {output.rutg_fa}

        echo "Converted initial GFA outputs to FASTA for {wildcards.species}" > {log}
        """
