rule assembly_stats_initial:
    input:
        fa="results/{species}/assembly/initial/collapsed/{species}.fa"
    output:
        stats="results/{species}/qc/assembly_stats/{species}.tsv"
    log:
        "logs/assembly_stats/{species}.log"
    conda:
        "../../envs/assembly_stats.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p results/{wildcards.species}/qc/assembly_stats logs/assembly_stats

        assembly-stats {input.fa} > {output.stats} 2> {log}
        """

rule compleasm_initial:
    input:
        fa="results/{species}/assembly/initial/collapsed/{species}.fa"
    output:
        done="results/{species}/qc/compleasm/{species}/.done"
    log:
        "logs/compleasm/{species}.log"
    conda:
        "../../envs/compleasm.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p results/{wildcards.species}/qc/compleasm/{wildcards.species} logs/compleasm

        # fill in real compleasm command after checking --help
        compleasm.py --help > {log} 2>&1

        touch {output.done}
        """
