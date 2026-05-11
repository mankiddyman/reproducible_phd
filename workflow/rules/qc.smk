rule assembly_stats_initial:
    input:
        fa="results/{species}/assembly/initial/{role}/{species}.fa"
    output:
        stats="results/{species}/qc/assembly_stats/initial/{role}.tsv"
    log:
        "logs/assembly_stats/{species}_{role}.log"
    conda:
        "../../envs/assembly_stats.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p results/{wildcards.species}/qc/assembly_stats/initial logs/assembly_stats

        assembly-stats {input.fa} > {output.stats} 2> {log}
        """

rule compleasm_initial:
    input:
        fa="results/{species}/assembly/initial/{role}/{species}.fa"
    output:
        summary="results/{species}/qc/compleasm/initial/{role}/summary.txt"
    log:
        "logs/compleasm/{species}_{role}.log"
    threads: 16
    conda:
        "../../envs/compleasm.yaml"
    params:
        outdir="results/{species}/qc/compleasm/initial/{role}",
        lineage=config["compleasm"]["lineage"],
        library_path=config["compleasm"].get("library_path", "resources/compleasm_lineages")
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir} logs/compleasm

        compleasm run \
            -a {input.fa} \
            -o {params.outdir} \
            -l {params.lineage} \
            -L {params.library_path} \
            -t {threads} \
            > {log} 2>&1
        """



rule jellyfish_count_genomescope:
    input:
        reads=get_fastqs
    output:
        jf="results/{species}/qc/read/genomescope2/jellyfish/{species}.jf"
    params:
        k=lambda wc: int(species_df.loc[wc.species, "k"]),
        hash_size=config["genomescope2"].get("hash_size", 20000000000)
    threads:
        config["genomescope2"].get("threads", 16)
    log:
        "logs/genomescope2/jellyfish_count/{species}.log"
    conda:
        "../../envs/genomescope2.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
          results/{wildcards.species}/qc/read/genomescope2/jellyfish \
          logs/genomescope2/jellyfish_count

        jellyfish count \
          -C \
          -m {params.k} \
          -s {params.hash_size} \
          -t {threads} \
          -o {output.jf} \
          <(zcat {input.reads}) \
          > {log} 2>&1
        """


rule jellyfish_histo_genomescope:
    input:
        jf="results/{species}/qc/read/genomescope2/jellyfish/{species}.jf"
    output:
        histo="results/{species}/qc/read/genomescope2/jellyfish/{species}.histo"
    threads:
        config["genomescope2"].get("threads", 16)
    log:
        "logs/genomescope2/jellyfish_histo/{species}.log"
    conda:
        "../../envs/genomescope2.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p logs/genomescope2/jellyfish_histo

        jellyfish histo \
          -t {threads} \
          {input.jf} \
          > {output.histo} 2> {log}
        """


rule genomescope2_initial:
    input:
        histo="results/{species}/qc/read/genomescope2/jellyfish/{species}.histo"
    output:
        summary="results/{species}/qc/read/genomescope2/model/summary.txt"
    params:
        k=lambda wc: int(species_df.loc[wc.species, "k"]),
        outdir="results/{species}/qc/read/genomescope2/model"
    log:
        "logs/genomescope2/model/{species}.log"
    conda:
        "../../envs/genomescope2.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} logs/genomescope2/model
        genomescope2 -i {input.histo} -o {params.outdir} -k {params.k} > {log} 2>&1
        """
