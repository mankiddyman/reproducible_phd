# Stage 2: per-species decisions.yaml
#
# This is a Snakemake `checkpoint` (not a rule) because future stages will
# branch their DAG on the contents of decisions.yaml (e.g., assembly_strategy
# driving whether to run haphic on hap1/hap2 vs the ragtag path).
#
# For the current stage there are no downstream consumers — you can run just
# this far with `snakemake decisions --cores N` to produce and inspect the
# yaml files before building anything downstream.

checkpoint compile_decisions:
    input:
        genomescope_summary="results/{species}/qc/read/genomescope2/model/summary.txt",
    output:
        yaml="results/{species}/decisions.yaml",
    params:
        exp_ploidy=lambda wc: int(species_df.loc[wc.species, "exp_ploidy"]),
        centromere=lambda wc: str(species_df.loc[wc.species, "centromere"]),
        chr_number_2n=lambda wc: int(species_df.loc[wc.species, "chr_number_2n"]),
        notes=lambda wc: str(species_df.loc[wc.species, "notes"] or ""),
    log:
        "logs/decisions/{species}.log",
    conda:
        "../../envs/decisions.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.yaml})" logs/decisions

        python scripts/compile_decisions.py \
          --species "{wildcards.species}" \
          --exp-ploidy {params.exp_ploidy} \
          --centromere "{params.centromere}" \
          --chr-number-2n {params.chr_number_2n} \
          --notes "{params.notes}" \
          --genomescope-summary "{input.genomescope_summary}" \
          --out "{output.yaml}" \
          > {log} 2>&1
        """


# Convenience alias target: `snakemake decisions --cores N` produces
# decisions.yaml for every species.
rule decisions:
    input:
        expand("results/{species}/decisions.yaml", species=SPECIES),
