# Finalization: freeze curated scaffolds and apply chrN_hapM naming.
#
#   freeze_scaffolds  : latest Juicebox review.assembly -> juicer post ->
#                       {species}.frozen.fa  (pre-rename, tracked)
#   rename_scaffolds  : apply config/naming/{species}.csv ->
#                       {species}_chr.fa (named chromosomes) +
#                       {species}_unplaced.fa (everything else)
#
# stage wildcard = pass1 | onepass | pass2 (which scaffolding path was curated).

rule freeze_scaffolds:
    """Generate the final curated FASTA from the latest Juicebox review.assembly
    via juicer post. Frozen pre-rename assembly; rename_scaffolds names it after."""
    input:
        review="results/{species}/scaffolding/{stage}/04.build/out_JBAT.review.assembly",
        liftover="results/{species}/scaffolding/{stage}/04.build/out_JBAT.liftover.agp",
        ref=lambda wc: _glob.glob(f"results/{wc.species}/scaffolding/{wc.stage}/work/*.fa")[0],
    output:
        frozen_fa="results/{species}/assembly_final/{stage}/{species}.frozen.fa",
    params:
        prefix="results/{species}/assembly_final/{stage}/{species}.frozen",
        juicer=HAPHIC_JUICER_POST,
    wildcard_constraints:
        stage=r"pass1|onepass|pass2",
    threads: 2
    resources:
        mem_mb=12000,
    conda:
        HAPHIC_ENV
    log:
        "logs/freeze_scaffolds/{species}_{stage}.log"
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        mkdir -p "$(dirname {output.frozen_fa})" "$(dirname {log})"
        export PATH=/usr/bin:$PATH   # system java for juicer
        "$REPO_ROOT/{params.juicer}" post -o "{params.prefix}" \
            "{input.review}" "{input.liftover}" "{input.ref}" 2>> {log}
        mv "{params.prefix}.FINAL.fa" "{output.frozen_fa}"
        echo "Frozen: {output.frozen_fa}" >> {log}
        """



rule rename_scaffolds:
    """Apply config/naming/{species}.csv: mapped scaffolds -> chrN_hapM in
    {species}_chr.fa (map order); all others -> unplaced_NNN (by length) in
    {species}_unplaced.fa. The naming map is hand-authored per species since
    chromosome/haplotype homology assignment is not recoverable from order."""
    input:
        fa="results/{species}/assembly_final/{stage}/{species}.frozen.fa",
        namemap="config/naming/{species}.csv",
    output:
        chr_fa="results/{species}/assembly_final/{stage}/{species}_chr.fa",
        unplaced_fa="results/{species}/assembly_final/{stage}/{species}_unplaced.fa",
    wildcard_constraints:
        stage=r"pass1|onepass|pass2",
    threads: 1
    resources:
        mem_mb=8000,
    log:
        "logs/rename_scaffolds/{species}_{stage}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.chr_fa})" "$(dirname {log})"
        python3 workflow/scripts/rename_scaffolds.py \
            --fasta {input.fa} --map {input.namemap} \
            --out-chr {output.chr_fa} --out-unplaced {output.unplaced_fa} \
            > {log} 2>&1
        cat {log}
        """


# ---- Freeze decontaminated (no-HiC) assemblies ----
# For species without HiC scaffolding: promote the decontaminated hifiasm output
# to assembly_final/{stage}/{species}_chr.fa, the same stable path the HiC path
# produces via rename_scaffolds. Two modes (frozen_stage in annotation.csv):
#   dualhap   : concat decontaminated/hap1 + hap2 (hifiasm h1*/h2* names are
#               already disjoint, so cat is collision-free) -> diploids
#   collapsed : copy decontaminated/p_utg -> polyploids (primary unitigs)
def _freeze_decontam_inputs(wc):
    sp = wc.species
    if wc.stage == "dualhap":
        return [f"results/{sp}/assembly/decontaminated/hap1/{sp}.fa",
                f"results/{sp}/assembly/decontaminated/hap2/{sp}.fa"]
    elif wc.stage == "collapsed":
        return [f"results/{sp}/assembly/decontaminated/p_utg/{sp}.fa"]
    raise ValueError(f"freeze_decontaminated: bad stage {wc.stage}")

rule freeze_decontaminated:
    """Promote decontaminated (no-HiC) assembly to the stable frozen path.
    dualhap -> concat hap1+hap2 (collision-free h1*/h2* contig names);
    collapsed -> copy p_utg. Output is contig-level (no chromosome naming;
    these species lack HiC), which is fine for annotation."""
    input:
        fas=_freeze_decontam_inputs,
    output:
        chr_fa="results/{species}/assembly_final/{stage}/{species}_chr.fa",
    wildcard_constraints:
        stage=r"dualhap|collapsed",
    log:
        "logs/freeze_decontaminated/{species}_{stage}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.chr_fa})" "$(dirname {log})"
        cat {input.fas} > "{output.chr_fa}"
        echo "froze {wildcards.species} ({wildcards.stage}): $(grep -c '^>' {output.chr_fa}) contigs from {input.fas}" > "{log}"
        """
