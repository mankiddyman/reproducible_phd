rule run_hifiasm:
    input:
        hifi_reads=get_fastqs,
    output:
        done="results/{species}/hifiasm/{species}.done",
        collapsed_gfa="results/{species}/hifiasm/{species}.p_ctg.gfa",
        rutg_gfa="results/{species}/hifiasm/{species}.r_utg.gfa",
        putg_gfa="results/{species}/hifiasm/{species}.p_utg.gfa",
    params:
        out_prefix="results/{species}/hifiasm/{species}",
        hic_flags=lambda wc: hifiasm_hic_flags(wc.species),
        # Scorpioides (autotetraploid): --n-hap 4 for cleaner p_utg graph.
        # Doesn't affect correction/overlap stages, so cached .bin files
        # remain reusable.
        extra=lambda wc: (
            "--n-hap 4" if wc.species == "Drosera_scorpioides"
            else config["hifiasm"].get("extra", "")
        ),
    threads: 32
    resources:
        # Scorpioides is autotetraploid — needs 950GB for haplotype-resolved
        # partition step. .ec.bin and .ovlp.*.bin caches from previous runs
        # survive and let restarts skip ~30h of correction/overlap work.
        mem_mb=lambda wc, attempt: (
            950000 if wc.species == "Drosera_scorpioides"
            else 500000 * attempt
        ),
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

        # p_ctg, r_utg, p_utg are ALWAYS produced — link them.
        ln -srf "${{prefix}}.${{mode_infix}}.p_ctg.gfa"       "${{prefix}}.p_ctg.gfa"
        ln -srf "${{prefix}}.${{mode_infix}}.r_utg.gfa"       "${{prefix}}.r_utg.gfa"
        ln -srf "${{prefix}}.${{mode_infix}}.p_utg.gfa"       "${{prefix}}.p_utg.gfa"

        # Symlink every haplotype hifiasm produced (hap1..hapN). N depends on
        # --n-hap: diploids get 2, scorpioides (--n-hap 4) gets 4. The .done
        # sentinel is the rule's contract; downstream standardize_initial keys
        # off it and constructs per-hap source paths via hifiasm_role_gfa().
        shopt -s nullglob
        for src in "${{prefix}}.${{mode_infix}}".hap*.p_ctg.gfa; do
            hap=$(basename "$src" | sed -E 's/.*\.(hap[0-9]+)\.p_ctg\.gfa/\1/')
            ln -srf "$src" "${{prefix}}.${{hap}}.p_ctg.gfa"
        done
        shopt -u nullglob


        touch {output.done}
        """

rule standardize_initial:
    """Symlink one hifiasm GFA to a standardized per-role path.
    One job per (species, role); role ∈ {collapsed, hapN, r_utg, p_utg}.
    Depends on the hifiasm .done sentinel, so all hap outputs are guaranteed
    present regardless of how many haplotypes --n-hap produced.
    """
    input:
        done="results/{species}/hifiasm/{species}.done",
    output:
        gfa="results/{species}/assembly/initial/{role}/{species}.gfa",
    wildcard_constraints:
        role=r"collapsed|hap\d+|r_utg|p_utg",
    params:
        src=lambda wc: hifiasm_role_gfa(wc.species, wc.role),
    log:
        "logs/assembly_standardization/{species}_{role}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gfa})" logs/assembly_standardization
        ln -sf "$(realpath {params.src})" {output.gfa}
        echo "Standardized {wildcards.role} for {wildcards.species}" > {log}
        """


rule gfa_to_fasta_initial:
    """Convert a standardized GFA to FASTA. One job per (species, role)."""
    input:
        gfa="results/{species}/assembly/initial/{role}/{species}.gfa",
    output:
        fa="results/{species}/assembly/initial/{role}/{species}.fa",
    wildcard_constraints:
        role=r"collapsed|hap\d+|r_utg|p_utg",
    log:
        "logs/gfa_to_fasta/{species}_{role}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p logs/gfa_to_fasta
        awk '/^S/{{print ">"$2"\n"$3}}' {input.gfa} > {output.fa}
        echo "Converted {wildcards.role} GFA→FASTA for {wildcards.species}" > {log}
        """
