# Synteny: GENESPACE (module-provided, v1.3.1) across chromosome-level genomes.
# NB: GENESPACE bundles its OWN OrthoFinder (v2.5.5, module). Deliberately kept
# separate from rules/comparative.smk's OrthoFinder -- versions are incompatible
# and the two must never share an input or output directory.
import pandas as _pd

genespace_df = _pd.read_csv("config/genespace.csv").set_index("species_id", drop=False)
GENESPACE_SPECIES = list(genespace_df.index)


def _gs(species, field):
    return str(genespace_df.loc[species, field]).strip()


rule genespace_prep:
    """One haploid chromosome set per species -> GENESPACE bed/ + peptide/.

    Subsets the final GFF3 to the chromosomes matching chr_regex (config/
    genespace.csv), then derives BOTH the bed and the peptide fasta from that
    same subset so the ids are guaranteed to match -- mismatched ids are the
    usual way GENESPACE fails.

    Keeps ALL genes regardless of RNA support (_rnasupp and _norna alike);
    revisit if synteny blocks look noisy in the TE-rich genomes."""
    input:
        gff="results/{species}/annotation/final/{species}.final.gff3",
        genome=lambda wc: f"results/{wc.species}/assembly_final/{_gs(wc.species,'stage')}/{wc.species}_chr.fa",
    output:
        bed="results/comparative/genespace/bed/{species}.bed",
        pep="results/comparative/genespace/peptide/{species}.fa",
        gff=temp("results/comparative/genespace/tmp/{species}.haploid.gff3"),
    params:
        chr_regex=lambda wc: _gs(wc.species, "chr_regex"),
    log:
        "logs/genespace_prep/{species}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bed})" "$(dirname {output.pep})" \
                 "$(dirname {output.gff})" "$(dirname {log})"
        module load genespace 2>/dev/null || true

        # 1. subset GFF3 to the haploid chromosome set
        awk -F'\t' -v re='{params.chr_regex}' \
            '/^#/ {{print; next}} $1 ~ re' "{input.gff}" > "{output.gff}"
        echo "=== {wildcards.species} regex={params.chr_regex} ===" > "{log}"
        echo "chromosomes kept: $(awk -F'\t' '!/^#/{{print $1}}' {output.gff} | sort -u | tr '\n' ' ')" >> "{log}"
        echo "transcripts kept: $(awk -F'\t' '!/^#/ && $3=="transcript"' {output.gff} | wc -l)" >> "{log}"

        # 2. bed (4 col) and 3. peptide, both from the same subset.
        #    sed strips gffread's '.' intron/exon markers, which break DIAMOND.
        gffread "{output.gff}" --bed -o - 2>> "{log}" | cut -f1-4 > "{output.bed}"
        gffread "{output.gff}" -g "{input.genome}" -y - 2>> "{log}" \
            | sed '/^>/! s/\.//g' > "{output.pep}"

        NB=$(wc -l < "{output.bed}"); NP=$(grep -c '^>' "{output.pep}")
        echo "bed entries: $NB | peptides: $NP" >> "{log}"
        [ "$NB" -eq "$NP" ] || {{ echo "ERROR: bed/peptide count mismatch" >> "{log}"; exit 1; }}
        """


rule genespace_prep_all:
    input:
        expand("results/comparative/genespace/bed/{sp}.bed", sp=GENESPACE_SPECIES),
        expand("results/comparative/genespace/peptide/{sp}.fa", sp=GENESPACE_SPECIES),


rule genespace_run:
    """Run GENESPACE v1.3.1 (module) over the prepped haploid genomes.

    GENESPACE wants ONE working directory containing bed/ and peptide/, and it
    writes its results back into that same directory (orthofinder/, results/,
    riparian/, syntenicHits/). It runs its OWN bundled OrthoFinder v2.5.5 --
    NEVER point this at results/comparative/orthofinder/, the versions are
    incompatible.

    path2mcscanx is read from the loaded module (MCSCANX_HOME/mcscanx_root) and
    falls back to config['synteny']['path2mcscanx'] if the module does not
    export one -- the path moved between module versions (stretchApps/va8443a9
    -> bookwormApps/vb1ca533), so do not hardcode it."""
    input:
        beds=expand("results/comparative/genespace/bed/{sp}.bed", sp=GENESPACE_SPECIES),
        peps=expand("results/comparative/genespace/peptide/{sp}.fa", sp=GENESPACE_SPECIES),
        script="workflow/scripts/comparative/run_genespace.R",
    output:
        riparian=directory("results/comparative/genespace/riparian"),
        results=directory("results/comparative/genespace/results"),
    params:
        wd="results/comparative/genespace",
        mcscanx=config.get("synteny", {}).get("path2mcscanx", ""),
        # genomeIDs in config/genespace.csv row order; ref = first row
        ids=",".join(GENESPACE_SPECIES),
        ref=GENESPACE_SPECIES[0],
        # ploidy in the SAME order as ids -- init_genespace matches positionally
        ploidy=",".join(str(int(genespace_df.loc[sp, "ploidy"])) for sp in GENESPACE_SPECIES),
    threads: 48
    resources:
        mem_mb=200000,
        runtime=2880,
    log:
        "logs/genespace/run.log",
    shell:
        r"""
        set -euo pipefail
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"; mkdir -p "$(dirname $LOG)"
        module load genespace 2>/dev/null || {{ echo "cannot load genespace module" > "$LOG"; exit 1; }}

        MCX="{params.mcscanx}"
        if [ -z "$MCX" ]; then
            MCX="${{MCSCANX_HOME:-}}"
            [ -z "$MCX" ] && MCX=$(dirname "$(command -v MCScanX 2>/dev/null || echo /nonexistent)")
        fi
        [ -d "$MCX" ] || {{ echo "MCScanX path not found: '$MCX' -- set config synteny.path2mcscanx" >> "$LOG"; exit 1; }}
        echo "=== genespace $(date) ===" > "$LOG"
        echo "wd=$REPO_ROOT/{params.wd}"  >> "$LOG"
        echo "mcscanx=$MCX"               >> "$LOG"
        echo "orthofinder=$(command -v orthofinder)" >> "$LOG"

        IDS="{params.ids}"
        REF="{params.ref}"
        PLOIDY="{params.ploidy}"
        echo "ref=$REF ids=$IDS ploidy=$PLOIDY" >> "$LOG"
        Rscript "{input.script}" "$REPO_ROOT/{params.wd}" "$MCX" {threads} "$REF" "$IDS" "$PLOIDY" >> "$LOG" 2>&1
        """
