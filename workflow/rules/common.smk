import pandas as pd
from pathlib import Path
import yaml

species_df = pd.read_csv(config["species_table"]).fillna("")

species_df = species_df[
    species_df["run_hifiasm"].astype(str).str.strip().str.lower() == "yes"
].copy()

SPECIES = species_df["species_id"].tolist()
species_df = species_df.set_index("species_id", drop=False)

# Roles evaluated by QC rules.
# r_utg is standardized and fasta-converted (retained for polyploid ragtag path)
# but not QC'd — it is a raw unitig graph, not a candidate assembly.
QC_ROLES = ["hap1", "hap2", "collapsed"]

def get_fastqs(wildcards):
    raw = species_df.loc[wildcards.species, "genomic_fastq"]
    fastqs = [x.strip() for x in str(raw).split(";") if x.strip()]
    if not fastqs:
        raise ValueError(f"No FASTQ files found for species {wildcards.species}")
    return fastqs


def load_decisions(species: str) -> dict:
    """Load the per-species decisions.yaml produced by the compile_decisions
    checkpoint. Used by downstream rules (future stages) that need to branch
    on assembly_strategy, is_heterozygous, etc.

    Raises FileNotFoundError if called before the checkpoint has run — in a
    Snakemake checkpoint-aware context, callers must reach this via
    checkpoints.compile_decisions.get(species=...).output.yaml, not directly.
    """
    path = Path(f"results/{species}/decisions.yaml")
    if not path.is_file():
        raise FileNotFoundError(
            f"decisions.yaml for {species} not yet produced: {path}"
        )
    with path.open() as f:
        return yaml.safe_load(f)


# ---- HiC library helpers (consumed by run_hifiasm) -------------------------
import glob as _glob

def _load_hic_libraries(hic_csv: str = "config/hic_libraries.csv") -> pd.DataFrame:
    """Load the per-species hic library manifest. Empty DataFrame if missing."""
    path = Path(hic_csv)
    if not path.is_file():
        return pd.DataFrame(columns=["species_id", "library_id", "r1_glob", "r2_glob"])
    return pd.read_csv(path).fillna("")


_HIC_LIBS = _load_hic_libraries()


def get_hic_pairs(species: str) -> list[tuple[list[str], list[str]]]:
    """Return [(r1_files, r2_files), ...] per hic library for the species.

    Each tuple corresponds to one CSV row; r1_files / r2_files are
    sorted glob-expansions of the r1_glob / r2_glob fields.
    Pairs are skipped if either side has zero matches.
    Returns [] if the species has no hic data.
    """
    if _HIC_LIBS.empty:
        return []
    subset = _HIC_LIBS[_HIC_LIBS["species_id"] == species]
    pairs = []
    for _, row in subset.iterrows():
        r1 = sorted(_glob.glob(str(row["r1_glob"])))
        r2 = sorted(_glob.glob(str(row["r2_glob"])))
        if r1 and r2:
            pairs.append((r1, r2))
    return pairs


def has_hic(species: str) -> bool:
    """True if species has at least one resolvable hic library."""
    return len(get_hic_pairs(species)) > 0


def hic_r1_files(species: str) -> list[str]:
    """All R1 hic files for a species, flattened across libraries."""
    return [f for pair in get_hic_pairs(species) for f in pair[0]]


def hic_r2_files(species: str) -> list[str]:
    """All R2 hic files for a species, flattened across libraries."""
    return [f for pair in get_hic_pairs(species) for f in pair[1]]


def hifiasm_hic_flags(species: str) -> str:
    """Return '--h1 r1_a.fq,r1_b.fq --h2 r2_a.fq,r2_b.fq' or '' if no hic.

    Comma-separated, no spaces — hifiasm's expected format.
    """
    r1 = hic_r1_files(species)
    r2 = hic_r2_files(species)
    if not (r1 and r2):
        return ""
    return f"--h1 {','.join(r1)} --h2 {','.join(r2)}"


# Roles evaluated by QC rules.
# r_utg/p_utg are graph outputs; QC them on hap1/hap2/collapsed which are linear
QC_ROLES = ["hap1", "hap2", "collapsed"]

# Assembly outputs that exist for every species (post-hifiasm + standardize)
ALL_ASSEMBLY_TYPES = ["collapsed", "hap1", "hap2", "r_utg", "p_utg"]


def decontamination_targets(species: str) -> list[str]:
    """Return which assembly file(s) blobtoolkit should be run on for a species.

    Logic:
      - ploidy <= 2 (or unknown): use hap1, hap2 (hifiasm's haplotype split is
        meaningful for diploids).
      - ploidy > 2 (polyploids): use p_utg (hifiasm's hap1/hap2 split forces
        4+ alleles into 2 bins and is biologically meaningless; the processed
        unitig graph retains the full information for downstream phasing).
    """
    ploidy_raw = species_df.loc[species, "exp_ploidy"]
    try:
        ploidy = int(ploidy_raw) if str(ploidy_raw).strip() else 2
    except (ValueError, TypeError):
        ploidy = 2  # safe default
    if ploidy > 2:
        return ["p_utg"]
    return ["hap1", "hap2"]


def all_decontamination_done(species_list: list[str]) -> list[str]:
    """Expand .done files across (species, assembly) per-species decontamination targets."""
    out = []
    for sp in species_list:
        for tgt in decontamination_targets(sp):
            out.append(f"results/{sp}/blobtoolkit/initial/{tgt}/.done")
    return out


# ---- Decontamination target selection --------------------------------------

def decontamination_targets(species: str) -> list[str]:
    """Return which assembly file(s) blobtoolkit should be run on for a species.

    Logic:
      - ploidy <= 2 (or unknown): use hap1, hap2 (hifiasm's haplotype split is
        meaningful for diploids).
      - ploidy > 2 (polyploids): use p_utg (hifiasm's hap1/hap2 split forces
        4+ alleles into 2 bins and is biologically meaningless; the processed
        unitig graph retains full information for downstream phasing).
    """
    ploidy_raw = species_df.loc[species, "exp_ploidy"]
    try:
        ploidy = int(ploidy_raw) if str(ploidy_raw).strip() else 2
    except (ValueError, TypeError):
        ploidy = 2  # safe default
    if ploidy > 2:
        return ["p_utg"]
    return ["hap1", "hap2"]


def all_decontamination_done(species_list: list[str]) -> list[str]:
    """Expand .done files across (species, assembly) per-species decontamination targets."""
    return [
        f"results/{sp}/blobtoolkit/initial/{tgt}/.done"
        for sp in species_list
        for tgt in decontamination_targets(sp)
    ]


# ---- Decontamination config & helpers --------------------------------------

def load_decontamination_config(path: str = "config/decontamination.yaml") -> dict:
    """Load the decontamination filter config. Returns dict with 'defaults' and
    'species_overrides'."""
    p = Path(path)
    if not p.is_file():
        return {"defaults": {"keep_classes": ["Magnoliopsida", "no-hit"]},
                "species_overrides": {}}
    with p.open() as f:
        cfg = yaml.safe_load(f) or {}
    cfg.setdefault("defaults", {"keep_classes": ["Magnoliopsida", "no-hit"]})
    cfg.setdefault("species_overrides", {})
    if cfg["species_overrides"] is None:
        cfg["species_overrides"] = {}
    return cfg


DECONTAMINATION_CONFIG = load_decontamination_config()


def decontamination_targets_files(species_list: list[str], extension: str) -> list[str]:
    """Expand per-(species, assembly) decontamination output files.
    
    For each species, uses decontamination_targets(species) to pick the right
    assembly types (hap1/hap2 for diploids, p_utg for polyploids).
    """
    return [
        f"results/{sp}/assembly/decontaminated/{tgt}/{sp}.{extension}"
        for sp in species_list
        for tgt in decontamination_targets(sp)
    ]


def decontamination_qc_files(species_list: list[str], qc_type: str) -> list[str]:
    """Expand per-(species, assembly) QC files for decontaminated assemblies.
    
    qc_type: 'assembly_stats' or 'compleasm'
    """
    paths = []
    for sp in species_list:
        for tgt in decontamination_targets(sp):
            if qc_type == "assembly_stats":
                paths.append(f"results/{sp}/qc/assembly_stats/decontaminated/{tgt}.tsv")
            elif qc_type == "compleasm":
                paths.append(f"results/{sp}/qc/compleasm/decontaminated/{tgt}/summary.txt")
    return paths
