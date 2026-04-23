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
