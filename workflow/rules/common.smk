import pandas as pd

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
