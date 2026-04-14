# reproducible_phd

A Snakemake-based workflow repository for reproducible comparative genomics analyses.

## Current scope
- Initial HiFi assembly with hifiasm
- Standardized initial assembly outputs
- Planned QC modules:
  - assembly-stats
  - compleasm
  - BlobToolKit
  - GenomeScope
  - optional Merqury

## Structure
- `config/` metadata and workflow configuration
- `workflow/` Snakemake workflow
- `envs/` per-tool conda environment YAMLs
- `scripts/` helper scripts
- `results/` generated outputs (not tracked)
- `logs/` generated logs (not tracked)
- `benchmarks/` generated benchmark files (not tracked)
