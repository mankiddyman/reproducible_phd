# Methylation BAM inventory

Goal: locate the source BAM (with MM/ML methylation tags) for each FASTQ used in
hifiasm assembly. Required for downstream comparative epigenomics.

## How to verify a BAM has methylation tags

```bash
samtools view <bam> | head -1 | tr '\t' '\n' | grep -E "^M[MmL]:|^Mm:|^Ml:"
```
Should see something like `MM:Z:C+m,5,12,...` and `ML:B:C,180,200,...`. If
present → has methylation. If absent → BAM was generated without modification
calling and is no useful for methylation analysis.

## Search hints

```bash
# Find any BAM with a specific PacBio run identifier
find /biodata -name "*<run_id>*.bam" 2>/dev/null

# Find any BAM that contains a specific cell ID in its header
samtools view -H <some.bam> | grep "@RG"

# Find all big BAMs in /biodata under marques
find /biodata/dep_mercier/grp_marques -name "*.bam" -size +1G 2>/dev/null
```

---

## Inventory (mark BAM path or NOT FOUND as you locate each)

### Drosera_paradoxa (ploidy=2)

- [x] `6298.A.fastq.gz`
  - dir: `/biodata/dep_mercier/grp_marques/marques/Hologen/Drosera/Drosera_paradoxa_6298A/`
  - BAM: `6298.A.bam` (14.2GB) ✓ same dir
  - methylation tags verified: TODO

- [ ] `6298_A.run801.hifi_reads.fastq.gz`
  - dir: `/biodata/dep_mercier/grp_marques/marques/Hologen/Drosera/Drosera_paradoxa_6298A/`
  - probable BAM: `6298.A.bam` (might be the same source as above? verify by header)
  - BAM path: ___
  - methylation tags verified: ___

- [ ] `r84128_20231119.fastq.gz`
  - dir: `/biodata/dep_mercier/grp_marques/marques/Hologen/Drosera/Drosera_paradoxa_6298A/6036_B/01.Data_result/`
  - BAM path: ___
  - methylation tags verified: ___

### Drosera_binata (ploidy=2)

- [ ] `Drosera_binata.fastq.gz`
  - dir: `/biodata/.../23022LRa012_5546D_75pM-Cell1/`
  - probable BAM: `with_5mC.bam` (9.5GB) — likely same data, renamed
  - BAM path: ___
  - methylation tags verified: ___

- [ ] `deduplicated.6718.B.0_0.fastq.gz`
  - dir: `/biodata/.../Drosera_binata_6718B/`
  - probable BAM: `6718.B.bam` (18.2GB) — might be pre-dedup version
  - BAM path: ___
  - methylation tags verified: ___

### Drosera_roseana (ploidy=2)

- [x] `demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.fastq.gz`
  - BAM: `demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam` (13.6GB) ✓
  - methylation tags verified: TODO

- [x] `6298.C.fastq.gz`
  - BAM: `6298.C.bam` (31.5GB) ✓
  - methylation tags verified: TODO

### Drosera_scorpioides (ploidy=4)

- [ ] `6298.B.m84191_240214_103740_s4.6298.B.7_7.fastq.gz`
  - dir: `/biodata/.../Drosera_scorpioides_6298B/`
  - probable BAM: `6298.B.bam` (21.4GB)
  - BAM path: ___
  - methylation tags verified: ___

- [ ] `6298_B.run800.hifi_reads_m84191_240123_151731_s2.6298.B.7_7.fastq.gz`
  - dir: same as above
  - probable BAM: same `6298.B.bam`? verify by run ID in header
  - BAM path: ___
  - methylation tags verified: ___

- [ ] `r84128_20231119.fastq.gz`
  - dir: `/biodata/.../Drosera_scorpioides_6298B/6036_A/01.Data_result/`
  - BAM path: ___
  - methylation tags verified: ___

- [ ] `demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.fastq.gz`
  - dir: `/biodata/.../Drosera_scorpioides_6298B/Drosera_scorpioides/`
  - BAM path: ___
  - methylation tags verified: ___

- [ ] `demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads2.fastq.gz`
  - dir: same as above
  - BAM path: ___
  - methylation tags verified: ___

- [ ] `m64078e_220722_120155.hifi_reads.fastq.gz`
  - dir: same as above
  - BAM path: ___
  - methylation tags verified: ___

### Drosera_aliciae (ploidy=5)

- [ ] `Drosera_aliciae.fastq.gz`
  - dir: `/biodata/.../Drosera_aliciae_5546A/23022LRa010_5546A_85pM-Cell4/`
  - probable BAM: `with_5mC.bam` (19.8GB) — likely same data, renamed
  - BAM path: ___
  - methylation tags verified: ___

### Drosera_tokaiensis (ploidy=6)

- [ ] `23022LRa011_B01_5546C_90pM-Cell2.fastq.gz`
  - dir: `/biodata/.../23022LRa011_B01_5546C_90pM-Cell2/`
  - probable BAM: `with_5mC.bam` (11.6GB)
  - BAM path: ___
  - methylation tags verified: ___

### Drosera_filiformis (ploidy=?)

- [ ] `deduplicated.6718.A.0_0.fastq.gz`
  - dir: `/biodata/.../Drosera_filiformis_6718A/`
  - probable BAM: `6718.A.bam` (12.3GB) — might be pre-dedup
  - BAM path: ___
  - methylation tags verified: ___

---

## Notes on common patterns

- **`*_with_5mC.bam`** filenames strongly suggest CCS was run with methylation
  calling enabled. Likely good for our purposes.
- **`6298.X.bam`** style — these are the merged/processed bams from the data
  provider; check headers for what cells they include.
- **`deduplicated.X.fastq.gz`** — derived from a non-dedup BAM. The original
  BAM may or may not still have all reads. For methylation, even if the BAM has
  more reads than the FASTQ, that's fine (we just need methylation calls).
- **`r84128_*`** — Revio instrument. The hifi.bam should exist somewhere in the
  delivery directory; if a 01.Data_result folder only has FASTQ, look for a
  sibling 02.Analysis or similar.

