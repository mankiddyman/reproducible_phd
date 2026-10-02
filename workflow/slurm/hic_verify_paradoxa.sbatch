#!/bin/bash
#SBATCH --job-name=hicv_par
#SBATCH --cpus-per-task=64
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --output=logs/hic_verify_paradoxa_%j.out
# Hi-C check of the D. paradoxa joins: the first NPAIRS read pairs of one deep library,
# each mate mapped on its own to the final dual assembly (MAPQ >= 30), paired, counted, drawn.
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
MM=/opt/share/software/bin/minimap2
PY=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/python
REF=results/Drosera_paradoxa/assembly_final/pass2/Drosera_paradoxa_chr.fa
OUT=results/Drosera_paradoxa/qc/hic_verify
W=$OUT/work
NPAIRS=60000000
LIB=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/Files/hic_reads/reads_deep_22_jan_2025
mkdir -p $W
R1=$(ls $LIB/*R1*.fastq.gz | head -1)
R2=$LIB/$(basename $R1 | sed 's/R1/R2/')
[ -s "$R2" ] || R2=$(ls $LIB/*R2*.fastq.gz | head -1)
echo "R1 $R1"; echo "R2 $R2"; date
echo "== index (single part, -I 10G)"
[ -s $W/ref.sr.mmi ] || $MM -x sr -I 10G -t 64 -d $W/ref.sr.mmi $REF 2> $W/index.log
date; echo "== map the first $NPAIRS pairs, each mate on its own, keep MAPQ >= 30"
for m in 1 2; do
  if [ $m = 1 ]; then F=$R1; else F=$R2; fi
  ( zcat $F | head -n $((4 * NPAIRS)) \
      | $MM -x sr --secondary=no -t 30 $W/ref.sr.mmi - 2> $W/map_r$m.log \
      | awk 'BEGIN {OFS = "\t"} $1 != prev {prev = $1; if ($12 >= 30) {n = $1; sub(/\/[12]$/, "", n); print n, $6, int(($8 + $9) / 2)}}' \
      | LC_ALL=C sort -k1,1 -S 30G --parallel 4 -T $W > $W/r$m.tsv ) &
done
wait
date; echo "== pair the mates"
LC_ALL=C join -t $'\t' $W/r1.tsv $W/r2.tsv > $W/pairs.tsv
wc -l $W/r1.tsv $W/r2.tsv $W/pairs.tsv
date; echo "== count and draw"
$PY workflow/scripts/hic_verify.py $REF $W/pairs.tsv $OUT \
    --agp results/Drosera_paradoxa/assembly_final/pass2/Drosera_paradoxa.frozen.FINAL.agp \
    --ragtag results/Drosera_paradoxa/scaffolding/pass2/ragtag/ragtag_output/ragtag.scaffold.agp \
    --pollen /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk/qc/linkage/Dparadoxa_std
date
