#!/bin/bash
#SBATCH --job-name=hic_remap_par
#SBATCH --cpus-per-task=96
#SBATCH --mem=400G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/hic_remap_paradoxa_%j.out
# Every D. paradoxa Hi-C read pair under Files/hic_reads, mapped onto the final dual assembly
# (each mate on its own, minimap2 -x sr, primary hit, MAPQ >= 1), paired by read name, then
# workflow/scripts/hic_remap_verify.py: library stats, join tests, maps, bait profiles.
cd /netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd
MM=/opt/share/software/bin/minimap2
PY=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/python
REF=results/Drosera_paradoxa/assembly_final/pass2/Drosera_paradoxa_chr.fa
IDX=results/Drosera_paradoxa/qc/hic_verify/work/ref.sr.mmi
HICDIR=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/Files/hic_reads
OUT=results/Drosera_paradoxa/qc/hic_remap
W=$OUT/work
mkdir -p $W
DEC=zcat; command -v pigz > /dev/null && DEC="pigz -dc -p 8"
date; echo "== index"
[ -s $IDX ] || $MM -x sr -I 10G -t 96 -d $IDX $REF 2> $W/index.log
echo "== Hi-C file pairs under $HICDIR"
$PY - "$HICDIR" > $W/libraries.tsv << 'ENDPY'
import os, re, sys
root = sys.argv[1]
libs, unused = {}, []
for d, _, fs in sorted(os.walk(root)):
    for f in sorted(fs):
        if not re.search(r"\.f(ast)?q\.gz$", f):
            continue
        m1 = re.match(r"^(.*[._])R1([._]\d+)?(\.f(?:ast)?q\.gz)$", f)
        m2 = re.match(r"^(.*[._])1(\.f(?:ast)?q\.gz)$", f)
        if m1:
            mate = m1.group(1) + "R2" + (m1.group(2) or "") + m1.group(3)
        elif m2:
            mate = m2.group(1) + "2" + m2.group(2)
        else:
            if not re.search(r"([._]R2([._]\d+)?|[._]2)\.f(ast)?q\.gz$", f):
                unused.append(os.path.join(d, f))
            continue
        if os.path.exists(os.path.join(d, mate)):
            lib = os.path.relpath(d, root)
            libs.setdefault(lib, ([], []))
            libs[lib][0].append(os.path.join(d, f))
            libs[lib][1].append(os.path.join(d, mate))
        else:
            unused.append(os.path.join(d, f))
for lib, (r1, r2) in libs.items():
    print("%s\t%s\t%s" % (lib, ";".join(r1), ";".join(r2)))
for u in unused:
    print("   not used (no mate): %s" % u, file=sys.stderr)
ENDPY
awk -F'\t' '{n = split($2, a, ";"); print "  " $1 ": " n " file pairs"}' $W/libraries.tsv
[ -s $W/libraries.tsv ] || { echo "no Hi-C file pairs found"; exit 1; }
: > $W/pairs_all.tsv
while IFS=$'\t' read -r lib r1s r2s; do
  tag=${lib//\//_}
  date; echo "== $lib"
  for m in 1 2; do
    if [ $m = 1 ]; then FS=${r1s//;/ }; else FS=${r2s//;/ }; fi
    ( $DEC $FS | $MM -x sr --secondary=no -t 44 $IDX - 2> $W/map_${tag}_r$m.log \
        | awk 'BEGIN {OFS = "\t"} $1 != prev {prev = $1; if ($12 >= 1) {n = $1; sub(/\/[12]$/, "", n); print n, $6, int(($8 + $9) / 2), $12}}' \
        | LC_ALL=C sort -k1,1 -S 80G --parallel 8 -T $W > $W/${tag}.r$m.tsv ) < /dev/null &
  done
  wait
  LC_ALL=C join -t $'\t' $W/${tag}.r1.tsv $W/${tag}.r2.tsv | awk -v L="$lib" 'BEGIN {OFS = "\t"} {print $0, L}' >> $W/pairs_all.tsv
  echo "   mate 1 placed: $(wc -l < $W/${tag}.r1.tsv)   mate 2 placed: $(wc -l < $W/${tag}.r2.tsv)   pairs so far: $(wc -l < $W/pairs_all.tsv)"
  rm -f $W/${tag}.r1.tsv $W/${tag}.r2.tsv
done < $W/libraries.tsv
date; echo "== count, test, draw"
$PY workflow/scripts/hic_remap_verify.py $REF $W/pairs_all.tsv $OUT
date
