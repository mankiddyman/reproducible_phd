#!/bin/bash
# DR05b — six concatenated trees + ASTRAL.
#   all13_* : 13 tips, both subgenomes of every species. Places regia
#             relative to Dionaea -- the question DR03b could not settle.
#   A_*/B_* : 7 tips each. Same species, two independent partitions.
#             Same topology = the single-hybridisation prediction.
#   *_stable: mislabel sensitivity. 89% of segments never flip in 1000
#             bootstraps, so these differ by ~30 loci and should match.
set -u
IQ=/opt/share/software/bin/iqtree2
T=64
cd "$SUBG_BASE"
say() { printf "\n\033[1m=== %s ===\033[0m  [%s]\n" "$1" "$(date +%H:%M:%S)"; }

say "0. environment"
$IQ --version | head -2
echo "threads: $T"
ls -la DR/tree/*.phy
echo "gene alignments: $(ls DR/tree/genes/*.fna 2>/dev/null | wc -l)"

for M in all13_full all13_stable A_full B_full A_stable B_stable; do
  say "concatenated tree: $M"
  $IQ -s DR/tree/$M.phy -m GTR+F+I+G4 -B 1000 -alrt 1000 \
      -T $T --prefix DR/tree/concat/$M -redo 2>&1 | grep -E "Best-fit|BEST SCORE|Total wall"
  echo "  --- site concordance factors ---"
  $IQ -t DR/tree/concat/$M.treefile -s DR/tree/$M.phy \
      --scf 100 -T $T --prefix DR/tree/concat/${M}_scf -redo 2>&1 | grep -E "Total wall"
  echo "  TREE: $(cat DR/tree/concat/$M.treefile)"
done

say "gene trees for ASTRAL"
mkdir -p DR/tree/genetrees
ls DR/tree/genes/*.fna > /tmp/gtlist.txt
echo "  $(wc -l < /tmp/gtlist.txt) alignments at -j $T"
time parallel -j $T --bar "b=\$(basename {} .fna); \
  $IQ -s {} -m GTR+G -T 1 --prefix DR/tree/genetrees/\$b -redo > /dev/null 2>&1" \
  :::: /tmp/gtlist.txt
echo "  built: $(ls DR/tree/genetrees/*.treefile 2>/dev/null | wc -l)"

say "ASTRAL"
cat DR/tree/genetrees/*.treefile > DR/tree/astral/genetrees_all.tre 2>/dev/null
echo "  $(wc -l < DR/tree/astral/genetrees_all.tre) gene trees"
if command -v astral >/dev/null 2>&1; then
  astral -i DR/tree/astral/genetrees_all.tre -o DR/tree/astral/astral_all.tre \
    2> DR/tree/astral/astral_all.log
  echo "  ASTRAL: $(cat DR/tree/astral/astral_all.tre)"
  grep -iE "normalized quartet|final quartet" DR/tree/astral/astral_all.log
else
  echo "  astral NOT on PATH. Gene trees are in DR/tree/astral/genetrees_all.tre"
fi

say "SUMMARY"
for M in all13_full all13_stable A_full B_full A_stable B_stable; do
  echo "--- $M ---"; cat DR/tree/concat/$M.treefile 2>/dev/null; echo
done
cat <<'TXT'

READ IN THIS ORDER
  1. all13_full: is regia_A inside the Drosera A clade, or sister to
     Dionaea_A? Check regia_B independently -- both should agree.
     THE question DR03b could not settle.
  2. A_full vs B_full: same Drosera topology?
     same   -> one hybridisation, both subgenomes share a history
     differ -> the subgenomes have different histories. Big result.
  3. *_stable vs *_full: identical? then the 6-17% mislabel rate does
     not affect the conclusion and no threshold was ever needed.
  4. sCF in DR/tree/concat/*_scf.cf.branch. Bootstrap saturates at 100
     on 2 Mb of concatenated data; sCF says what fraction of SITES
     actually support each branch. THAT is the number to quote.
  5. ASTRAL vs concatenation: disagreement means real ILS.
  6. capensis_B occupancy is 0.322 -- two-thirds gaps. Its placement is
     the least reliable in the tree.
TXT
