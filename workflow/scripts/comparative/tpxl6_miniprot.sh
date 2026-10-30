#!/usr/bin/env bash
# tpxl6_miniprot.sh -- TPXL6 in every frozen genome by spliced protein-to-genome alignment (miniprot)
#
# Settles audit C1 (is the D. binata TPXL6 copy a pseudogene or an intact reading frame?) without any gene
# model: miniprot aligns the protein across introns and reports frameshifts and in-frame stop codons itself.
# The June calls came from tblastn hits, which split ANY multi-exon gene across reading frames whenever an
# intron length is not a multiple of 3 -- and binata's introns are 15-19 kb of transposon sequence.
#
# Queries (from tpx2/tpx2_family.faa):
#   regia_g4881.t1   D. regia TPXL6 (139 aa)
#   regia_g4900.t1   D. regia member of the TPXL1/TPXL5 group, TPXL6's closest relatives; each genomic locus is
#                    later assigned to whichever query aligns there with the higher score, so a locus of the
#                    closest paralog is not mistaken for TPXL6
# Genomes: every species in config/annotation.csv, frozen chromosome FASTA
#   results/{species}/assembly_final/{frozen_stage}/{species}_chr.fa ; missing ones are reported and skipped.
#   Order: D. binata and D. regia first (the C1 decision), then the rest.
# miniprot settings: defaults (max intron 200 kb) plus --outs=0.3 -p 0.3, so a second haplotype or a degraded
#   copy scoring well below the best hit is still reported.
# Output: results/comparative/holocentricity/tpx2/miniprot/{species}.gff, .log, .done
#   (read by workflow/scripts/comparative/tpxl6_miniprot_summary.R)
#
# usage: bash tpxl6_miniprot.sh [REPO] [THREADS]
REPO=${1:-/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd}
THREADS=${2:-32}
ENV=${PHYLO_ENV:-/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/phylo}
B=$ENV/bin
IN=$REPO/results/comparative/holocentricity/tpx2
OUT=$IN/miniprot
die() { echo "!! $*" >&2; exit 1; }
[ -x "$B/miniprot" ] || die "$B/miniprot missing: create the env from envs/phylo.yaml"
[ -s "$IN/tpx2_family.faa" ] || die "missing $IN/tpx2_family.faa"
[ -s "$REPO/config/annotation.csv" ] || die "missing $REPO/config/annotation.csv"
mkdir -p "$OUT" && cd "$OUT" || die "cannot enter $OUT"

awk 'BEGIN { want["regia_g4881.t1"]; want["regia_g4900.t1"] }
     /^>/ { id = substr($1, 2); keep = (id in want); if (keep) print ">" id; next }
     keep { gsub(/[ \t\r*]/, ""); print }' "$IN/tpx2_family.faa" > queries.faa
if ! grep -q '^>regia_g4900.t1$' queries.faa && [ -s "$IN/tpx2_short.faa" ]; then     # fallback: the short-tree input
  awk '/^>/ { keep = (substr($1, 2) == "regia_g4900.t1"); if (keep) print ">regia_g4900.t1"; next }
       keep { gsub(/[ \t\r*]/, ""); print }' "$IN/tpx2_short.faa" >> queries.faa
fi
[ "$(grep -c '^>' queries.faa)" -eq 2 ] || die "regia_g4881.t1 and regia_g4900.t1 not both found in tpx2_family.faa / tpx2_short.faa"
awk '/^>/ { if (id) printf "   query %s %d aa\n", id, L; id = substr($1, 2); L = 0; next } { L += length($0) }
     END { printf "   query %s %d aa\n", id, L }' queries.faa
echo "   $("$B/miniprot" --version 2>&1 | head -1 | sed 's/^/miniprot /')"

# species -> frozen genome path, D. binata and D. regia first
awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "frozen_stage") c = i; next }
         c && $c != "" { printf "%s\tresults/%s/assembly_final/%s/%s_chr.fa\n", $1, $1, $c, $1 }' \
    "$REPO/config/annotation.csv" \
  | awk -F'\t' '{ o = ($1 == "Drosera_binata") ? 0 : ($1 == "Drosera_regia") ? 1 : 2; print o "\t" $0 }' \
  | sort -k1,1n -s | cut -f2- > genomes.tsv
echo "   $(wc -l < genomes.tsv) genomes listed in config/annotation.csv"

while IFS=$'\t' read -r sp rel; do
  G=$REPO/$rel
  if [ -s "$sp.done" ] && [ -s "$sp.gff" ]; then echo "== $sp: done earlier, skipped"; continue; fi
  if [ ! -s "$G" ]; then echo "== $sp: genome not found ($G), skipped"; continue; fi
  t0=$(date +%s)
  echo "== $sp   $(date '+%F %T')   $G ($(( $(stat -c %s "$G") / 1000000 )) MB)"
  if "$B/miniprot" -t "$THREADS" --gff --outs=0.3 -p 0.3 "$G" queries.faa > "$sp.gff.tmp" 2> "$sp.log"; then
    mv "$sp.gff.tmp" "$sp.gff"
    echo "$G" > "$sp.done"
    echo "   $(grep -c $'\tmRNA\t' "$sp.gff") alignments, $(( ($(date +%s) - t0) / 60 )) min"
  else
    echo "   !! miniprot failed for $sp (see $OUT/$sp.log)"
  fi
done < genomes.tsv
echo "== finished   $(date '+%F %T')"
