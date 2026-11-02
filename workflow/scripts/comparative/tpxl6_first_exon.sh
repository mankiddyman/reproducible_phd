#!/usr/bin/env bash
# tpxl6_first_exon.sh -- find the first exon and the promoter of a TPXL6 copy whose miniprot alignment does not start at
# amino acid 1, and measure CpG methylation there
#
# D. binata hap1 starts at amino acid 12: miniprot did not find its first exon (aa 1-11, ~35 bp), so x = 0 of its row
# in the methylation figures is exon 2, and its promoter was never measured. For every such copy of the species:
#   1 locate  (tpxl6_first_exon.R locate) take the start of the complete copy of the same species (the copy aligned
#             from amino acid 1): 3 kb upstream of its first exon, the first exon, intron 1 and exon 2, with CARP
#             repeats masked. Align it (minimap2 -x asm20) to the 200 kb upstream of the incomplete copy, and also scan
#             that DNA base by base for the first exon itself (too short for minimap2 on its own). The start of the
#             incomplete copy is where that first exon matches, or where the complete copy's start codon aligns. If
#             neither is found, the complete copy's start is aligned to the whole genome, to see whether it lies elsewhere
#   2 reads   if the HiFi reads already collected by tpxl6_methyl_locus.sh cover the found start (>= 80% of its CpGs from
#             -2 kb to +0.5 kb, median >= 5 reads, from the BAMs whose TEs are methylated), they are used. Otherwise
#             those BAMs are scanned again for reads carrying 31-mers of the found start +- 5 kb (repeats masked), and
#             the reads are aligned and piled up as in tpxl6_methyl_locus.sh (~10 min for the 9.5 GB D. binata BAM)
#   3 report  (tpxl6_first_exon.R report) CpG methylation around the found start in the windows of
#             tpxl6_methyl_tss_check.R, next to the complete copy and the related gene (same BAM), and a figure: how
#             the complete copy's start aligns to the incomplete copy, and methylation -3 kb to +1.5 kb around each start
# Output: results/comparative/holocentricity/tpx2/methyl/{species}/first_exon/ (anchor.tsv, blocks.tsv, exon1_hits.tsv,
#         and the new pileup if one was needed) and results/comparative/holocentricity/figures/tpxl6_first_exon_{species}.png/.pdf/.tsv
# Needs: tpxl6_methyl_locus.sh and plot_tpxl6_methyl_locus.R run first (regions, pileups, methyl_bam_check.tsv)
# usage: bash tpxl6_first_exon.sh [REPO] [SPECIES]   (default species: Drosera_binata)
REPO=${1:-/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd}
SP=${2:-Drosera_binata}
ENV=${METHYL_ENV:-/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/methyl}
B=$ENV/bin; SAM=$B/samtools
RS=${RSCRIPT:-/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/Rscript}
[ -x "$RS" ] || RS=$(command -v Rscript)
K=31; STRIDE=4; T=${THREADS:-16}
H=$REPO/results/comparative/holocentricity
MP=$H/tpx2/miniprot; D=$H/tpx2/methyl/$SP; O=$D/first_exon
SDIR=$(cd "$(dirname "$0")" && pwd)
die() { echo "!! $*" >&2; exit 1; }
for t in samtools minimap2 aligned_bam_to_cpg_scores; do [ -x "$B/$t" ] || die "$B/$t missing: create the env from envs/methyl.yaml"; done
[ -n "$RS" ] || die "no Rscript found"
[ -s "$D/regions.tsv" ] || die "$D/regions.tsv missing: run tpxl6_methyl_locus.sh first"
[ -s "$H/tpx2/methyl/methyl_bam_check.tsv" ] || die "$H/tpx2/methyl/methyl_bam_check.tsv missing: run plot_tpxl6_methyl_locus.R first"
cd "$REPO" || exit 1
mkdir -p "$O" "$H/figures"
G=$(awk -F'\t' -v sp="$SP" '$1 == sp { print $2 }' "$MP/genomes.tsv")
[ -s "$G" ] || die "genome of $SP not found in $MP/genomes.tsv"
[ -s "$G.fai" ] || $SAM faidx "$G"

# K-mers starting every STRIDE bp of the sequences on stdin (one per line), both strands; k-mers with N (masked) or
# one base > 70% dropped (as in tpxl6_methyl_locus.sh)
kmers() {
  awk -v k="$K" -v st="$STRIDE" '{ s = toupper($0); for (i = 1; i + k - 1 <= length(s); i += st) { m = substr(s, i, k); if (m !~ /N/) print m } }' \
    | awk '{ print; r = ""; for (i = length($0); i > 0; i--) { b = substr($0, i, 1)
             r = r (b == "A" ? "T" : b == "C" ? "G" : b == "G" ? "C" : b == "T" ? "A" : "N") }; print r }' \
    | awk -v k="$K" '$0 !~ /N/ { m = 0; for (b in c) delete c[b]; for (i = 1; i <= k; i++) c[substr($0, i, 1)]++
                                 for (b in c) if (c[b] > m) m = c[b]; if (m <= 0.7 * k) print }' | LC_ALL=C sort -u
}

# one BAM: reads carrying a k-mer of the found start -> $O/reads.$n.fq (MM/ML, HP = BAM number), as in tpxl6_methyl_locus.sh
scan_bam() {
  local n=$1 bam=$2 hp
  if [ -s "$O/scan.$n.tsv" ] && [ -e "$O/reads.$n.fq" ]; then echo "   BAM $n: already scanned"; return 0; fi
  [ -s "$bam" ] || { echo "   !! BAM $n not found: $bam"; return 1; }
  hp=$(( n <= 2 ? n : 0 ))
  $SAM view -@ 4 --keep-tag MM,ML,Mm,Ml "$bam" | LC_ALL=C grep -F -f "$O/regions.k" | awk -F'\t' -v hp="$hp" -v n="$n" '
    BEGIN { cp["A"] = "T"; cp["C"] = "G"; cp["G"] = "C"; cp["T"] = "A"; cp["N"] = "N" }
    int($2 / 256) % 2 || int($2 / 2048) % 2 { next }                    # secondary / supplementary
    { seq = $10; q = $11
      if (q == "*") { q = sprintf("%" length(seq) "s", ""); gsub(/ /, "~", q) }     # no stored qualities
      if (int($2 / 16) % 2) {                                           # aligned on the reverse strand: restore the read
        r = ""; for (i = length(seq); i > 0; i--) r = r cp[substr(seq, i, 1)]; seq = r
        r = ""; for (i = length(q); i > 0; i--) r = r substr(q, i, 1); q = r }
      tags = ""
      for (i = 12; i <= NF; i++) { t = $i; sub(/^Mm:/, "MM:", t); sub(/^Ml:/, "ML:", t); if (t ~ /^(MM|ML):/) tags = tags "\t" t }
      print "@" $1 "_b" n tags "\tHP:i:" hp; print seq; print "+"; print q }' > "$O/reads.$n.fq.part" \
    && mv "$O/reads.$n.fq.part" "$O/reads.$n.fq" || { echo "   !! BAM $n: scan failed"; return 1; }
  echo -e "$n\t$bam\t$(( $(wc -l < "$O/reads.$n.fq") / 4 ))" > "$O/scan.$n.tsv"
  echo "   BAM $n: $(( $(wc -l < "$O/reads.$n.fq") / 4 )) reads carry a k-mer of the found start   ($(basename "$bam"), done $(date '+%T'))"
}

# 1 locate
echo "== $SP: locate the missing first exon   $(date '+%F %T')"
rm -f "$O/anchor.tsv"
"$RS" "$SDIR/tpxl6_first_exon.R" locate "$REPO" "$SP" "$G" "$O" "$B" "$T" || die "locate step failed"
[ -s "$O/anchor.tsv" ] || die "locate step wrote no $O/anchor.tsv"

# 2 reads: rescan only when the existing pileup covers the found start too thinly
need=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } $c["need_rescan"] == "TRUE" { n++ } END { print n + 0 }' "$O/anchor.tsv")
if [ "$need" -gt 0 ]; then
  echo "== $SP: the reads collected earlier cover the found start too thinly; scanning the usable BAMs again   $(date '+%F %T')"
  CARP=results/$SP/repeats/carp/Repeat_Annotation_Unified.gff3
  while IFS=$'\t' read -r ctg lo hi; do
    echo "   region $ctg:$lo-$hi"
    if [ -s "$CARP" ]; then
      awk -F'\t' -v c="$ctg" -v lo="$lo" -v hi="$hi" '!/^#/ && $1 == c && $9 !~ /(^|;)Parent=/ && $4 <= hi && $5 >= lo' "$CARP" > "$O/repeats.tmp"
    else : > "$O/repeats.tmp"; fi
    $SAM faidx "$G" "$ctg:$lo-$hi" | grep -v '^>' | tr -d '\n' | awk -v c="$ctg" -v lo="$lo" -v rf="$O/repeats.tmp" '
      FILENAME == rf { if ($1 == c) { n++; A[n] = $4; Z[n] = $5 }; next }
      { s = $0; for (i = 1; i <= n; i++) { a = A[i] - lo + 1; b = Z[i] - lo + 1; if (a < 1) a = 1; if (b > length(s)) b = length(s)
          if (a <= b) { x = sprintf("%" (b - a + 1) "s", ""); gsub(/ /, "N", x); s = substr(s, 1, a - 1) x substr(s, b + 1) } }
        print s }' "$O/repeats.tmp" -
  done < <(tail -n +2 "$O/rescan_region.tsv") | kmers > "$O/regions.k.new"
  rm -f "$O/repeats.tmp"
  if ! cmp -s "$O/regions.k.new" "$O/regions.k"; then rm -f "$O"/scan.*.tsv "$O"/reads.*.fq; fi
  mv "$O/regions.k.new" "$O/regions.k"
  echo "   $(wc -l < "$O/regions.k") k-mers (both strands, one every $STRIDE bp) from the unmasked part"
  tail -n +2 "$O/usable_bams.tsv" > "$O/.jobs"
  while IFS=$'\t' read -r n bam; do scan_bam "$n" "$bam" < /dev/null & done < "$O/.jobs"
  wait
  rm -f "$O/.jobs"
  ls "$O"/reads.*.fq 2>/dev/null | sort -V | xargs -r cat > "$O/reads.fq"
  [ -s "$O/reads.fq" ] || die "no reads found for the found start"
  "$B/minimap2" -ax map-hifi -y -Y --secondary=no -t "$T" "$G" "$O/reads.fq" 2> "$O/minimap2.log" \
    | $SAM sort -@ 4 -o "$O/aligned.bam" - 2>> "$O/minimap2.log" && $SAM index "$O/aligned.bam" || die "alignment failed, see $O/minimap2.log"
  echo "   aligned: $($SAM view -c -F 0x904 "$O/aligned.bam") primary alignments, $($SAM view -c -F 0x904 -q 1 "$O/aligned.bam") with MAPQ >= 1"
  "$B/aligned_bam_to_cpg_scores" --bam "$O/aligned.bam" --output-prefix "$O/$SP" --pileup-mode count \
    --modsites-mode reference --ref "$G" --min-coverage 3 --min-mapq 1 --threads "$T" > "$O/pb_cpg.log" 2>&1 \
    || die "pb-CpG-tools failed, see $O/pb_cpg.log"
  rm -f "$O/reads.fq"
  echo "   pileup: $(ls "$O/$SP".*.bed.gz | xargs -n1 basename | tr '\n' ' ')"
else
  rm -f "$O/$SP".* "$O"/aligned.bam*                     # an older rescan would otherwise be read by the report
  echo "== $SP: the reads collected earlier cover the found start; no new scan needed"
fi

# 3 report
echo "== $SP: methylation at the found start   $(date '+%F %T')"
"$RS" "$SDIR/tpxl6_first_exon.R" report "$REPO" "$SP" "$G" "$O" "$B" "$T" "$H/figures/tpxl6_first_exon_$SP" 2>&1 | grep -v 'Fontconfig warning'
echo "== done $(date '+%F %T')"
