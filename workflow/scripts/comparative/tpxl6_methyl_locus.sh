#!/usr/bin/env bash
# tpxl6_methyl_locus.sh -- CpG methylation over every TPXL6 copy and one expressed relative, from PacBio HiFi 5mC
# tags (MM/ML), without aligning the whole read set
#
# For each species given (default: Drosera_binata Drosera_paradoxa), with its HiFi BAMs from config/hifi_5mc.csv
# (column hifi_bam, ';'-separated, as in config/species.csv):
#   1 regions  every TPXL6 copy (assigned TPXL6 in tpx2/miniprot/tpxl6_miniprot_loci.tsv) and the best TPXL1/5-group
#              locus (the control gene: expressed in leaves in every species with RNA), each +- FLANK bp
#   2 k-mers   31-mers of those regions with every CARP repeat masked, so TE k-mers do not pull in reads from the TE's
#              other copies; both strands; low-complexity k-mers removed
#   3 reads    every read of each BAM that carries one of those k-mers, kept with its 5mC tags. A BAM is skipped (and
#              named) when fewer than half of its first 1000 reads carry 5mC tags. Each read is labelled with the
#              number of its BAM in the HP tag (1 or 2; 0 from a third BAM on), so pb-CpG-tools also writes one
#              pileup per BAM ("hap1" = BAM 1, "hap2" = BAM 2): the plot script checks each BAM separately, because
#              an amplified (PCR) library carries 5mC tags but no methylation
#   4 align    minimap2 -ax map-hifi -y (tags kept) -Y (supplementary alignments soft-clipped, as pb-CpG-tools needs)
#              to the whole frozen genome, so reads from elsewhere land elsewhere
#   5 pileup   pb-CpG-tools, count mode, reference CpG sites: % of reads methylated at each CpG (coverage >= 3)
# Output: results/comparative/holocentricity/tpx2/methyl/{species}/  regions.tsv, repeats.gff, bams.tsv, aligned.bam,
#         {species}.combined.bed.gz, {species}.hap1/.hap2.bed.gz (per BAM)
# usage: bash tpxl6_methyl_locus.sh [REPO] [SPECIES ...]
REPO=${1:-/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd}; shift
SPECIES=${*:-Drosera_binata Drosera_paradoxa}
ENV=${METHYL_ENV:-/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/methyl}
B=$ENV/bin; SAM=$B/samtools
FL=5000; K=31; T=${THREADS:-16}
MP=$REPO/results/comparative/holocentricity/tpx2/miniprot
OUT=$REPO/results/comparative/holocentricity/tpx2/methyl
CFG=$REPO/config/hifi_5mc.csv
Q6=regia_g4881.t1; Q5=regia_g4900.t1
die() { echo "!! $*" >&2; exit 1; }
for t in samtools minimap2 aligned_bam_to_cpg_scores; do [ -x "$B/$t" ] || die "$B/$t missing: create the env from envs/methyl.yaml"; done
[ -s "$CFG" ] || die "$CFG missing"
[ -s "$MP/tpxl6_miniprot_loci.tsv" ] || die "$MP/tpxl6_miniprot_loci.tsv missing: run tpxl6_miniprot.sh first"
cd "$REPO" || exit 1

# all K-mers of the sequences on stdin (one per line), both strands; k-mers with N (masked) or one base > 70% dropped
kmers() {
  awk -v k="$K" '{ s = toupper($0); for (i = 1; i + k - 1 <= length(s); i++) { m = substr(s, i, k); if (m !~ /N/) print m } }' \
    | awk '{ print; r = ""; for (i = length($0); i > 0; i--) { b = substr($0, i, 1)
             r = r (b == "A" ? "T" : b == "C" ? "G" : b == "G" ? "C" : b == "T" ? "A" : "N") }; print r }' \
    | awk -v k="$K" '$0 !~ /N/ { m = 0; for (b in c) delete c[b]; for (i = 1; i <= k; i++) c[substr($0, i, 1)]++
                                 for (b in c) if (c[b] > m) m = c[b]; if (m <= 0.7 * k) print }' | LC_ALL=C sort -u
}

for sp in $SPECIES; do
  d=$OUT/$sp; mkdir -p "$d"
  G=$(awk -F'\t' -v sp="$sp" '$1 == sp { print $2 }' "$MP/genomes.tsv")
  [ -s "$G" ] || { echo "== $sp: genome not found in $MP/genomes.tsv, skipped"; continue; }
  [ -s "$G.fai" ] || $SAM faidx "$G"
  echo "== $sp   $(date '+%F %T')"
  # 1 regions: name, contig, start, end, strand (the TPXL6 copies, then the best TPXL1/5-group locus)
  awk -F'\t' -v sp="$sp" -v q6="$Q6" -v q5="$Q5" -v fl="$FL" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    $c["species"] == sp && $c["query"] == q6 && $c["assigned"] == "TPXL6" { print "TPXL6\t" $c["contig"] "\t" $c["start"] "\t" $c["end"] "\t" $c["strand"] }
    $c["species"] == sp && $c["query"] == q5 && $c["score"] > best { best = $c["score"]; ctl = "control\t" $c["contig"] "\t" $c["start"] "\t" $c["end"] "\t" $c["strand"] }
    END { if (ctl != "") print ctl }' "$MP/tpxl6_miniprot_loci.tsv" > "$d/regions.tsv"
  grep -q '^TPXL6' "$d/regions.tsv" || { echo "   no TPXL6 copy, skipped"; continue; }
  while IFS=$'\t' read -r set ctg s e str; do echo "   $set $ctg:$s-$e ($str)"; done < "$d/regions.tsv"
  # repeats (all CARP Level-1 features) within FLANK of each region, for masking and for the plot
  CARP=results/$sp/repeats/carp/Repeat_Annotation_Unified.gff3
  if [ -s "$CARP" ]; then
    awk -F'\t' -v fl="$FL" 'NR == FNR { n++; C[n] = $2; L[n] = $3 - fl; H[n] = $4 + fl; want[$2] = 1; next }
      /^#/ || !($1 in want) || $9 ~ /(^|;)Parent=/ { next }
      { for (i = 1; i <= n; i++) if ($1 == C[i] && $4 <= H[i] && $5 >= L[i]) { print; break } }' "$d/regions.tsv" "$CARP" > "$d/repeats.gff"
    echo "   $(wc -l < "$d/repeats.gff") CARP repeats masked"
  else
    : > "$d/repeats.gff"; echo "   !! no CARP annotation: nothing masked, so TE k-mers will pull in extra reads"
  fi
  # 2 k-mers of each region with its repeats masked to N
  while IFS=$'\t' read -r set ctg s e str; do
    lo=$((s - FL)); [ "$lo" -ge 1 ] || lo=1; hi=$((e + FL))
    $SAM faidx "$G" "$ctg:$lo-$hi" | grep -v '^>' | tr -d '\n' | awk -v c="$ctg" -v lo="$lo" -v rf="$d/repeats.gff" '
      FILENAME == rf { if ($1 == c) { n++; A[n] = $4; Z[n] = $5 }; next }
      { s = $0; for (i = 1; i <= n; i++) { a = A[i] - lo + 1; b = Z[i] - lo + 1; if (a < 1) a = 1; if (b > length(s)) b = length(s)
          if (a <= b) { x = sprintf("%" (b - a + 1) "s", ""); gsub(/ /, "N", x); s = substr(s, 1, a - 1) x substr(s, b + 1) } }
        print s }' "$d/repeats.gff" -
  done < "$d/regions.tsv" | kmers > "$d/regions.k"
  echo "   $(wc -l < "$d/regions.k") k-mers (both strands) from the unmasked part of the regions"
  # 3 reads carrying a region k-mer, from every BAM with 5mC tags; FASTQ with MM/ML (and HP = BAM number) in the header
  : > "$d/reads.fq"; echo -e "n\tbam\ttagged_of_first_1000\treads_kept" > "$d/bams.tsv"; n=0
  while IFS= read -r bam; do
    [ -n "$bam" ] || continue; n=$((n + 1))
    if [ ! -s "$bam" ]; then echo "   !! BAM $n not found: $bam"; echo -e "$n\t$bam\tmissing\t0" >> "$d/bams.tsv"; continue; fi
    tagged=$($SAM view "$bam" 2>/dev/null | head -n 1000 | grep -cE $'\t(MM|Mm):Z:C\\+m')
    if [ "$tagged" -lt 500 ]; then
      echo "   !! BAM $n: only $tagged of its first 1000 reads carry 5mC tags, skipped: $bam"
      echo -e "$n\t$bam\t$tagged\t0" >> "$d/bams.tsv"; continue
    fi
    hp=$(( n <= 2 ? n : 0 )); before=$(( $(wc -l < "$d/reads.fq") / 4 ))
    $SAM view -@ 4 --keep-tag MM,ML,Mm,Ml "$bam" | LC_ALL=C grep -F -f "$d/regions.k" | awk -F'\t' -v hp="$hp" -v n="$n" '
      BEGIN { cp["A"] = "T"; cp["C"] = "G"; cp["G"] = "C"; cp["T"] = "A"; cp["N"] = "N" }
      int($2 / 256) % 2 || int($2 / 2048) % 2 { next }                  # secondary / supplementary
      { seq = $10; q = $11
        if (q == "*") { q = sprintf("%" length(seq) "s", ""); gsub(/ /, "~", q) }   # no stored qualities
        if (int($2 / 16) % 2) {                                         # aligned on the reverse strand: restore the read
          r = ""; for (i = length(seq); i > 0; i--) r = r cp[substr(seq, i, 1)]; seq = r
          r = ""; for (i = length(q); i > 0; i--) r = r substr(q, i, 1); q = r }
        tags = ""
        for (i = 12; i <= NF; i++) { t = $i; sub(/^Mm:/, "MM:", t); sub(/^Ml:/, "ML:", t); if (t ~ /^(MM|ML):/) tags = tags "\t" t }
        print "@" $1 "_b" n tags "\tHP:i:" hp; print seq; print "+"; print q }' >> "$d/reads.fq"
    kept=$(( $(wc -l < "$d/reads.fq") / 4 - before ))
    echo "   BAM $n: $tagged/1000 reads tagged; $kept reads carry a region k-mer   ($(basename "$bam"))"
    echo -e "$n\t$bam\t$tagged\t$kept" >> "$d/bams.tsv"
  done < <(awk -F, -v sp="$sp" 'NR > 1 && $1 == sp { print $2 }' "$CFG" | tr ';' '\n')
  [ -s "$d/reads.fq" ] || { echo "   !! no reads with 5mC tags: nothing to align"; continue; }
  # 4 align to the whole genome, keeping the tags
  "$B/minimap2" -ax map-hifi -y -Y --secondary=no -t "$T" "$G" "$d/reads.fq" 2> "$d/minimap2.log" \
    | $SAM sort -@ 4 -o "$d/aligned.bam" - 2>> "$d/minimap2.log" && $SAM index "$d/aligned.bam" || { echo "   !! alignment failed, see $d/minimap2.log"; continue; }
  echo "   aligned: $($SAM view -c -F 0x904 "$d/aligned.bam") primary alignments, $($SAM view -c -F 0x904 -q 1 "$d/aligned.bam") with MAPQ >= 1"
  # 5 per-CpG methylation (combined, and per BAM through the HP label)
  "$B/aligned_bam_to_cpg_scores" --bam "$d/aligned.bam" --output-prefix "$d/$sp" --pileup-mode count \
    --modsites-mode reference --ref "$G" --min-coverage 3 --min-mapq 1 --threads "$T" > "$d/pb_cpg.log" 2>&1 \
    || { echo "   !! pb-CpG-tools failed, see $d/pb_cpg.log"; continue; }
  rm -f "$d/reads.fq"
  echo "   pileup: $(ls "$d/$sp".*.bed.gz | xargs -n1 basename | tr '\n' ' ')"
done
echo "== done $(date '+%F %T')"
