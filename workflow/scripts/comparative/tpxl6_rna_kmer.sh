#!/usr/bin/env bash
# tpxl6_rna_kmer.sh -- is TPXL6 transcribed in leaf RNA? Reads searched directly for exact 31-mers of the gene
#
# Why this way: no RNA alignment exists to count from (rule star_align was never run, and BRAKER4 cleans up its
# HISAT2 alignments), so the reads themselves are searched. For every species with RNA in config/annotation.csv
# (has_rna = yes, rna_dir set) whose genome has a TPXL6 copy:
#   TPXL6     every TPXL6 copy miniprot found in that genome (tpx2/miniprot/tpxl6_miniprot_loci.tsv, assigned TPXL6)
#   control   the best-scoring TPXL1/TPXL5-group locus in the same genome (TPXL6's closest relatives), so each
#             species' library is checked against a related gene, whatever its depth
# k-mers (31 bp, both strands; low-complexity ones, with one base > 70%, removed):
#   exon      every 31-mer inside the locus' coding exons (from the miniprot CDS coordinates). Reads carrying these
#             can come from mRNA, pre-mRNA or genomic DNA in the library.
#   junction  31-mers spanning an exon-exon junction with >= 8 bp on each side, minus any that also occur in the
#             locus' genomic sequence. Only spliced mRNA carries these: this is the evidence of expression.
#   k-mers shared between TPXL6 and the control are dropped from both.
# reads: every *_R1_* and *_R2_* fastq.gz in rna_dir, files scanned in parallel; a read counts once per k-mer set.
# Output: results/comparative/holocentricity/tpx2/rna_kmer/{species}/... and summary.tsv
# usage: bash tpxl6_rna_kmer.sh [REPO] [PARALLEL_FILES]
REPO=${1:-/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd}
PAR=${2:-16}
ENV=${PHYLO_ENV:-/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/phylo}
SAM=$ENV/bin/samtools
MP=$REPO/results/comparative/holocentricity/tpx2/miniprot
OUT=$REPO/results/comparative/holocentricity/tpx2/rna_kmer
Q6=regia_g4881.t1; Q5=regia_g4900.t1; K=31
die() { echo "!! $*" >&2; exit 1; }
[ -x "$SAM" ] || die "$SAM missing"
[ -s "$MP/tpxl6_miniprot_loci.tsv" ] || die "run tpxl6_miniprot.sh and tpxl6_miniprot_summary.R first"
mkdir -p "$OUT" || die "cannot create $OUT"
echo -e "species\tset\texon_kmers\tjunction_kmers\treads_scanned\texon_reads\tjunction_reads\texon_per_million\tjunction_per_million" > "$OUT/summary.tsv"

# exon coordinates (contig, start, end) of one locus, genomic order, from the miniprot GFF
exons_of() {  # gff contig start end query
  awk -F'\t' -v c="$2" -v s="$3" -v e="$4" -v q="$5" '
    $3 == "mRNA" && $1 == c && $4 == s && $5 == e && index($9, "Target=" q " ") { match($9, /ID=[^;]+/); id = substr($9, RSTART + 3, RLENGTH - 3) }
    $3 == "CDS" && id != "" && index($9, "Parent=" id ";") { print $1 "\t" $4 "\t" $5 }' "$1" | sort -k2,2n
}
# all K-mers of the sequences on stdin (one per line), both strands, low-complexity removed
kmers() {
  awk -v k="$K" '{ s = toupper($0); for (i = 1; i + k - 1 <= length(s); i++) print substr(s, i, k) }' \
    | awk '{ print; r = ""; for (i = length($0); i > 0; i--) { b = substr($0, i, 1)
             r = r (b == "A" ? "T" : b == "C" ? "G" : b == "G" ? "C" : b == "T" ? "A" : "N") }; print r }' \
    | awk -v k="$K" '$0 !~ /N/ { m = 0; for (b in c) delete c[b]; for (i = 1; i <= k; i++) c[substr($0, i, 1)]++
                                 for (b in c) if (c[b] > m) m = c[b]; if (m <= 0.7 * k) print }' | sort -u
}
# K-mers spanning the middle of 60-bp junction sequences (30 + 30), >= 8 bp on each side, both strands
jkmers() {
  awk -v k="$K" '{ s = toupper($0); for (i = 8; i <= 23; i++) print substr(s, i, k) }' \
    | awk '{ print; r = ""; for (i = length($0); i > 0; i--) { b = substr($0, i, 1)
             r = r (b == "A" ? "T" : b == "C" ? "G" : b == "G" ? "C" : b == "T" ? "A" : "N") }; print r }' \
    | awk '$0 !~ /N/' | sort -u
}
seq_of() { "$SAM" faidx "$1" "$2:$3-$4" | grep -v '^>' | tr -d '\n'; echo; }

while IFS=, read -r sp has_rna rna_dir _h _a _t stage; do
  [ "$sp" = species_id ] && continue
  [ "$has_rna" = yes ] && [ -n "$rna_dir" ] || continue
  G=$REPO/results/$sp/assembly_final/$stage/${sp}_chr.fa
  gff=$MP/$sp.gff; d=$OUT/$sp; mkdir -p "$d"
  [ -s "$gff" ] && [ -s "$G" ] || { echo "== $sp: no miniprot GFF or genome, skipped"; continue; }
  awk -F'\t' -v sp="$sp" -v q6="$Q6" -v q5="$Q5" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    $c["species"] == sp && $c["query"] == q6 && $c["assigned"] == "TPXL6" { print "TPXL6\t" $c["contig"] "\t" $c["start"] "\t" $c["end"] "\t" q6 }
    $c["species"] == sp && $c["query"] == q5 { print "control\t" $c["contig"] "\t" $c["start"] "\t" $c["end"] "\t" q5 "\t" $c["score"] }' \
    "$MP/tpxl6_miniprot_loci.tsv" > "$d/loci.all"
  grep -q '^TPXL6' "$d/loci.all" || { echo "== $sp: no TPXL6 copy in the genome, skipped"; continue; }
  { grep '^TPXL6' "$d/loci.all"; grep '^control' "$d/loci.all" | sort -t$'\t' -k6,6gr | head -1 | cut -f1-5; } > "$d/loci.tsv"
  echo "== $sp   $(date '+%F %T')"
  for set in TPXL6 control; do : > "$d/$set.exon.seq"; : > "$d/$set.junction.seq"; done
  : > "$d/locus.seq"
  while IFS=$'\t' read -r set ctg s e q; do
    exons_of "$gff" "$ctg" "$s" "$e" "$q" > "$d/tmp.exons"
    seq_of "$G" "$ctg" "$s" "$e" >> "$d/locus.seq"
    echo "   $set $ctg:$s-$e  $(wc -l < "$d/tmp.exons") exons"
    while IFS=$'\t' read -r c a b; do seq_of "$G" "$c" "$a" "$b"; done < "$d/tmp.exons" >> "$d/$set.exon.seq"
    awk 'NR > 1 { print prev_end "\t" $2 } { prev_end = $3; c = $1 }' "$d/tmp.exons" | while IFS=$'\t' read -r pe ns; do
      ctg1=$(head -1 "$d/tmp.exons" | cut -f1)
      echo "$(seq_of "$G" "$ctg1" $((pe - 29)) "$pe")$(seq_of "$G" "$ctg1" "$ns" $((ns + 29)))"
    done >> "$d/$set.junction.seq"
  done < "$d/loci.tsv"
  kmers < "$d/locus.seq" > "$d/locus.k"                      # every 31-mer of the genomic loci, introns included
  for set in TPXL6 control; do
    kmers < "$d/$set.exon.seq" > "$d/$set.exon.k"
    jkmers < "$d/$set.junction.seq" | comm -23 - "$d/locus.k" > "$d/$set.junction.k"
  done
  comm -12 "$d/TPXL6.exon.k" "$d/control.exon.k" > "$d/shared.k"
  for set in TPXL6 control; do
    comm -23 "$d/$set.exon.k" "$d/shared.k" > "$d/$set.exon.k.u" && mv "$d/$set.exon.k.u" "$d/$set.exon.k"
  done
  cat "$d"/TPXL6.exon.k "$d"/TPXL6.junction.k "$d"/control.exon.k "$d"/control.junction.k | sort -u > "$d/all.k"
  echo "   k-mers: TPXL6 exon $(wc -l < "$d/TPXL6.exon.k"), junction $(wc -l < "$d/TPXL6.junction.k") | control exon $(wc -l < "$d/control.exon.k"), junction $(wc -l < "$d/control.junction.k")"
  ls "$rna_dir"/*_R1_*.fastq.gz "$rna_dir"/*_R2_*.fastq.gz 2>/dev/null > "$d/fastq.list"
  [ -s "$d/fastq.list" ] || { echo "   !! no fastq.gz in $rna_dir"; continue; }
  echo "   scanning $(wc -l < "$d/fastq.list") fastq files ($(du -chL $(cat "$d/fastq.list") | tail -1 | cut -f1)) with $PAR in parallel"
  rm -rf "$d/scan"; mkdir -p "$d/scan"
  export d
  xargs -P "$PAR" -I{} bash -c 'b=$(basename "{}" .fastq.gz)
    zcat "{}" | awk '"'"'NR % 4 == 2 { n++; print } END { print n > "/dev/stderr" }'"'"' 2> "$d/scan/$b.n" \
      | LC_ALL=C grep -F -f "$d/all.k" > "$d/scan/$b.hits"' < "$d/fastq.list"
  N=$(cat "$d"/scan/*.n | awk '{ s += $1 } END { print s + 0 }')
  { echo -e "file\treads\tTPXL6_exon_reads\tTPXL6_junction_reads\tcontrol_junction_reads"
    for h in "$d"/scan/*.hits; do b=$(basename "$h" .hits)
      t6e=$(LC_ALL=C grep -c -F -f "$d/TPXL6.exon.k" "$h"); t6j=0; cj=0
      [ -s "$d/TPXL6.junction.k" ] && t6j=$(LC_ALL=C grep -c -F -f "$d/TPXL6.junction.k" "$h")
      [ -s "$d/control.junction.k" ] && cj=$(LC_ALL=C grep -c -F -f "$d/control.junction.k" "$h")
      echo -e "$b\t$(cat "$d/scan/$b.n")\t$t6e\t$t6j\t$cj"; done; } > "$d/per_file.tsv"
  cat "$d"/scan/*.hits > "$d/hits.txt"
  for set in TPXL6 control; do
    n=$(LC_ALL=C grep -c -F -f "$d/$set.exon.k" "$d/hits.txt"); j=0
    [ -s "$d/$set.junction.k" ] && j=$(LC_ALL=C grep -c -F -f "$d/$set.junction.k" "$d/hits.txt")
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n" "$sp" "$set" "$(wc -l < "$d/$set.exon.k")" "$(wc -l < "$d/$set.junction.k")" "$N" "$n" "$j" \
      "$(echo "$n $N" | awk '{ print ($2 > 0) ? 1e6 * $1 / $2 : 0 }')" "$(echo "$j $N" | awk '{ print ($2 > 0) ? 1e6 * $1 / $2 : 0 }')" >> "$OUT/summary.tsv"
    echo "   $set: $n of $N reads carry an exon k-mer; $j carry a splice-junction k-mer (= spliced mRNA)"
  done
done < "$REPO/config/annotation.csv"
echo; column -t -s$'\t' "$OUT/summary.tsv" 2>/dev/null || cat "$OUT/summary.tsv"
echo "== done $(date '+%F %T')"
