# Tag BRAKER GFF3 genes by evidence support instead of dropping unsupported ones.
#
# Pass 1  (gene_support.tsv): a TRANSCRIPT is supported if
#            multi-exon  & >=1 intron supported, OR
#            single-exon & its exon supported.
#          A GENE is supported if ANY of its transcripts is.
# Pass 2  (BRAKER GFF3): print EVERY record, appending a gene-level suffix
#            _rnasupp  supported
#            _norna    not supported
#          to the gene id everywhere it appears: ID=, Parent=, gene_id=,
#          transcript_id= (transcript ids are <gene>.tN, so the same rewrite
#          covers them). AGAT-generated ids (agat-cds-*, agat-exon-*) are left
#          alone; they are Parent-linked, not gene-derived.
#
# mode=any     use RNA-or-protein evidence   (cols 10 introns_sup_any, 17 exons_sup_any)
# mode=rnaseq  use RNA-Seq evidence only     (cols  6 introns_sup_rnaseq, 13 exons_sup_rnaseq)
#
# Was: dropped unsupported records. Changed 2026-08 so downstream gene scans can
# SEE evidence status rather than silently losing genes (see D. binata TPXL6:
# an intact ORF dropped for lack of leaf RNA looked like a real gene absence).
BEGIN { FS="\t"; OFS="\t" }

# ---- pass 1: support table ----
FNR==NR {
    if ($0 ~ /^#/) next
    if ($1=="gene_id") next
    gid=$1; ni=$5+0
    if (mode=="rnaseq") { isup=$6+0;  esup=$13+0 }
    else                { isup=$10+0; esup=$17+0 }
    if ( (ni>0 && isup>0) || (ni==0 && esup>0) ) supgene[gid]=1
    seengene[gid]=1
    next
}

# ---- pass 2: GFF3 ----
/^#/ { print; next }
{
    gid=""
    if (match($9, /gene_id=[^;]+/)) gid=substr($9, RSTART+8, RLENGTH-8)
    if (gid=="") { print; next }          # no gene_id -> leave untouched
    sfx = (gid in supgene) ? "_rnasupp" : "_norna"
    new = gid sfx
    # rewrite the gene id wherever it appears as a whole token or as <gid>.tN
    $9 = rewrite($9, "ID", gid, new)
    $9 = rewrite($9, "Parent", gid, new)
    $9 = rewrite($9, "gene_id", gid, new)
    $9 = rewrite($9, "transcript_id", gid, new)
    print
}

function rewrite(attrs, key, old, new,   n, i, parts, out, k, v, p) {
    n = split(attrs, parts, ";")
    out = ""
    for (i=1; i<=n; i++) {
        p = parts[i]
        if (p == "") continue
        k = p; sub(/=.*/, "", k)
        v = p; sub(/^[^=]*=/, "", v)
        if (k == key) {
            if (v == old)                       v = new
            else if (index(v, old ".") == 1)     v = new substr(v, length(old)+1)
        }
        out = (out=="" ? k "=" v : out ";" k "=" v)
    }
    return out
}
