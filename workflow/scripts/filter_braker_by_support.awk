BEGIN { FS="\t" }
FNR==NR {
    if ($0 ~ /^#/) next
    if (FNR==1 && $1=="gene_id") next
    tid=$2; ni=$5+0
    if (mode=="rnaseq") { isup=$6+0; esup=$13+0 } else { isup=$10+0; esup=$17+0 }
    if ( (ni>0 && isup>0) || (ni==0 && esup>0) ) { keep[tid]=1; keptgene[$1]=1 }
    next
}
/^#/ { print; next }
{
    tid=""
    if (match($9, /transcript_id=[^;]+/)) tid=substr($9, RSTART+14, RLENGTH-14)
    if ($3=="gene") {
        gid=""
        if (match($9, /ID=[^;]+/)) gid=substr($9, RSTART+3, RLENGTH-3)
        if (gid in keptgene) print
        next
    }
    if (tid != "" && (tid in keep)) print
}
