suppressMessages(library(data.table))
setwd("/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/comparative/holocentricity/tpx2")

d <- fread("tpxl6_map.tsv", col.names=c("sp","ctg","qstart","qend","sstart","send",
                                        "frame","evalue","pident","bits","inside"))
d[, gene := gsub("^ID=|;.*$", "", inside)]
d[, gmin := pmin(sstart, send)]
setorder(d, sp, ctg, gmin)
d[, newloc := as.integer(ctg != shift(ctg, fill="") | (gmin - shift(gmin, fill=-1e9)) > 50000), by=sp]
d[, locus := cumsum(newloc), by=sp]

locsum <- d[, .(tot=round(sum(bits)), nhsp=.N, nframe=uniqueN(frame),
                qmin=min(qstart), qmax=max(qend), best_e=min(evalue),
                ctg=ctg[1], genes=paste(unique(gene[gene!="-"]), collapse=",")), by=.(sp,locus)]
setorder(locsum, sp, -tot)
cat("=== ALL loci per species (best first) ===\n"); print(locsum, nrows=60)

best <- locsum[, .SD[1], by=sp]
cat("\n=== BEST locus per species ===\n"); print(best)

hs <- merge(d, best[, .(sp, locus)], by=c("sp","locus"))
setorder(hs, sp, qstart)
cat("\n=== HSPs of the best locus (query coords = which part of the protein) ===\n")
print(hs[, .(sp, qstart, qend, frame, pident=round(pident), evalue, gene)], nrows=60)
