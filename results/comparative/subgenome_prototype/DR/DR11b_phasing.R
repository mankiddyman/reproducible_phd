#!/usr/bin/env Rscript
# ============================================================================
# DR11b — CAN A1/A2 BE PHASED?  (segments, not chromosomes)
#
# DR11 section 3 used CHROMOSOMES as the unit. That is wrong and DR01 already
# established why: Drosera chromosomes are fusion mosaics, which is the reason
# DR02 defined SEGMENTS (species x ancestral region x chromosome) at all.
#
# THE RIGHT TEST
#   Within one ancestral region a species has exactly TWO A segments. For two
#   species, match their A copies locus by locus by minimum total distance.
#   Then ask: does the matching pick the SAME orientation at every locus?
#     consistently one orientation -> the two A segments correspond across
#       species, so A1 and A2 are real inherited lineages and PHASING WORKS
#     50/50                        -> the pairings are noise, phasing fails
#   Exactly two orientations, so the null is genuinely 0.5.
#
#   DR07 already showed the per-locus pairings exist for the core four
#   (0.000-0.099 "independent", p to 6e-22) and NOT for regia (0.385, at the
#   null). What was never tested is whether they are spatially coherent, which
#   is what phasing requires -- subgenomes are inherited as tracts.
#
# IN   DR/out/DR02_gene_delta.csv, DR/out/DR02_segments.csv,
#      DR/out/pairwise_ks.csv
# OUT  DR/out/DR11b_phasing.csv
# FIG  DR/fig/DR11b_phasing.pdf
# ============================================================================
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})
setwd(Sys.getenv("SUBG_BASE", getwd()))
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",70), t, strrep("=",70)))

K <- read_csv("DR/out/pairwise_ks.csv", show_col_types=FALSE) %>%
     filter(!is.na(dS), dS>=0, dS<5, codons>=100)
DD <- bind_rows(K %>% transmute(locus=anchor, a=seq1, b=seq2, d=dS),
                K %>% transmute(locus=anchor, a=seq2, b=seq1, d=dS))
dm <- setNames(DD$d, paste(DD$locus, DD$a, DD$b))
gd <- function(l,x,y) ifelse(x==y, 0, unname(dm[paste(l,x,y)]))

hr("1. ASSIGN EVERY A COPY TO ITS SEGMENT")
G  <- read_csv("DR/out/DR02_gene_delta.csv", show_col_types=FALSE)
SS <- read_csv("DR/out/DR02_segments.csv", show_col_types=FALSE)
lab <- SS %>% filter(label=="A") %>%
       select(genome, region, chr, mb_lo, mb_hi, segment)
GL <- G %>% mutate(mb=mid/1e6) %>%
  inner_join(lab, by=c("genome","region","chr"), relationship="many-to-many") %>%
  filter(mb>=mb_lo, mb<=mb_hi) %>% group_by(tip) %>% slice(1) %>% ungroup() %>%
  mutate(sp=sub("Drosera_","",genome))
cat(sprintf("  A copies with a segment: %d\n", nrow(GL)))
cat("\n  A segments per (species, region) -- should be 2 for a hexaploid:\n")
print(as.data.frame(GL %>% distinct(sp, region, segment) %>%
  count(sp, region) %>% count(sp, n) %>%
  pivot_wider(names_from=n, values_from=nn, values_fill=0,
              names_prefix="segments=")), row.names=FALSE)

hr("2. ORIENTATION CONSISTENCY WITHIN EACH REGION")
cat("  Within a region, each species has 2 A segments. Match copies by minimum\n")
cat("  total distance and record which orientation wins. Null = 0.5.\n\n")
SP <- sort(unique(GL$sp))
two <- GL %>% count(sp, locus, name="k") %>% filter(k==2)
res <- bind_rows(lapply(seq_along(SP), function(i)
  bind_rows(lapply(SP[-(1:i)], function(s2) {
    s1 <- SP[i]
    L <- intersect(two$locus[two$sp==s1], two$locus[two$sp==s2])
    if (length(L) < 20) return(NULL)
    r <- bind_rows(lapply(L, function(l) {
      x <- GL %>% filter(locus==l, sp==s1) %>% arrange(segment)
      y <- GL %>% filter(locus==l, sp==s2) %>% arrange(segment)
      if (nrow(x)!=2 || nrow(y)!=2) return(NULL)
      if (x$region[1] != x$region[2] || y$region[1] != y$region[2]) return(NULL)
      d11 <- gd(l,x$tip[1],y$tip[1]); d22 <- gd(l,x$tip[2],y$tip[2])
      d12 <- gd(l,x$tip[1],y$tip[2]); d21 <- gd(l,x$tip[2],y$tip[1])
      if (any(is.na(c(d11,d22,d12,d21)))) return(NULL)
      tibble(region=x$region[1],
             segpair = paste(x$segment[1], x$segment[2], "|",
                             y$segment[1], y$segment[2]),
             orient  = if (d11+d22 < d12+d21) "straight" else "crossed")
    }))
    if (!nrow(r)) return(NULL)
    r %>% count(segpair, region, orient) %>%
      pivot_wider(names_from=orient, values_from=n, values_fill=0) %>%
      mutate(straight = if ("straight" %in% names(.)) straight else 0L,
             crossed  = if ("crossed"  %in% names(.)) crossed  else 0L,
             n = straight + crossed) %>%
      filter(n >= 10) %>%
      mutate(pair=paste(s1,s2,sep="-"),
             consistency = pmax(straight, crossed)/n,
             p = mapply(function(a,b) binom.test(max(a,b), a+b, 0.5)$p.value,
                        straight, crossed))
  }))))
if (!nrow(res)) { cat("  no segment pairs with >=10 loci\n"); quit(status=0) }
print(as.data.frame(res %>% transmute(pair, region, n,
  straight, crossed, consistency=round(consistency,3), p=signif(p,3)) %>%
  arrange(desc(consistency))), row.names=FALSE)

hr("3. SUMMARY BY SPECIES PAIR")
sm <- res %>% group_by(pair) %>%
  summarise(segment_pairs=n(), loci=sum(n),
            median_consistency=round(median(consistency),3),
            n_signif=sum(p < 0.05), .groups="drop") %>%
  arrange(desc(median_consistency))
print(as.data.frame(sm), row.names=FALSE)
cat("\n  0.5 = the matching is random, phasing IMPOSSIBLE.\n")
cat("  >0.8 = the two A segments correspond across species, phasing WORKS.\n")
write_csv(res, "DR/out/DR11b_phasing.csv")

p <- ggplot(res, aes(reorder(pair, consistency), consistency)) +
  geom_hline(yintercept=0.5, linetype="dashed", colour="grey40") +
  geom_boxplot(fill="#378ADD", alpha=0.5, outlier.size=0.8) +
  geom_jitter(width=0.12, size=1.4, alpha=0.7, colour="grey25") +
  annotate("text", x=0.6, y=0.52, hjust=0, size=3, colour="grey35",
           label="0.5 = random") +
  coord_flip(ylim=c(0.4,1)) +
  labs(title="DR11b - can A1 and A2 be phased across species?",
       subtitle=paste0("Within an ancestral region each species has two A ",
         "segments. Consistency = fraction of loci\nwhose minimum-distance ",
         "matching picks the same orientation. Two orientations, so null = 0.5."),
       x=NULL, y="orientation consistency") +
  theme_minimal(11) + theme(plot.subtitle=element_text(size=8.5, colour="grey30"))
suppressWarnings(ggsave("DR/fig/DR11b_phasing.pdf", p, width=9, height=5))

hr("VERDICT")
cat("  Core-core pairs high (>0.8)  -> A1/A2 are real inherited lineages and\n")
cat("    phasing is possible for those species.\n")
cat("  Core-core pairs near 0.5     -> DR07's pairings are not spatially\n")
cat("    coherent and phasing is impossible despite p = 6e-22.\n")
cat("  regia pairs should be LOWER than core-core either way, since DR07 found\n")
cat("  no shared duplication to align regia against.\n")
