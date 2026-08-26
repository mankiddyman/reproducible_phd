#!/usr/bin/env Rscript
.p <- grep("micromamba_envs/smk", .libPaths(), value=TRUE); if (length(.p)) .libPaths(.p)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
BASE <- Sys.getenv("SUBG_BASE", getwd()); setwd(BASE); NPERM <- 9999
hr <- function(t) cat(sprintf("\n%s\n%s\n%s\n", strrep("=",66), t, strrep("=",66)))
seg <- read_csv("trees/qc/segment_delta_calibrated.csv", show_col_types=FALSE)

hr("1. per-region composition vs a binomial null with the SAME marginal skew")
cat("  Constitution -> exactly 2A:1B wherever 3 segments survive (under-dispersed).\n")
cat("  Marginal skew only -> Binomial(n, p_species).\n\n")
set.seed(1)
out <- bind_rows(lapply(c(0.05, 0.08, 0.12), function(B) {
  cc <- seg %>% mutate(side = case_when(med_delta < -B ~ "A",
                                        med_delta >  B ~ "B", TRUE ~ "amb")) %>%
        filter(side != "amb")
  bind_rows(lapply(split(cc, cc$genome), function(d) {
    rg <- d %>% group_by(region) %>%
          summarise(n=n(), a=sum(side=="A"), .groups="drop") %>% filter(n>=2)
    if (nrow(rg) < 4) return(NULL)
    p <- sum(rg$a)/sum(rg$n)
    obs_v <- var(rg$a/rg$n)
    nul <- vapply(seq_len(NPERM), function(i)
      var(rbinom(nrow(rg), rg$n, p)/rg$n), numeric(1))
    three <- rg %>% filter(n==3)
    tibble(band=B, genome=d$genome[1], regions=nrow(rg), p_A=round(p,3),
           obs_var=round(obs_v,4), null_var=round(mean(nul),4),
           disp=round(obs_v/mean(nul),2),
           p_under=(1+sum(nul <= obs_v))/(NPERM+1),
           n3=nrow(three), n3_2to1=sum(three$a==2),
           exp3=round(3*p^2*(1-p),3))
  }))
}))
print(as.data.frame(out %>% arrange(band, genome)), row.names=FALSE)
cat("\n  disp < 1 with small p_under = under-dispersed = fixed constitution.\n")
cat("  disp ~ 1 = the 2:1 ratio is nothing but the marginal skew.\n")

hr("2. three-segment regions pooled across species")
for (B in c(0.05, 0.08, 0.12)) {
  cc <- seg %>% mutate(side=case_when(med_delta < -B ~ "A", med_delta > B ~ "B",
                                      TRUE ~ "amb")) %>% filter(side != "amb")
  rg <- cc %>% group_by(genome, region) %>%
        summarise(n=n(), a=sum(side=="A"), .groups="drop") %>% filter(n==3)
  p <- sum(cc$side=="A")/nrow(cc)
  k <- sum(rg$a==2)
  cat(sprintf("  band %.2f : %d of %d three-segment regions are 2A:1B",
              B, k, nrow(rg)))
  cat(sprintf("   (expected %.2f)  binom p = %.3g\n", 3*p^2*(1-p),
              binom.test(k, nrow(rg), 3*p^2*(1-p))$p.value))
  cat(sprintf("             A-count distribution: %s\n",
              paste(names(table(rg$a)), table(rg$a), sep="x", collapse=" ")))
}

hr("3. capensis — is it a doubled AAB?")
cat("  Doubled AAB = AAAABB = 4:2 = ratio 2.0. Observed 3.9-4.8.\n\n")
for (B in c(0.05, 0.08, 0.12)) {
  cp <- seg %>% filter(genome=="Drosera_capensis") %>%
    mutate(side=case_when(med_delta < -B ~ "A", med_delta > B ~ "B", TRUE ~ "amb")) %>%
    group_by(region) %>%
    summarise(A=sum(side=="A"), B=sum(side=="B"), amb=sum(side=="amb"), .groups="drop")
  cat(sprintf("  band %.2f\n", B)); print(as.data.frame(cp), row.names=FALSE); cat("\n")
}
