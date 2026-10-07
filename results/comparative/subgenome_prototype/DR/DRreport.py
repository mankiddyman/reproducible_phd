#!/usr/bin/env python3
"""DRreport.py -- subgenome phasing and dating of Drosera: what was done, how good it is, what is open.

One self-contained HTML page (figures embedded) plus numbers.txt, built from the outputs of the DR pipeline in
results/comparative/subgenome_prototype. Every number, table and figure is computed here from those outputs, and
every threshold quoted in the text is read from the R script that applies it. A section that cannot be built is
listed in red at the top with its error; missing inputs are listed at the end.

  0 summary  1 material and terms  2 pipeline at a glance  3 how the labels are made  4 do the labels hold up
  5 constitution, trees, riparian  6 dating  7 limitations  8 next steps and asks  appendix

Computed here rather than by an existing pipeline script (definitions in the code below):
  - the Nepenthes control and the per-gene noise of Δ, re-derived from pairwise_ks.csv exactly as DR02_label.R
    does and checked against the offset DR02 subtracted
  - autocorrelation of Δ along chromosomes against within-chromosome shuffles
  - the cross-species dS test, per subgenome, with a permuted-label control
  - MCMCtree convergence (ESS as coda::effectiveSize; pooled = sum over replicate chains, as coda does for an
    mcmc.list), posterior means and 95% HPD intervals, rooting and root-calibration checks

Usage (from results/comparative/subgenome_prototype):
  python DR/DRreport.py [--out DR/report] [--genespace ../genespace/results] [--rscript PATH] [--no_tests]
Writes OUT/report.html, OUT/numbers.txt and OUT/fig/*.png.
"""
import argparse
import glob
import hashlib
import html
import math
import os
import re
import shutil
import subprocess

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from report_kit import Report, Html, fi, ff, ok, binom_p, norm_sf, wilson, need, code_const

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="DR/report")
ap.add_argument("--genespace", default="../genespace/results")
ap.add_argument("--rscript", default="/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/Rscript")
ap.add_argument("--no_tests", action="store_true", help="do not run DR/DRtest_phasing.R")
ap.add_argument("--perm", type=int, default=200, help="shuffles behind each chance envelope")
A = ap.parse_args()

R = Report("Subgenome phasing and dating of <i>Drosera</i>: what was done, how good it is, what is open",
           A.out, "DR/DRreport.py")
rng = np.random.default_rng(1)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

# ---- constants of this report (everything else is read from the pipeline's code below)
MIN_UNIT = 8          # measured genes a (species, region, chromosome) unit needs to enter the noise estimate
MAX_LAG = 60          # lags (in measured genes) for the autocorrelation of Δ
N_CTRL = 20           # permuted-label repeats behind the cross-species control
ESS_MIN = 200         # effective sample size below which posterior ages are not cited
TIME_UNIT_MA = 100.0  # MCMCtree time unit: calibrations such as B(0.672,1.003) are in units of 100 Myr

DR02, DR02C, DR14, DR21, DR22 = ("DR/DR02_label.R", "DR/DR02c_blocks.R", "DR/DR14_phasing_quality.R",
                                 "DR/DR21_copynumber.R", "DR/DR22_truncation.R")
RIP, DR05F, DR10 = "DR/DRrip6_riparian_FINAL.R", "DR/DR05f_control.R", "DR/DR10_figure.R"
P = dict(meta="DR/locus_meta.tsv", locsum="DR/out/DR00_locus_summary.csv", ks="DR/out/pairwise_ks.csv",
         gd="DR/out/DR02_gene_delta.csv", seg="DR/out/DR02_segments.csv", agree="DR/out/DR02_agreement.csv",
         pb="DR/out/DR02_propagated_blocks.csv", pur="DR/out/DR14_segment_quality.csv",
         conf="DR/out/DR14_gene_confidence.csv", cn="DR/out/DR21_copynumber.csv", trunc="DR/out/DR22_truncation.csv",
         c05f="DR/out/DR05f_control.csv", q06="DR/out/DR06_quartet_counts.csv", tracts="DR/out/DR25_tracts.csv",
         dropped="DR/out/DR25_dropped.csv", frac="fractionation_by_chrpair.csv",
         bed=os.path.join(A.genespace, "combBed.txt"))
PIPE_FIGS = dict(astralpro="DR/fig/DR10_MAIN_astralpro.png", concat13="DR/fig/DR05_13tip_tree.png",
                 riparian="DR/fig/DR25_riparian_AB.png", dr24="DR/fig/DR24_all13_dated.pdf",
                 modelfit="DR/fig/DR10_SUPP_modelfit.pdf")


# ============================================================================ code facts
def code_text(path):
    return open(path).read() if R.have(path) else ""


def code_str(path, pattern, default=None, flags=0):
    m = re.search(pattern, code_text(path), flags)
    return m.group(1) if m else default


_dros = code_str(DR02, r"DROS <- c\(([^)]*)\)", "", re.S)
DROS = re.findall(r'"([^"]+)"', _dros) or ["Drosera_regia", "Drosera_binata", "Drosera_paradoxa",
                                             "Drosera_scorpioides", "Drosera_capensis"]
NEP, DIO = "Nepenthes_gracilis", "Dionaea_muscipula"
DSMAX, = code_const(DR02, r"DSMAX <- (\d+)")
DCONV, = code_const(DR02, r"DCONV <- ([\d.]+)")
NPERM, = code_const(DR02, r"NPERM <- (\d+)L")
MINSEG, = code_const(DR02, r"MINSEG <- (\d+)L")
SEG_LO, = code_const(DR02, r"med_delta < (-[\d.]+)")
SEG_HI, = code_const(DR02, r"med_delta >\s+([\d.]+)")
MINDENS, = code_const(DR02, r"MINDENS <- ([\d.]+)")
PRED_A, PRED_B = code_const(DR02, r"A ~ (-[\d.]+), B ~ \+([\d.]+)", 2)
SEP, = code_const(DR02, r"separation ([\d.]+)\.")
NOISE_PRED, = code_const(DR02, r"noise sd ~ ([\d.]+)")
PUR_MIN, = code_const(DR02C, r"BLK\$purity >= ([\d.]+)")
NVOTE_MIN, = code_const(DR02C, r"BLK\$nvote >= (\d+)")
NEAR_MB, = code_const(DR14, r"d_nearest < ([\d.]+)\)")
FAR_MB, = code_const(DR14, r"d_nearest > ([\d.]+)\)")
GENE_AMBIG, = code_const(DR14, r"gene_call = ifelse\(abs\(delta\) < ([\d.]+)")
PUR_MARG, = code_const(DR14, r"marginal=sum\(purity < ([\d.]+)\)")
BRAID_PUR, = code_const(RIP, r"/length\(z\) >= ([\d.]+)\)")
TRACT_MIN, = code_const(RIP, r"TR\$n >= (\d+)")
Q1_RED, = code_const(DR10, r"short = q1 < ([\d.]+)")
F05_CLEAR, = code_const(DR05F, r"rg > max\(ot\) \+ ([\d.]+)")
F05_NONE, = code_const(DR05F, r"rg <= max\(ot\) \+ ([\d.]+)")
COL = {"A": code_str(RIP, r'BARA <- "(#[0-9A-Fa-f]{6})"', "#1D9E75"),
       "B": code_str(RIP, r'BARB <- "(#[0-9A-Fa-f]{6})"', "#D85A30"),
       "ambiguous": "#B4B2A9", "unresolved": "#D3D1C7", "control": "#888780"}
_pred = re.search(r"PRED <- c\((.*?)\)", code_text(DR22), re.S)
PRED = {k: int(v) for k, v in re.findall(r"(Drosera_\w+)=(\d+)", _pred.group(1))} if _pred else {}
SPC = dict(zip(DROS, ["#534AB7", "#378ADD", "#D4537E", "#5F5E5A", "#0C447C"]))


def sp(g):
    return g.replace("Drosera_", "").replace("_gracilis", "").replace("_muscipula", "")


def isp(g):
    return Html("<i>%s</i>" % g.replace("_", " ").replace("Drosera ", "D. "))


def html_text(s):
    return html.escape(str(s))


# ============================================================================ inputs
D = {}


def load(key, fn):
    p = P[key]
    if R.have(p):
        v = R.guarded("read " + p, lambda: fn(p))
        if v is not None:
            D[key] = v


load("meta", lambda p: need(pd.read_csv(p, sep="\t", dtype=str), ["locus", "tip", "gene", "genome", "chr"], p))
load("locsum", lambda p: pd.read_csv(p, dtype={"locus": str}))
load("ks", lambda p: need(pd.read_csv(p, usecols=["anchor", "seq1", "seq2", "codons", "dS"],
                                      dtype={"anchor": str, "seq1": str, "seq2": str}),
                          ["anchor", "seq1", "seq2", "codons", "dS"], p))
load("gd", lambda p: need(pd.read_csv(p, dtype={"locus": str, "gene": str, "chr": str, "region": str}),
                          ["locus", "genome", "tip", "gene", "chr", "mid", "region", "delta_raw", "delta"], p))
load("seg", lambda p: need(pd.read_csv(p, dtype={"chr": str, "region": str}),
                           ["genome", "region", "chr", "segment", "n", "med_delta", "label"], p))
load("agree", lambda p: need(pd.read_csv(p, dtype={"locus": str}), ["genome", "delta", "fp"], p))
load("pb", lambda p: need(pd.read_csv(p, dtype={"gene": str, "chr": str, "region": str}),
                          ["genome", "gene", "chr", "region", "label"], p))
load("pur", lambda p: need(pd.read_csv(p), ["genome", "label", "votes", "purity"], p))
load("conf", lambda p: need(pd.read_csv(p, dtype={"gene": str}), ["genome", "d_nearest"], p))
load("cn", lambda p: need(pd.read_csv(p, dtype={"locus": str}), ["genome", "locus", "A", "B"], p))
load("trunc", lambda p: need(pd.read_csv(p, dtype={"locus": str}), ["genome", "nA", "nA_max", "unlab"], p))
load("c05f", lambda p: need(pd.read_csv(p), ["focal", "loci", "focal_Dionaea"], p))
load("q06", lambda p: need(pd.read_csv(p), ["test", "n", "top_n", "mid_n", "low_n"], p))
load("tracts", lambda p: need(pd.read_csv(p), ["genome", "chr", "region", "label"], p))
load("dropped", lambda p: pd.read_csv(p))
load("frac", lambda p: need(pd.read_csv(p), ["exp_pair", "chrA", "chrB", "retained_more"], p))


def read_bed(p):
    b = pd.read_csv(p, sep="\t", quoting=3, low_memory=False, dtype={"chr": str, "id": str},
                    usecols=lambda c: c in ("genome", "id", "chr", "start", "end", "isArrayRep"))
    need(b, ["genome", "id", "chr", "isArrayRep"], p)
    rep = b.isArrayRep.astype(str).str.upper()
    return b[rep.isin(["TRUE", "T", "NA", "NAN", ""])]     # as DR02: is.na(isRep) | isRep


load("bed", read_bed)


def have_all(*keys):
    missing = [P[k] for k in keys if k not in D]
    if missing:
        raise FileNotFoundError("needs " + ", ".join(missing))


# ============================================================================ distances
def dist_table():
    """pairwise_ks.csv after DR02's filter, both orientations, first row per (locus, a, b) as DR02's lookup."""
    k = D["ks"]
    k = k[k.dS.notna() & (k.dS >= 0) & (k.dS < DSMAX) & (k.codons >= 100)]
    kk = pd.concat([k.rename(columns={"anchor": "locus", "seq1": "a", "seq2": "b"}),
                    k.rename(columns={"anchor": "locus", "seq1": "b", "seq2": "a"})], ignore_index=True)
    kk = kk[["locus", "a", "b", "dS"]].drop_duplicates(["locus", "a", "b"], keep="first")
    # one hashed index, built once: lookups are then a reindex rather than a merge per call
    return pd.Series(kk.dS.values, index=kk.locus.astype(str) + "\t" + kk.a.astype(str) + "\t" + kk.b.astype(str))


def look(df, x, y, name, K):
    """Add column `name` = dS between tips df[x] and df[y] at df.locus (NaN when not measured)."""
    out = df.copy()
    key = out.locus.astype(str) + "\t" + out[x].astype(str) + "\t" + out[y].astype(str)
    out[name] = K.reindex(key.values).values
    return out


S = {}   # every computed result, by topic


def compute_distances():
    have_all("ks")
    S["K"] = dist_table()
    S["ks_rows"], S["ks_loci"] = len(D["ks"]), D["ks"].anchor.nunique()
    S["ks_kept"] = len(S["K"]) // 2


R.guarded("distance table", compute_distances)


# ============================================================================ material and pipeline counts
def material():
    have_all("bed", "meta")
    bed, meta = D["bed"], D["meta"]
    rows = []
    for g in [NEP, DIO] + DROS:
        b, m = bed[bed.genome == g], meta[meta.genome == g]
        role = ("outgroup; its chromosomes name the ancestral regions" if g == NEP else
                "reference: its two homeologous chromosome sets give Dio_A and Dio_B" if g == DIO else
                "phased into subgenomes A and B")
        rows.append([isp(g), role, fi(len(b)), fi(b.chr.nunique()), fi(len(m)), fi(m.locus.nunique())])
        R.put("genes." + sp(g), len(b))
    S["material"] = rows
    S["dros_genes"] = int(bed.genome.isin(DROS).sum())
    S["loci"], S["tips"] = meta.locus.nunique(), meta.tip.nunique()
    R.put("loci", S["loci"])
    R.put("drosera_genes", S["dros_genes"])


R.guarded("material", material)


# ============================================================================ Δ: Nepenthes control, noise
def axis_and_nepenthes():
    """DR02_label.R sections 1 and 4b, re-derived: the Dionaea axis loci and Δ for Nepenthes (true Δ = 0)."""
    have_all("meta", "ks", "frac")
    K, frac = S["K"], D["frac"]
    loc = D["meta"].drop_duplicates("tip", keep="first").sort_values("tip", kind="mergesort")  # group_by(tip) slice(1)
    side = {}
    for a, b, more in zip(frac.chrA, frac.chrB, frac.retained_more):
        side[a], side[b] = ("A", "B") if more == a else ("B", "A")
    dio = loc[loc.genome == DIO].copy()
    dio["side"] = dio.chr.map(side)
    dio = dio[dio.side.notna()]
    cnt = dio.groupby(["locus", "side"]).tip.agg(["size", "first"]).unstack("side")
    one = (cnt[("size", "A")] == 1) & (cnt[("size", "B")] == 1)
    ax = pd.DataFrame({"locus": cnt.index[one], "DA": cnt[("first", "A")][one].values,
                       "DB": cnt[("first", "B")][one].values})
    S["axis_one_each"] = len(ax)
    ax = look(ax, "DA", "DB", "dAB", K)
    ax = ax[ax.dAB.notna()]
    S["axis_measured"] = len(ax)
    S["axis_conv"] = int((ax.dAB < DCONV).sum())
    ax = ax[ax.dAB >= DCONV]
    S["axis"] = len(ax)
    nep = loc[loc.genome == NEP].groupby("locus").tip.first().rename("N").reset_index()
    n = nep.merge(ax, on="locus")
    n = look(look(n, "N", "DA", "dRA", K), "N", "DB", "dRB", K)
    n = n[n.dRA.notna() & n.dRB.notna()].copy()
    n["delta"] = (n.dRA - n.dRB) / n.dAB
    S["nep"] = n.delta.values
    x = S["nep"]
    S["nep_median"], S["nep_mean"], S["nep_sd"] = float(np.median(x)), float(np.mean(x)), float(np.std(x, ddof=1))
    S["nep_rsd"] = float(1.4826 * np.median(np.abs(x - np.median(x))))
    S["nep_p"] = wilcoxon_p(x)
    if "gd" in D:
        off = (D["gd"].delta_raw - D["gd"].delta)
        S["off_file"] = float(off.median())
        S["off_spread"] = float(off.max() - off.min())
    for k in ("axis_one_each", "axis_measured", "axis_conv", "axis", "nep_median", "nep_sd", "nep_rsd", "nep_p"):
        R.put("delta." + k, S[k])
    R.put("delta.nepenthes_n", len(x))


def wilcoxon_p(x):
    """Two-sided Wilcoxon signed-rank test of median 0, normal approximation with continuity and tie correction
    (R's wilcox.test with exact = FALSE)."""
    x = np.asarray(x, float)
    x = x[x != 0]
    n = len(x)
    if n < 10:
        return float("nan")
    r = pd.Series(np.abs(x)).rank().values
    w = r[x > 0].sum()
    mu = n * (n + 1) / 4
    _, t = np.unique(np.abs(x), return_counts=True)
    var = n * (n + 1) * (2 * n + 1) / 24 - (t ** 3 - t).sum() / 48
    z = w - mu
    z = (z - 0.5 * np.sign(z)) / math.sqrt(var)
    return 2 * norm_sf(abs(z))


R.guarded("Nepenthes control", axis_and_nepenthes)


def unit_noise():
    """Spread of Δ inside (species, region, chromosome) units with >= MIN_UNIT measured genes: the per-gene noise
    if each unit is one ancestry (a unit spanning a boundary would inflate it, which the size check addresses)."""
    have_all("gd")
    g = D["gd"].groupby(["genome", "region", "chr"]).delta
    u = pd.DataFrame({"n": g.size(), "sd": g.std(ddof=1)}).reset_index()
    u = u[u.n >= MIN_UNIT]
    S["units"] = u
    S["unit_sd"] = float(math.sqrt(((u.n - 1) * u.sd ** 2).sum() / (u.n - 1).sum()))
    S["unit_sd_median"] = float(u.sd.median())
    S["unit_rho"] = float(u.n.corr(u.sd, method="spearman"))
    S["unit_sd_by_sp"] = {gn: float(math.sqrt(((x.n - 1) * x.sd ** 2).sum() / (x.n - 1).sum()))
                          for gn, x in u.groupby("genome")}
    # DR02's design arithmetic (header of DR02_label.R), redone with the measured noise: how many sd apart the
    # A and B centres are for one gene
    S["sep_sd_gene"] = SEP / S["unit_sd"]
    S["sep_sd_gene_design"] = SEP / NOISE_PRED
    for k in ("unit_sd", "unit_sd_median", "unit_rho", "sep_sd_gene", "sep_sd_gene_design"):
        R.put("delta." + k, S[k])
    R.put("delta.units", len(u))
    for gn, v in S["unit_sd_by_sp"].items():
        R.put("delta.unit_sd." + sp(gn), v)


R.guarded("per-gene noise", unit_noise)


def segments():
    have_all("seg", "gd")
    seg, gd = D["seg"], D["gd"]
    tr = seg[["genome", "region", "chr"]].drop_duplicates()
    S["tracks"], S["segments"] = len(tr), len(seg)
    S["seg_counts"] = seg.groupby(["genome", "label"]).size().unstack(fill_value=0)
    called = seg[seg.label.isin(["A", "B"])]
    S["called"] = called
    if len(called):
        # DR02's own per-segment score: conf = |median Δ| / (MAD / sqrt(n)), distance of the median from 0 in
        # standard errors. A median's standard error is about 1.2533 times larger than MAD / sqrt(n).
        conf = called["conf"] if "conf" in called else (called.med_delta.abs() / (called.mad_delta / np.sqrt(called.n)))
        S["seg_n_median"] = float(called.n.median())
        S["conf3"] = float((conf >= 3).mean())
        S["conf3_strict"] = float((conf / 1.2533 >= 3).mean())
        if "sep_sd_gene" in S:
            S["sep_at_median"] = S["sep_sd_gene"] * math.sqrt(S["seg_n_median"])
            S["sep_at_median_design"] = S["sep_sd_gene_design"] * math.sqrt(S["seg_n_median"])
        for k in ("seg_n_median", "conf3", "conf3_strict", "sep_at_median", "sep_at_median_design"):
            if k in S:
                R.put("segments." + k, S[k])
    S["measured_genes"] = gd.groupby("genome").gene.nunique()
    S["measured_total"] = int(gd[["genome", "gene"]].drop_duplicates().shape[0])
    R.put("segments.tracks", S["tracks"])
    R.put("segments.total", S["segments"])
    R.put("delta.measured_genes", S["measured_total"])
    for (g, l), v in S["seg_counts"].stack().items():
        R.put("segments.%s.%s" % (sp(g), l), v)


R.guarded("segments", segments)


def coverage():
    have_all("pb", "bed", "gd")
    pb, bed = D["pb"], D["bed"]
    S["pb_dup"] = int(pb.duplicated(["genome", "gene"]).sum())
    tot = bed[bed.genome.isin(DROS)].groupby("genome").size()
    lab = pb[pb.label.isin(["A", "B"])].groupby("genome").gene.nunique()
    meas = S["measured_genes"]
    cov = pd.DataFrame({"total": tot, "labelled": lab, "measured": meas}).reindex(DROS).fillna(0)
    cov["pct_measured"] = 100 * cov.measured / cov.total
    cov["pct_labelled"] = 100 * cov.labelled / cov.total
    S["cov"] = cov
    S["labelled_total"] = int(cov.labelled.sum())
    S["pct_measured_all"] = 100 * cov.measured.sum() / cov.total.sum()
    S["pct_labelled_all"] = 100 * cov.labelled.sum() / cov.total.sum()
    for g, r in cov.iterrows():
        R.put("coverage.%s.pct_measured" % sp(g), r.pct_measured)
        R.put("coverage.%s.pct_labelled" % sp(g), r.pct_labelled)
    R.put("coverage.pct_measured_all", S["pct_measured_all"])
    R.put("coverage.pct_labelled_all", S["pct_labelled_all"])
    if "conf" in D:
        # DR14 measures, for every labelled gene, the distance to the nearest measured gene of its track; measured
        # genes are at distance 0 from themselves, so the carried labels are the ones with distance > 0
        c = D["conf"]
        S["conf_zero"] = float((c.d_nearest <= 0).mean())
        c = c[c.d_nearest > 0]
        S["conf"] = c
        S["near"] = c.groupby("genome").d_nearest.apply(lambda v: float((v < NEAR_MB).mean())).to_dict()
        S["far"] = c.groupby("genome").d_nearest.apply(lambda v: float((v > FAR_MB).mean())).to_dict()
        S["near_all"], S["far_all"] = float((c.d_nearest < NEAR_MB).mean()), float((c.d_nearest > FAR_MB).mean())
        R.put("coverage.carried_within_near_all", S["near_all"])
        R.put("coverage.carried_beyond_far_all", S["far_all"])
    if "pur" in D:
        pu = D["pur"]
        pu = pu[pu.label.isin(["A", "B"])]
        S["purity_n"] = len(pu)
        S["purity_median"] = float(pu.purity.median())
        S["purity_marginal"] = int((pu.purity < PUR_MARG).sum())
        if "p_binom" in pu:
            S["purity_tested"] = int(pu.p_binom.notna().sum())
            S["purity_sig"] = int((pu.p_binom < 0.05).sum())
        if "unit_sd" in S:
            # expected purity of a correctly called segment under the measured per-gene noise: genes beyond the
            # ambiguity band on the right side / genes beyond it on either side, centres at the predicted ±c
            c0, s0 = abs(PRED_A), S["unit_sd"]
            right, wrong = 1 - norm_sf((c0 - GENE_AMBIG) / s0), norm_sf((c0 + GENE_AMBIG) / s0)
            S["purity_expected"] = right / (right + wrong)
            R.put("purity.expected_if_correct", S["purity_expected"])
        R.put("purity.median", S["purity_median"])
        R.put("purity.below_marker", S["purity_marginal"])
        R.put("purity.binomial_sig", S.get("purity_sig"))


R.guarded("coverage", coverage)


# ============================================================================ validation 1: autocorrelation
def autocorr():
    """Pearson correlation of Δ between measured genes k apart on one chromosome (pairs pooled over chromosomes),
    against two shuffles: within each species (destroys all positional information) and within each chromosome
    (keeps each chromosome's mean, so it asks for blocks shorter than a chromosome). Uses no labels."""
    have_all("gd")
    g = D["gd"].dropna(subset=["mid", "delta"]).sort_values(["genome", "chr", "mid"], kind="mergesort")
    x, pos = g.delta.values.astype(float), g.mid.values.astype(float)
    sid = pd.factorize(g.genome + "|" + g.chr)[0]      # rows of one chromosome are contiguous, codes increasing
    spi, spn = pd.factorize(g.genome)
    lags = np.arange(1, MAX_LAG + 1)
    masks = [sid[:-k] == sid[k:] for k in lags]

    def ac(v):
        out = np.full(len(lags), np.nan)
        for i, k in enumerate(lags):
            m = masks[i]
            if m.sum() > 20:
                out[i] = np.corrcoef(v[:-k][m], v[k:][m])[0, 1]
        return out

    def lag1(v, sel):
        m = masks[0] & sel[:-1]
        return np.corrcoef(v[:-1][m], v[1:][m])[0, 1] if m.sum() > 20 else np.nan

    def shuffled(groups):
        return x[np.lexsort((rng.random(len(x)), groups))]
    obs = ac(x)
    band = {}
    for key, grp in (("species", spi), ("chromosome", sid)):
        null = np.array([ac(shuffled(grp)) for _ in range(A.perm)])
        lo, hi = np.nanpercentile(null, 2.5, 0), np.nanpercentile(null, 97.5, 0)
        above = obs > hi
        reach = int(np.argmin(above)) if not above.all() else len(lags)
        band[key] = dict(lo=lo, hi=hi, reach=reach)
    dmb = np.array([np.median(pos[k:][masks[i]] - pos[:-k][masks[i]]) / 1e6 for i, k in enumerate(lags)])
    per = {}
    nper = max(50, A.perm // 2)
    for j, gn in enumerate(spn):
        sel = spi == j
        hs = np.nanpercentile([lag1(shuffled(spi), sel) for _ in range(nper)], 97.5)
        hc = np.nanpercentile([lag1(shuffled(sid), sel) for _ in range(nper)], 97.5)
        per[gn] = dict(obs=lag1(x, sel), hi_species=hs, hi_chrom=hc, n=int(sel.sum()))
    S["ac"] = dict(lags=lags, obs=obs, band=band, dmb=dmb, per=per)
    R.put("autocorr.lag1", obs[0])
    for key, b in band.items():
        R.put("autocorr.%s_shuffle.lag1_hi" % key, b["hi"][0])
        R.put("autocorr.%s_shuffle.reach_lags" % key, b["reach"])
        R.put("autocorr.%s_shuffle.reach_mb" % key, dmb[b["reach"] - 1] if b["reach"] else float("nan"))
    for gn, v in per.items():
        R.put("autocorr.%s.lag1" % sp(gn), v["obs"])


R.guarded("autocorrelation of Δ", autocorr)


# ============================================================================ validation 2: cross-species dS
def cross_species():
    """For two species at one locus, each carrying an A- and a B-labelled copy: is the first species' A copy
    closer (Drosera-to-Drosera dS, never used to build Δ) to the second species' A copy than to its B copy?
    And likewise for B. Chance is 50% whatever the copy numbers. The first copy of each label in locus_meta
    order represents the species. Control: labels permuted among the copies of each species at each locus."""
    have_all("meta", "pb", "ks")
    K = S["K"]
    lab = D["pb"][D["pb"].label.isin(["A", "B"])].drop_duplicates(["genome", "gene"])[["genome", "gene", "label"]]
    L = D["meta"][D["meta"].genome.isin(DROS)][["locus", "tip", "gene", "genome"]].copy()
    L["order"] = np.arange(len(L))
    L = L.merge(lab, on=["genome", "gene"])
    rank = {g: i for i, g in enumerate(DROS)}

    def run(labels):
        X = L.assign(lb=labels).sort_values("order")
        f = X.groupby(["locus", "genome", "lb"]).tip.first().unstack("lb")
        if "A" not in f or "B" not in f:
            return None
        f = f.dropna(subset=["A", "B"]).reset_index()
        f["r"] = f.genome.map(rank)
        m = f.merge(f, on="locus", suffixes=("1", "2"))
        m = m[m.r1 < m.r2]
        m = look(m, "A1", "A2", "dAA", K)
        m = look(m, "A1", "B2", "dAB", K)
        m = look(m, "B1", "B2", "dBB", K)
        m = look(m, "B1", "A2", "dBA", K)
        m["okA"] = m.dAA.notna() & m.dAB.notna() & (m.dAA != m.dAB)
        m["okB"] = m.dBB.notna() & m.dBA.notna() & (m.dBB != m.dBA)
        m["hitA"], m["hitB"] = m.dAA < m.dAB, m.dBB < m.dBA
        return m
    m = run(L.label.values)
    if m is None or not len(m):
        raise ValueError("no locus carries A and B copies in two species")
    res = {}
    for s in ("A", "B"):
        t = m[m["ok" + s]]
        k, n = int(t["hit" + s].sum()), len(t)
        res[s] = dict(k=k, n=n, f=k / n, ci=wilson(k, n), p=binom_p(k, n), ties=int(((m["d%s%s" % (s, s)] ==
                                                                                       m["d%s%s" % (s, "B" if s == "A" else "A")])).sum()))
    pa, pbb = res["A"]["f"], res["B"]["f"]
    pp = (res["A"]["k"] + res["B"]["k"]) / (res["A"]["n"] + res["B"]["n"])
    z = (pbb - pa) / math.sqrt(pp * (1 - pp) * (1 / res["A"]["n"] + 1 / res["B"]["n"]))
    res["AvsB_p"] = 2 * norm_sf(abs(z))
    # per species pair
    rows = []
    for (g1, g2), t in m.groupby(["genome1", "genome2"]):
        row = {"g1": g1, "g2": g2}
        for s in ("A", "B"):
            u = t[t["ok" + s]]
            row[s + "k"], row[s + "n"] = int(u["hit" + s].sum()), len(u)
        rows.append(row)
    pairs = pd.DataFrame(rows)
    pairs["order"] = pairs.g1.map(rank) * 10 + pairs.g2.map(rank)
    # control: labels permuted within (locus, species)
    grp = pd.factorize(L.locus + "|" + L.genome)[0]
    ctl = {"A": [], "B": []}
    ctl_pairs = []
    for _ in range(N_CTRL):
        i1 = np.lexsort((L.order.values, grp))
        i2 = np.lexsort((rng.random(len(L)), grp))
        perm = np.empty(len(L), dtype=object)
        perm[i1] = L.label.values[i2]
        mc = run(perm)
        for s in ("A", "B"):
            t = mc[mc["ok" + s]]
            ctl[s].append(float(t["hit" + s].mean()))
        ctl_pairs.append(mc.groupby(["genome1", "genome2"]).apply(
            lambda t: pd.Series({"cA": t.loc[t.okA, "hitA"].mean(), "cB": t.loc[t.okB, "hitB"].mean()})))
    cp = pd.concat(ctl_pairs).groupby(level=[0, 1]).mean().reset_index().rename(columns={"genome1": "g1",
                                                                                         "genome2": "g2"})
    pairs = pairs.merge(cp, on=["g1", "g2"], how="left").sort_values("order")
    res["ctl"] = {s: (float(np.mean(v)), float(np.min(v)), float(np.max(v))) for s, v in ctl.items()}
    res["pairs"] = pairs
    res["loci"] = int(m.locus.nunique())
    S["xs"] = res
    for s in ("A", "B"):
        R.put("crossspecies.%s.correct" % s, res[s]["f"])
        R.put("crossspecies.%s.tests" % s, res[s]["n"])
        R.put("crossspecies.%s.p" % s, res[s]["p"])
        R.put("crossspecies.%s.control_mean" % s, res["ctl"][s][0])
    R.put("crossspecies.A_vs_B_p", res["AvsB_p"])
    for _, r in pairs.iterrows():
        for s in ("A", "B"):
            R.put("crossspecies.pair.%s-%s.%s" % (sp(r.g1), sp(r.g2), s), r[s + "k"] / max(r[s + "n"], 1))


R.guarded("cross-species dS test", cross_species)


# ============================================================================ validation 3: Δ vs four-point
def mcnemar():
    have_all("agree")
    a = D["agree"]
    a = a[a.fp.isin(["A", "B"])].copy()
    a["dl"] = np.where(a.delta < 0, "A", "B")
    rows = []
    for g, t in list(a.groupby("genome")) + [("all", a)]:
        aa, ab = int(((t.dl == "A") & (t.fp == "A")).sum()), int(((t.dl == "A") & (t.fp == "B")).sum())
        ba, bb = int(((t.dl == "B") & (t.fp == "A")).sum()), int(((t.dl == "B") & (t.fp == "B")).sum())
        chi = (abs(ab - ba) - 1) ** 2 / (ab + ba) if ab + ba else float("nan")
        p = math.erfc(math.sqrt(chi / 2)) if ok(chi) else float("nan")
        n = aa + ab + ba + bb
        rows.append([isp(g) if g != "all" else "all species", fi(n), ff(100 * (aa + bb) / max(n, 1)), fi(ab), fi(ba),
                     "%.2g" % p if ok(p) else "–"])
        R.put("mcnemar.%s.agreement" % (sp(g) if g != "all" else "all"), (aa + bb) / max(n, 1))
        R.put("mcnemar.%s.p" % (sp(g) if g != "all" else "all"), p)
    S["mcn"] = rows
    S["outside"] = float((D["agree"].fp == "outside").mean())
    R.put("fourpoint.outside_share", S["outside"])


R.guarded("Δ vs four-point", mcnemar)


# ============================================================================ validation 4: the test suite
def run_tests():
    if A.no_tests:
        S["tests"] = None
        return
    rs = A.rscript if A.rscript and os.path.exists(A.rscript) else shutil.which("Rscript")
    if not rs:
        raise FileNotFoundError("Rscript not found (--rscript %s)" % A.rscript)
    if not R.have("DR/DRtest_phasing.R"):
        raise FileNotFoundError("DR/DRtest_phasing.R")
    env = dict(os.environ, SUBG_BASE=os.getcwd())
    pr = subprocess.run([rs, "DR/DRtest_phasing.R"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        universal_newlines=True, timeout=3600, env=env)
    out = pr.stdout
    rows = []                      # (assertion, PASS/FAIL, detail): DRtest prints a FAIL's detail on the next line
    for line in out.splitlines():
        m = re.match(r"^\s*(PASS|FAIL)\s+(.*)$", line)
        if m:
            rows.append([m.group(2).strip(), m.group(1), ""])
        elif rows and rows[-1][1] == "FAIL" and line.startswith("        ") and line.strip():
            rows[-1][2] = (rows[-1][2] + " " + line.strip()).strip()
    S["tests"] = dict(rc=pr.returncode, out=out, rows=rows, passed=[r[0] for r in rows if r[1] == "PASS"],
                      failed=[r[0] for r in rows if r[1] == "FAIL"], rscript=rs)
    R.put("tests.passed", len(S["tests"]["passed"]))
    R.put("tests.failed", len(S["tests"]["failed"]))


R.guarded("test suite DRtest_phasing.R", run_tests)


# ============================================================================ constitution
def constitution():
    have_all("cn")
    cn = D["cn"]
    rows, S["cn_q"] = [], {}
    tr = D.get("trunc")
    for g in DROS:
        x = cn[cn.genome == g]
        if not len(x):
            continue
        q = dict(loci=len(x), A95=x.A.quantile(.95), A99=x.A.quantile(.99), Amax=x.A.max(),
                 B95=x.B.quantile(.95), B99=x.B.quantile(.99), Bmax=x.B.max())
        pa = PRED.get(g)
        q["predA"], q["predB"] = pa, (pa / 2 if pa else None)
        q["at_ceiling"] = float((x.A >= pa).mean()) if pa else float("nan")
        if tr is not None:
            t = tr[tr.genome == g]
            q["A95_worst"] = t.nA_max.quantile(.95) if len(t) else float("nan")
            sh = t[t.nA < pa] if pa else t.iloc[0:0]
            q["short"], q["short_explained"] = len(sh), float((sh.unlab == 0).mean()) if len(sh) else float("nan")
        S["cn_q"][g] = q
        for k, v in q.items():
            if v is not None:
                R.put("constitution.%s.%s" % (sp(g), k), v)
    S["cn_ok"] = {g: (q["A95"] == q["predA"] and q["B95"] == q["predB"]) for g, q in S["cn_q"].items()
                  if q["predA"]}


R.guarded("constitution", constitution)


def regia_control():
    have_all("c05f")
    c = D["c05f"].copy()
    lo, hi = [], []
    for s in c.get("CI", pd.Series([""] * len(c))):        # DR05f writes the binomial CI as "[lo,hi]"
        m = re.findall(r"[-\d.eE]+", str(s))
        lo.append(float(m[0]) if len(m) == 2 else np.nan)
        hi.append(float(m[1]) if len(m) == 2 else np.nan)
    c["lo"], c["hi"] = lo, hi
    rg = c[c.focal == "regia"]
    oth = c[c.focal != "regia"]
    S["c05f"] = c
    if len(rg) and len(oth):
        gap = float(rg.focal_Dionaea.iloc[0] - oth.focal_Dionaea.max())
        S["c05f_gap"] = gap
        S["c05f_rg"] = rg.iloc[0]
        S["c05f_oth"] = (float(oth.focal_Dionaea.min()), float(oth.focal_Dionaea.max()),
                         oth.focal[oth.focal_Dionaea.idxmax()])
        # DR05f's own verdict rule (thresholds read from DR05f_control.R)
        S["c05f_verdict"] = ("regia clearly highest" if gap > F05_CLEAR else
                             "regia not distinguishable from the others" if gap <= F05_NONE else "suggestive only")
        R.put("regia_control.regia", rg.focal_Dionaea.iloc[0])
        R.put("regia_control.gap_to_highest_other", gap)
        R.put("regia_control.verdict", S["c05f_verdict"])


R.guarded("regia four-point control", regia_control)


def riparian():
    if "tracts" in D:
        t = D["tracts"]
        t = t[t.region != "Dionaea"]
        S["tracts_n"] = len(t)
        S["tracts_lab"] = t.label.value_counts().to_dict()
        S["tracts_both"] = int((t.groupby(["genome", "chr"]).label.nunique() > 1).sum())
        R.put("riparian.tracts", len(t))
        R.put("riparian.chromosomes_with_both", S["tracts_both"])
    if "dropped" in D:
        S["dropped_n"] = len(D["dropped"])
        R.put("riparian.braids_dropped", S["dropped_n"])


R.guarded("riparian", riparian)


# ============================================================================ dating
def parse_ctl(path):
    out = {}
    for line in open(path):
        line = line.split("*", 1)[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_tree(text):
    """MCMCtree tree file -> (tips, nodes, root_label). Tips are numbered 1..s in order of appearance and internal
    nodes s+1.. in the order of their opening parenthesis (pre-order, root = s+1): MCMCtree's t_n<k> numbering."""
    lines = [l for l in text.splitlines() if l.strip()]
    if lines and re.match(r"^\s*\d+\s+\d+\s*$", lines[0]):
        lines = lines[1:]
    s = " ".join(lines)
    s = s[:s.index(";")] if ";" in s else s
    pos = [0]
    tips, nodes = [], []

    def ws():
        while pos[0] < len(s) and s[pos[0]].isspace():
            pos[0] += 1

    def label():
        ws()
        out = ""
        while pos[0] < len(s):
            ch = s[pos[0]]
            if ch in "'\"":
                j = s.index(ch, pos[0] + 1)
                out += s[pos[0] + 1:j]
                pos[0] = j + 1
            elif ch in ",():;" or ch.isspace():
                break
            else:
                out += ch
                pos[0] += 1
        ws()
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            while pos[0] < len(s) and s[pos[0]] not in ",();":
                pos[0] += 1
        return out.strip()

    def node():
        ws()
        if s[pos[0]] == "(":
            k = len(nodes)
            nodes.append({"children": [], "label": ""})
            pos[0] += 1
            nodes[k]["children"].append(node())
            ws()
            while s[pos[0]] == ",":
                pos[0] += 1
                nodes[k]["children"].append(node())
                ws()
            if s[pos[0]] != ")":
                raise ValueError("unbalanced parentheses at character %d" % pos[0])
            pos[0] += 1
            nodes[k]["label"] = label()
            return ("n", k)
        name = label()
        m = re.match(r"^([^'<>@]+?)([<>@'].*)?$", name)
        tips.append(m.group(1) if m else name)
        return ("t", len(tips) - 1)
    node()

    def tipset(ch):
        return {tips[ch[1]]} if ch[0] == "t" else set().union(*[tipset(c) for c in nodes[ch[1]]["children"]])
    ns = len(tips)
    for k, nd in enumerate(nodes):
        nd["num"] = ns + 1 + k
        nd["tips"] = tipset(("n", k))
    return tips, nodes


def calib(lab):
    lab = (lab or "").strip().strip("'\"")
    num = r"\s*([-\d.eE+]+)\s*"
    for pat, kind in ((r"^B\(" + num + "," + num, "B"), (r"^>" + num + "<" + num, "B"), (r"^L\(" + num, "L"),
                      (r"^>" + num, "L"), (r"^U\(" + num, "U"), (r"^<" + num, "U")):
        m = re.match(pat, lab)
        if m:
            v = [float(x) for x in m.groups()]
            return (v[0], v[1]) if kind == "B" else (v[0], float("nan")) if kind == "L" else (float("nan"), v[0])
    return None


def tip_label(t):
    """'A' for 'binata_A', '' for 'Nepenthes' (tip names as in DR05h: species_A / species_B)."""
    m = re.match(r"^.+_([AB])$", t)
    return m.group(1) if m else ""


def node_name(nd, tips, nodes):
    """A name for an internal node from its tips: root, A/B split (children all-A and all-B), Droseraceae crown
    (everything but Nepenthes), Drosera crown, subgenome crowns; otherwise the tips themselves."""
    T, S_ = set(tips), set(nd["tips"])
    nep = {t for t in T if "nepenthes" in t.lower()}
    dio = {t for t in T if re.search("dionaea|muscipula", t, re.I)}
    dros = T - nep - dio
    labs = {tip_label(t) for t in T} - {""}
    uniform = len(labs) <= 1

    def child_tips(c):
        return {tips[c[1]]} if c[0] == "t" else set(nodes[c[1]]["tips"])

    def short(t):
        return sp(re.sub(r"_[AB]$", "", t)) if uniform else sp(t)
    if S_ == T:
        return "root"
    kids = [{tip_label(t) for t in child_tips(c)} for c in nd["children"]]
    labelled = {t for t in T if tip_label(t)}
    if (len(kids) == 2 and all(len(k) == 1 and "" not in k for k in kids) and kids[0] != kids[1]
            and labelled <= S_):
        return "A/B split"
    if S_ == T - nep:
        return "Droseraceae crown"
    if S_ == dros:
        return "Drosera crown"
    if not uniform:
        for X in sorted(labs):
            if S_ == {t for t in T if tip_label(t) == X}:
                return "subgenome %s crown" % X
            if S_ == {t for t in dros if tip_label(t) == X}:
                return "Drosera crown, subgenome %s" % X
    order = sorted(S_, key=tips.index)
    if S_ < dros and len(dros - S_) == 1 and uniform:
        return "Drosera without %s" % short(next(iter(dros - S_)))
    return " + ".join(short(t) for t in order)


def ess(x):
    """Effective sample size as coda::effectiveSize: n var(x) / spectral density at 0 from an AR model fitted by
    Yule-Walker with the order chosen by AIC (R's ar(), order.max = 10 log10 n)."""
    x = np.asarray(x, float)
    n = len(x)
    if n < 10 or np.std(x) == 0:
        return 0.0
    om = int(min(n - 1, math.floor(10 * math.log10(n))))
    xc = x - x.mean()
    r = np.array([np.dot(xc[:n - k], xc[k:]) / n for k in range(om + 1)])
    phi = np.zeros((om + 1, om + 1))
    v = np.zeros(om + 1)
    v[0] = r[0]
    for k in range(1, om + 1):
        acc = r[k] - np.dot(phi[k - 1, 1:k], r[k - 1:0:-1])
        kap = acc / v[k - 1]
        phi[k, k] = kap
        phi[k, 1:k] = phi[k - 1, 1:k] - kap * phi[k - 1, k - 1:0:-1]
        v[k] = v[k - 1] * (1 - kap ** 2)
    aic = n * np.log(np.maximum(v, 1e-300)) + 2 * np.arange(om + 1) + 2
    o = int(np.argmin(aic))
    vp = v[o] * n / (n - (o + 1))
    spec = vp / (1 - phi[o, 1:o + 1].sum()) ** 2
    return float(n * np.var(x, ddof=1) / spec)


def hpd(x, p=0.95):
    s = np.sort(np.asarray(x, float))
    n = len(s)
    m = int(math.ceil(p * n))
    if m >= n:
        return float(s[0]), float(s[-1])
    w = s[m - 1:] - s[:n - m + 1]
    i = int(np.argmin(w))
    return float(s[i]), float(s[i + m - 1])


def stem_rep(name):
    m = re.match(r"^(.*_[A-Za-z]+)(\d+)$", name)
    return (m.group(1), int(m.group(2))) if m else (name + "_run", 1)


def dating():
    runs = []
    for ctlp in sorted(glob.glob("DR/dating/*/mcmctree.ctl")):
        R.have(ctlp)
        d = os.path.dirname(ctlp)
        name = os.path.basename(d)
        cfg = parse_ctl(ctlp)
        st, rep = stem_rep(name)
        info = dict(name=name, dir=d, ctl=cfg, ctlp=ctlp, family=name.split("_")[0], stem=st, rep=rep,
                    usedata=int(re.findall(r"-?\d+", cfg.get("usedata", "-1"))[0]))
        tf = cfg.get("treefile", "")
        tpath = tf if os.path.isabs(tf) else os.path.join(d, tf)
        info["treepath"] = tpath
        if tf and os.path.exists(tpath):
            R.have(tpath)
            info["treetext"] = open(tpath).read()
            tips, nodes = parse_tree(info["treetext"])
            info["tips"], info["nodes"] = tips, nodes
            for nd in nodes:
                nd["name"] = node_name(nd, tips, nodes)
            root = nodes[0]
            ch = root["children"]
            nep = [t for t in tips if "nepenthes" in t.lower()]
            info["rooted_ok"] = bool(nep) and any(c[0] == "t" and tips[c[1]] in nep for c in ch)
            info["root_cal_tree"] = root["label"]
            info["calibs"] = {nd["num"]: calib(nd["label"]) for nd in nodes if calib(nd["label"])}
        else:
            info["tips"] = None
            info["rooted_ok"] = None
            info["root_cal_tree"] = ""
            info["calibs"] = {}
        ra = cfg.get("RootAge", "").strip()
        info["root_cal_ctl"] = ra
        info["double_root"] = bool(ra) and bool(info["root_cal_tree"])
        if info["tips"] is not None and ra and calib(ra) and not info["root_cal_tree"]:
            info["calibs"][len(info["tips"]) + 1] = calib(ra)
        info["inBV"] = os.path.exists(os.path.join(d, "in.BV")) and os.path.getsize(os.path.join(d, "in.BV")) > 0
        mf = os.path.join(d, cfg.get("mcmcfile", "mcmc.txt"))
        info["mcmc"] = mf if os.path.exists(mf) and os.path.getsize(mf) > 0 else None
        if info["mcmc"]:
            R.have(mf)
            try:
                m = pd.read_csv(mf, sep=r"\s+").apply(pd.to_numeric, errors="coerce").dropna()  # a running chain
            except Exception as e:                                                     # may end mid-line
                info["mcmc"], info["mcmc_error"] = None, "%s: %s" % (type(e).__name__, e)
                m = None
        if info["mcmc"]:
            tcols = [c for c in m.columns if re.match(r"^t_n\d+$", c)]
            info["samples"] = len(m)
            info["tcols"] = tcols
            info["ages"] = {int(c[3:]): m[c].values * TIME_UNIT_MA for c in tcols}
            info["ess"] = {c: ess(m[c].values) for c in m.columns if c != "Gen"}
            info["ess_min_nodes"] = min(info["ess"][c] for c in tcols) if tcols else float("nan")
            info["ess_min_all"] = min(info["ess"].values()) if info["ess"] else float("nan")
        logs = []
        for pat in ("out.txt", "*.log", "*.out", "*.err", "slurm*"):
            logs += glob.glob(os.path.join(d, pat))
        info["logs"] = sorted(set(logs), key=os.path.getmtime, reverse=True)      # newest first
        runs.append(info)
        R.put("dating.run.%s.samples" % name, info.get("samples", 0))
        if "ess_min_nodes" in info:
            R.put("dating.run.%s.ess_min_nodes" % name, info["ess_min_nodes"])
    if not runs:
        raise FileNotFoundError("no DR/dating/*/mcmctree.ctl")
    S["runs"] = runs
    # replicate groups of posterior chains
    groups = {}
    for r in runs:
        if r["usedata"] != 0 and r["mcmc"] and r.get("tips") is not None and r["tcols"]:
            groups.setdefault((r["family"], r["stem"]), []).append(r)
    G = []
    for (fam, st), rs in sorted(groups.items()):
        rs = sorted(rs, key=lambda r: r["rep"])
        cols = sorted(set.intersection(*[set(r["tcols"]) for r in rs]), key=lambda c: int(c[3:]))
        pooled = {c: sum(r["ess"][c] for r in rs) for c in cols}
        g = dict(family=fam, stem=st, runs=rs, names=" + ".join(r["name"] for r in rs), cols=cols,
                 pooled_min=min(pooled.values()) if pooled else float("nan"), pooled=pooled,
                 rooted_ok=all(r["rooted_ok"] for r in rs), double_root=any(r["double_root"] for r in rs))
        if len(rs) >= 2:
            g["rep_diff"] = max(abs(np.mean(rs[0]["ages"][int(c[3:])]) - np.mean(rs[1]["ages"][int(c[3:])]))
                                for c in cols)
        G.append(g)
        R.put("dating.group.%s.pooled_ess_min" % g["names"].replace(" ", ""), g["pooled_min"])
    S["groups"] = G
    # best replicate group per tree: correctly rooted groups first, then the highest pooled ESS
    best = {}
    for g in G:
        cur = best.get(g["family"])
        if cur is None or (bool(g["rooted_ok"]), g["pooled_min"]) > (bool(cur["rooted_ok"]), cur["pooled_min"]):
            best[g["family"]] = g
    S["best"] = best
    for fam, g in best.items():
        g["citable"] = bool(g["rooted_ok"]) and g["pooled_min"] >= ESS_MIN
        g["why_not"] = ("" if g["citable"] else "mis-rooted" if not g["rooted_ok"] else
                        "pooled ESS %s, below %d" % (fi(g["pooled_min"]), ESS_MIN))
        rs = g["runs"]
        tree = rs[0]
        g["summary"] = []
        for nd in tree["nodes"]:
            k = nd["num"]
            xs_ = [r["ages"][k] for r in rs if k in r["ages"]]
            if not xs_:
                continue                        # tree file and mcmc.txt disagree on the node count
            x = np.concatenate(xs_)
            lo, hi = hpd(x)
            g["summary"].append(dict(num=k, name=nd["name"], mean=float(np.mean(x)), lo=lo, hi=hi,
                                     ess=g["pooled"].get("t_n%d" % k, float("nan")), cal=tree["calibs"].get(k)))
            if g["citable"]:          # ages of unconverged or mis-rooted chains never reach numbers.txt either
                R.put("dating.family.%s.t_n%d.mean_ma" % (fam, k), float(np.mean(x)))
                R.put("dating.family.%s.t_n%d.hpd95_ma" % (fam, k), "%.1f-%.1f" % (lo, hi))
        R.put("dating.family.%s.best" % fam, g["names"])
        R.put("dating.family.%s.citable" % fam, g["citable"])
    S["priors"] = {r["family"]: r for r in runs if r["usedata"] == 0 and r["mcmc"]}
    # logs the job scripts write outside the run folders (slurm's default is the submission directory)
    S["dating_logs"] = sorted(set(sum((glob.glob(os.path.join("DR/dating", p)) for p in
                                       ("*.out", "*.err", "*.log", "slurm*", "logs/*")), [])),
                              key=os.path.getmtime, reverse=True)       # newest first


R.guarded("dating", dating)


# ============================================================================ known issues, provenance
def known_issues():
    out = []
    if "q06" in D:
        q = D["q06"].fillna({"mid_n": 0, "low_n": 0})
        bad = q[(q.top_n + q.mid_n + q.low_n) != q.n]
        S["q06_bad"] = bad
        R.put("dr06.rows_not_summing", len(bad))
    S["dr24"] = os.path.exists(PIPE_FIGS["dr24"])
    S["modelfit"] = os.path.exists(PIPE_FIGS["modelfit"])
    # which script draws the all-13 chronogram, and which run folders does it read?
    fig = os.path.basename(PIPE_FIGS["dr24"])
    S["dr24_src"] = []
    for f in sorted(glob.glob("DR/*.R") + glob.glob("DR/*.py") + glob.glob("DR/dating/*.R") +
                    glob.glob("DR/dating/*.py") + glob.glob("DR/dating/*.sh")):
        try:
            txt = open(f, errors="replace").read()
        except OSError:
            continue
        if fig in txt and os.path.basename(f) not in ("DRreport.py", "report_kit.py"):
            runs = sorted(set(re.findall(r"(all13\w*)/", txt)))
            S["dr24_src"].append((f, runs))
    return out


R.guarded("known issues", known_issues)


def provenance():
    def git(*a):
        pr = subprocess.run(["git"] + list(a), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            universal_newlines=True)
        return [l for l in pr.stdout.splitlines() if l.strip()] if pr.returncode == 0 else None
    # the same patterns the delivery block commits; shell globs (not git pathspec magic, which old git lacks)
    pats = ("DR/*.R", "DR/*.py", "DR/*.sh", "DR/*.md", "DR/*.sbatch", "DR/*/*.R", "DR/*/*.py", "DR/*/*.sh",
            "DR/*/*.sbatch", "DR/dating/*/mcmctree.ctl", "DR/dating/*/tree.trees")
    cand = sorted({f for p in pats for f in glob.glob(p) if os.path.isfile(f)})
    tracked = set()
    for i in range(0, len(cand), 500):
        t = git("ls-files", "--", *cand[i:i + 500])
        if t is None:
            raise RuntimeError("git ls-files failed (not a git work tree?)")
        tracked |= set(t)
    mod = git("status", "--porcelain", "--", "DR")
    mine = {"DR/DRreport.py", "DR/report_kit.py"}       # the report's own code, committed by its build block
    S["untracked"] = sorted(set(cand) - tracked - mine)
    S["modified"] = [l for l in (mod or []) if not l.startswith("??") and not l.startswith("!!")]
    R.put("provenance.untracked_code_files", len(S["untracked"]))


R.guarded("provenance", provenance)

# ============================================================================ figures
FIG = {}


def fig_delta():
    gd, seg = D["gd"], D["seg"]
    sps = [g for g in DROS if g in set(gd.genome)]
    fig, axs = plt.subplots(2, len(sps), figsize=(2.3 * len(sps), 4.6), sharex="row")
    axs = np.atleast_2d(axs)
    for j, g in enumerate(sps):
        ax = axs[0, j]
        x = gd.delta[gd.genome == g]
        ax.hist(x.clip(-1, 1), bins=np.linspace(-1, 1, 51), color="#888780")
        for v, s in ((PRED_A, "A"), (PRED_B, "B")):
            if ok(v):
                ax.axvline(v, color=COL[s], ls="--", lw=1)
        ax.set_title("%s\n%s measured genes" % (sp(g), fi(len(x))), fontsize=8, loc="left")
        if j == 0:
            ax.set_ylabel("genes")
        ax.set_xlabel("Δ per gene", fontsize=8)
        ax = axs[1, j]
        t = seg[(seg.genome == g) & (seg.n >= MINSEG)]
        bins = np.linspace(-0.5, 0.5, 41)
        ax.hist([t.med_delta[t.label == l].clip(-0.5, 0.5) for l in ("A", "B", "ambiguous")], bins=bins,
                stacked=True, color=[COL["A"], COL["B"], COL["ambiguous"]], label=["A", "B", "ambiguous"])
        for v in (SEG_LO, SEG_HI):
            ax.axvline(v, color="#444441", ls=":", lw=1)
        ax.set_title("%s segments" % fi(len(t)), fontsize=8, loc="left")
        if j == 0:
            ax.set_ylabel("segments")
            ax.legend(frameon=False, fontsize=7)
        ax.set_xlabel("segment median Δ", fontsize=8)
    fig.tight_layout()
    return fig


def fig_noise():
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), gridspec_kw={"width_ratios": [1, 1.3]})
    ax = axs[0]
    x = S["nep"]
    ax.hist(np.clip(x, -1.5, 1.5), bins=np.linspace(-1.5, 1.5, 61), color="#888780")
    ax.axvline(0, color="#2C2C2A", lw=1)
    ax.axvline(S["nep_median"], color=COL["B"], ls="--", lw=1)
    ax.set_xlabel("Δ for the Nepenthes copy (true value 0)")
    ax.set_ylabel("loci")
    ax.set_title("Nepenthes control: %s loci, median %+.3f, sd %.2f" % (fi(len(x)), S["nep_median"], S["nep_sd"]),
                 loc="left", fontsize=9)
    ax = axs[1]
    u = S["units"]
    sps = [g for g in DROS if g in set(u.genome)]
    for j, g in enumerate(sps):
        v = u.sd[u.genome == g].values
        ax.scatter(j + rng.uniform(-0.18, 0.18, len(v)), v, s=10, color=SPC[g], alpha=0.7, lw=0)
        ax.plot([j - 0.3, j + 0.3], [S["unit_sd_by_sp"][g]] * 2, color="#2C2C2A", lw=1.5)
    ax.axhline(NOISE_PRED, color=COL["B"], ls="--", lw=1)
    ax.text(len(sps) - 0.5, NOISE_PRED, " predicted in DR02 (%s)" % ff(NOISE_PRED, 2), color=COL["B"], fontsize=7,
            va="bottom", ha="right")
    ax.set_xticks(range(len(sps)))
    ax.set_xticklabels([sp(g) for g in sps])
    ax.set_ylabel("sd of Δ within a unit")
    ax.set_title("spread of Δ inside one species × region × chromosome (≥ %d genes); bar = pooled" % MIN_UNIT,
                 loc="left", fontsize=9)
    fig.tight_layout()
    return fig


def fig_coverage():
    cov = S["cov"]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))
    ax = axs[0]
    y = np.arange(len(cov))[::-1]
    m = cov.pct_measured.values
    pr = (cov.pct_labelled - cov.pct_measured).clip(lower=0).values
    ax.barh(y, m, color="#2C2C2A")
    ax.barh(y, pr, left=m, color="#B4B2A9")
    ax.set_yticks(y)
    ax.set_yticklabels([sp(g) for g in cov.index])
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of the species' genes")
    ax.set_title("how each gene got its label: black = own Δ measured,\ngrey = label carried from its block",
                 loc="left", fontsize=9)
    ax = axs[1]
    c = S.get("conf")
    if c is not None:
        for g in DROS:
            v = np.sort(c.d_nearest[c.genome == g].values)
            if len(v):
                ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=SPC[g], lw=1.2, label=sp(g))
        for v in (NEAR_MB, FAR_MB):
            ax.axvline(v, color="#444441", ls=":", lw=1)
        ax.set_xscale("log")
        ax.set_xlabel("distance to the nearest measured gene, Mb (log scale)")
        ax.set_ylabel("share of carried labels")
        ax.legend(frameon=False, fontsize=7, loc="upper left")
        ax.set_title("how far each label was carried (genes without their own Δ;\ndotted: %s kb and %s Mb)"
                     % (fi(1000 * NEAR_MB), ff(FAR_MB, 0)), loc="left", fontsize=9)
    fig.tight_layout()
    return fig


def fig_autocorr():
    a = S["ac"]
    bs, bc = a["band"]["species"], a["band"]["chromosome"]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.5), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axs[0]
    ax.fill_between(a["lags"], bs["lo"], bs["hi"], color="#E8E6DF", lw=0, label="shuffled within species (95%)")
    ax.fill_between(a["lags"], bc["lo"], bc["hi"], color="#B4B2A9", lw=0, alpha=0.6,
                    label="shuffled within chromosome (95%)")
    ax.plot(a["lags"], a["obs"], color="#2C2C2A", lw=1.5, marker="o", ms=2.5, label="observed")
    ax.axhline(0, color="#888780", lw=0.6)
    ax.set_xlabel("how many measured genes apart on one chromosome (lag)")
    ax.set_ylabel("correlation of their Δ")
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, -0.2), ncol=3)

    def reach_txt(b):
        return ("up to lag %d (%s Mb)" % (b["reach"], ff(a["dmb"][b["reach"] - 1], 1)) if b["reach"]
                else "not even at lag 1")
    ax.set_title("observed above species shuffles %s;\nabove chromosome shuffles %s" % (reach_txt(bs), reach_txt(bc)),
                 loc="left", fontsize=9)
    ax = axs[1]
    per = [(g, a["per"][g]) for g in DROS if g in a["per"]]
    y = np.arange(len(per))[::-1]
    ax.barh(y, [p[1]["obs"] for p in per], color=[SPC[p[0]] for p in per], height=0.6)
    ax.scatter([p[1]["hi_species"] for p in per], y, marker="|", s=260, color="#2C2C2A", lw=1.6)
    ax.scatter([p[1]["hi_chrom"] for p in per], y, marker="|", s=260, color="#888780", lw=1.6)
    ax.set_yticks(y)
    ax.set_yticklabels(["%s (%s genes)" % (sp(p[0]), fi(p[1]["n"])) for p in per], fontsize=8)
    ax.axvline(0, color="#888780", lw=0.6)
    ax.set_xlabel("lag-1 correlation of Δ")
    ax.set_title("per species; ticks = 97.5% of shuffles\n(black within species, grey within chromosome)",
                 loc="left", fontsize=9)
    fig.tight_layout()
    return fig


def fig_cross():
    xs = S["xs"]
    p = xs["pairs"]
    fig, ax = plt.subplots(figsize=(8.5, 0.42 * len(p) + 1.3))
    y = np.arange(len(p))[::-1]
    for s, dy in (("A", 0.13), ("B", -0.13)):
        f = p[s + "k"] / p[s + "n"].clip(lower=1)
        ci = np.array([wilson(k, n) for k, n in zip(p[s + "k"], p[s + "n"])])
        ax.errorbar(100 * f, y + dy, xerr=[100 * (f - ci[:, 0]), 100 * (ci[:, 1] - f)], fmt="o", ms=5, color=COL[s],
                    ecolor=COL[s], elinewidth=1, capsize=0, label="%s copies grouped with %s" % (s, s))
        ax.scatter(100 * p["c" + s], y + dy, marker="x", s=18, color=COL["control"],
                   label="labels permuted" if s == "A" else None)
    ax.axvline(50, color="#444441", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(["%s – %s  (%s / %s)" % (sp(a), sp(b), fi(na), fi(nb))
                        for a, b, na, nb in zip(p.g1, p.g2, p.An, p.Bn)], fontsize=8)
    ax.set_xlabel("% of tests where the same-label copy of the other species is the closer one (95% CI)")
    ax.set_xlim(0, 100)
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    ax.set_title("cross-species test, per species pair (tests A / B in brackets); dashed = chance", loc="left",
                 fontsize=9)
    fig.tight_layout()
    return fig


def fig_constitution():
    cn = D["cn"]
    sps = [g for g in DROS if g in S["cn_q"]]
    fig, axs = plt.subplots(1, len(sps), figsize=(2.2 * len(sps), 2.9), sharey=True)
    axs = np.atleast_1d(axs)
    top_a = int(min(8, cn.A.max())) if len(cn) else 4
    for ax, g in zip(axs, sps):
        x = cn[cn.genome == g]
        c = x.groupby(["A", "B"]).size().reset_index(name="n")
        c = c[(c.A <= 8) & (c.B <= 6)]
        ax.scatter(c.B, c.A, s=6 + 900 * c.n / c.n.max(), color=SPC[g], alpha=0.6, lw=0)
        big = c.loc[c.n.idxmax()]
        ax.annotate(fi(big.n), (big.B, big.A), fontsize=6.5, ha="center", va="center", color="white")
        q = S["cn_q"][g]
        if q["predA"]:
            ax.scatter([q["predB"]], [q["predA"]], marker="s", s=90, facecolor="none", edgecolor="#2C2C2A", lw=1.2)
        ax.plot([0, 3], [0, 6], color="#B4B2A9", ls="--", lw=0.8)
        ax.set_xlim(-0.5, 4.5)
        ax.set_ylim(-0.5, max(top_a, 4) + 0.5)
        ax.set_xlabel("B copies at a locus")
        ax.set_title("%s: 95%% of loci have\n≤ %s A and ≤ %s B copies" % (sp(g), ff(q["A95"], 0), ff(q["B95"], 0)),
                     loc="left", fontsize=8)
    axs[0].set_ylabel("A copies at a locus")
    fig.tight_layout()
    return fig


def fig_regia():
    c = S["c05f"].sort_values("focal_Dionaea")
    fig, ax = plt.subplots(figsize=(6.5, 0.35 * len(c) + 1.2))
    y = np.arange(len(c))
    f = c.focal_Dionaea.values
    spc = {sp(k): v for k, v in SPC.items()}
    ax.barh(y, f, color=[spc.get(v, "#888780") if v == "regia" else "#D3D1C7" for v in c.focal])
    ax.errorbar(f, y, xerr=[f - c.lo.values, c.hi.values - f], fmt="none", ecolor="#2C2C2A", lw=1)
    oth = c[c.focal != "regia"].focal_Dionaea
    if len(oth):
        ax.axvline(oth.max(), color="#888780", ls=":", lw=1)
        if ok(F05_CLEAR):
            ax.axvline(oth.max() + F05_CLEAR, color="#444441", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(["%s (%s loci)" % (a, fi(b)) for a, b in zip(c.focal, c.loci)])
    ax.set_xlim(0, max(0.6, float(np.nanmax(c.hi.values)) + 0.05))
    ax.set_xlabel("share of loci where the focal species' copy pairs with Dionaea's (95% CI)")
    ax.set_title("dotted: highest non-regia species; dashed: + %s, DR05f's bar for 'regia clearly highest'"
                 % ff(F05_CLEAR, 2), loc="left", fontsize=8.5)
    fig.tight_layout()
    return fig


def fam_color(fam):
    """Subgenome trees take the subgenome colours used everywhere in the report; other trees are grey."""
    return COL.get(fam, "#888780") if fam in ("A", "B") else "#888780"


def fig_ess():
    G = S["groups"]
    rows = []                 # (label, min ESS, family, pooled?, mis-rooted?)
    for g in G:
        mis = not g["rooted_ok"]
        for r in g["runs"]:
            rows.append((r["name"] + (" (mis-rooted)" if mis else ""), r["ess_min_nodes"], g["family"], False, mis))
        if len(g["runs"]) > 1:
            rows.append((g["names"] + " pooled" + (" (mis-rooted)" if mis else ""), g["pooled_min"], g["family"],
                         True, mis))
    fig, axs = plt.subplots(1, 2, figsize=(11, 0.26 * len(rows) + 1.8), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axs[0]
    y = np.arange(len(rows))[::-1]
    for yy, r in zip(y, rows):
        ax.barh(yy, r[1], color="white" if r[4] else fam_color(r[2]), alpha=1.0 if r[3] else 0.5,
                edgecolor=fam_color(r[2]), hatch="////" if r[4] else None, lw=0.8)
    ax.axvline(ESS_MIN, color="#A32D2D", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax.set_xlabel("smallest ESS over node ages (pooled = sum over replicate chains)")
    ax.set_title("convergence: dashed = ESS %d; dark = replicates pooled; hatched = mis-rooted tree" % ESS_MIN,
                 loc="left", fontsize=9)
    # replicate agreement as a relative difference, so no age of an unconverged chain is shown
    ax = axs[1]
    reps = [g for g in G if len(g["runs"]) >= 2]
    yr = np.arange(len(reps))[::-1]
    for yy, g in zip(yr, reps):
        a, b = g["runs"][0], g["runs"][1]
        d = []
        for c in g["cols"]:
            k = int(c[3:])
            d.append(100 * abs(np.mean(a["ages"][k]) - np.mean(b["ages"][k])) /
                     np.mean(np.concatenate([a["ages"][k], b["ages"][k]])))
        mis = not g["rooted_ok"]
        ax.scatter(d, yy + rng.uniform(-0.15, 0.15, len(d)), s=16, facecolor="none" if mis else fam_color(g["family"]),
                   edgecolor=fam_color(g["family"]), lw=0.9)
    ax.set_yticks(yr)
    ax.set_yticklabels([g["names"] + (" (mis-rooted)" if not g["rooted_ok"] else "") for g in reps], fontsize=7)
    ax.set_xlim(left=0)
    ax.set_xlabel("difference between the two replicates' posterior\nmean ages, % of the node's age (one dot per node)")
    ax.set_title("do replicate chains agree on every node?", loc="left", fontsize=9)
    fig.tight_layout()
    return fig


def draw_chrono(ax, g, color):
    tree = g["runs"][0]
    tips, nodes = tree["tips"], tree["nodes"]
    summ = {s["num"]: s for s in g["summary"]}
    ns = len(tips)
    order = []

    def walk(ch):
        if ch[0] == "t":
            order.append(ch[1])
        else:
            for c in nodes[ch[1]]["children"]:
                walk(c)
    walk(("n", 0))
    ypos = {t: len(order) - 1 - i for i, t in enumerate(order)}

    def y_of(ch):
        if ch[0] == "t":
            return ypos[ch[1]]
        ys = [y_of(c) for c in nodes[ch[1]]["children"]]
        return (min(ys) + max(ys)) / 2

    def age(ch):
        return 0.0 if ch[0] == "t" else summ[ns + 1 + ch[1]]["mean"]
    for k, nd in enumerate(nodes):
        x0 = -age(("n", k))
        ys = [y_of(c) for c in nd["children"]]
        ax.plot([x0, x0], [min(ys), max(ys)], color="#2C2C2A", lw=1)
        for c in nd["children"]:
            ax.plot([x0, -age(c)], [y_of(c)] * 2, color="#2C2C2A", lw=1)
        s = summ[ns + 1 + k]
        yy = y_of(("n", k))
        if s.get("cal"):
            lo, hi = s["cal"]
            lo = 0 if not ok(lo) else lo * TIME_UNIT_MA
            hi = lo + 1 if not ok(hi) else hi * TIME_UNIT_MA
            ax.add_patch(plt.Rectangle((-hi, yy - 0.32), hi - lo, 0.64, color="#F1EFE8", ec="#B4B2A9", lw=0.6, zorder=0))
        ax.plot([-s["hi"], -s["lo"]], [yy, yy], color=color, lw=6, alpha=0.45, solid_capstyle="butt")
        ax.text(x0, yy + 0.12, " %s" % ff(s["mean"], 0), fontsize=7, va="bottom", ha="left", color="#444441")
    root = summ[ns + 1]
    ax.plot([-root["mean"] * 1.04, -root["mean"]], [y_of(("n", 0))] * 2, color="#2C2C2A", lw=1)
    for t, yv in ypos.items():
        ax.text(1.5, yv, sp(tips[t]), va="center", fontsize=8)
    ax.set_xlim(-root["hi"] * 1.08, root["hi"] * 0.28)
    ax.set_ylim(-0.8, len(tips) - 0.2)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % abs(v) if v <= 0 else ""))
    ax.set_xlabel("Ma")


def fam_label(f):
    """'subgenome A' for the subgenome trees, 'the all13 runs (13-tip tree)' for the others."""
    if f in ("A", "B"):
        return "subgenome " + f
    g = S.get("best", {}).get(f)
    n = len(g["runs"][0]["tips"]) if g and g["runs"][0].get("tips") else None
    return "the %s runs%s" % (f, " (%d-tip tree)" % n if n else "")


def stalled():
    """Families of posterior runs (usedata != 0) of which none has produced samples."""
    out = []
    fams = {}
    for r in S.get("runs", []):
        if r["usedata"] != 0:
            fams.setdefault(r["family"], []).append(r)
    for f, rs in sorted(fams.items()):
        if any(r["mcmc"] for r in rs):
            continue
        tips = next((r["tips"] for r in rs if r.get("tips") is not None), None) or []
        # the nearest run that worked: the first run of a citable group, else any run with samples
        ok_runs = [g["runs"][0] for g in S.get("best", {}).values() if g.get("citable")] or \
            [r for r in S["runs"] if r["mcmc"] and r["usedata"] != 0]
        ref = ok_runs[0] if ok_runs else None
        diff = []
        if ref:
            r0 = rs[0]
            for k in sorted(set(r0["ctl"]) | set(ref["ctl"])):
                if k not in ("seed", "mcmcfile", "outfile") and r0["ctl"].get(k) != ref["ctl"].get(k):
                    diff.append("%s = %s (%s: %s)" % (k, r0["ctl"].get(k, "absent"), ref["name"],
                                                       ref["ctl"].get(k, "absent")))
            if r0.get("tips") is not None and ref.get("tips") is not None and len(r0["tips"]) != len(ref["tips"]):
                diff.append("%d tips (%s: %d)" % (len(r0["tips"]), ref["name"], len(ref["tips"])))
        out.append(dict(fam=f, runs=rs, ntip=len(tips) if tips else float("nan"),
                        rooted=all(r["rooted_ok"] for r in rs),
                        has_AB={tip_label(t) for t in tips} >= {"A", "B"},
                        double=[r["name"] for r in rs if r["double_root"]],
                        no_inbv=[r["name"] for r in rs if r["usedata"] == 2 and not r["inBV"]],
                        ref=ref, diff=diff))
    return out


def chrono_families():
    """Trees to draw: the subgenome trees and the all-13 tree, in that order, when they have a best group."""
    return [f for f in ["A", "B"] + sorted(k for k in S.get("best", {}) if k not in ("A", "B"))
            if f in S.get("best", {})]


def fig_chrono():
    """Chronograms only for citable trees (correct root, pooled ESS >= ESS_MIN); the rest are named, not drawn."""
    fams = [f for f in chrono_families() if S["best"][f]["citable"]]
    if not fams:
        return None
    fig, axs = plt.subplots(1, len(fams), figsize=(5.4 * len(fams), 3.6), squeeze=False)
    for ax, f in zip(axs[0], fams):
        g = S["best"][f]
        draw_chrono(ax, g, fam_color(f))
        ax.set_title("%s: %s (pooled ESS %s)" % ("subgenome " + f if f in ("A", "B") else f, g["names"],
                                                 fi(g["pooled_min"])), loc="left", fontsize=9)
    fig.tight_layout()
    return fig


def fig_prior():
    fams = [f for f in chrono_families() if f in S["priors"] and S["best"][f]["citable"]]
    cells = []
    for f in fams:
        g, pr = S["best"][f], S["priors"][f]
        for s in g["summary"]:
            if s["cal"] or s["name"] == "root":
                cells.append((f, s, pr))
    if not cells:
        return None
    nc = max(sum(1 for c in cells if c[0] == f) for f in fams)
    fig, axs = plt.subplots(len(fams), nc, figsize=(3.2 * nc, 2.4 * len(fams)), squeeze=False)
    for i, f in enumerate(fams):
        mine = [c for c in cells if c[0] == f]
        for j in range(nc):
            ax = axs[i, j]
            if j >= len(mine):
                ax.axis("off")
                continue
            _, s, pr = mine[j]
            post = np.concatenate([r["ages"][s["num"]] for r in S["best"][f]["runs"]])
            prior = pr["ages"].get(s["num"])
            lo = min(post.min(), prior.min() if prior is not None else post.min())
            hi = max(post.max(), prior.max() if prior is not None else post.max())
            bins = np.linspace(lo, hi, 60)
            if prior is not None:
                ax.hist(prior, bins=bins, density=True, color="#B4B2A9", alpha=0.7, label="prior (%s)" % pr["name"])
            ax.hist(post, bins=bins, density=True, color=fam_color(f), alpha=0.6, label="posterior")
            if s["cal"]:
                for v in s["cal"]:
                    if ok(v):
                        ax.axvline(v * TIME_UNIT_MA, color="#444441", ls=":", lw=1)
            ax.set_title("%s: %s (t_n%d)" % ("subgenome " + f if f in ("A", "B") else f, s["name"], s["num"]),
                         loc="left", fontsize=8)
            ax.set_xlabel("Ma")
            ax.set_yticks([])
            if j == 0:
                ax.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    return fig


for name, fn, need_keys in (("delta", fig_delta, ("gd", "seg")), ("noise", fig_noise, ()),
                            ("coverage", fig_coverage, ()), ("autocorr", fig_autocorr, ()),
                            ("cross", fig_cross, ()), ("constitution", fig_constitution, ("cn",)),
                            ("regia", fig_regia, ()), ("ess", fig_ess, ()), ("chrono", fig_chrono, ()),
                            ("prior", fig_prior, ())):
    def make(fn=fn, need_keys=need_keys):
        have_all(*need_keys)
        return fn()
    FIG[name] = R.guarded("figure " + name, make)


def pred_text():
    """DR22's predicted A : B ceilings, grouped: '2 A : 1 B for regia, binata, ...; 4 A : 2 B for capensis'."""
    by = {}
    for g in DROS:
        if PRED.get(g):
            by.setdefault(PRED[g], []).append(sp(g))
    return ("; ".join("%d A : %s B for %s" % (a, ff(a / 2, 0), ", ".join(v)) for a, v in sorted(by.items())) +
            "; the A ceiling is DR22's, the B ceiling half of it under the 2A : 1B constitution")


def show(name, look, width=100, absent=None):
    """Embed figure `name` with its what-to-look-for line. A figure function that returns None on purpose (nothing
    to draw) gets the `absent` sentence instead; a figure that failed points to the list at the top."""
    f = FIG.get(name)
    if f is None:
        failed = any(p[0] == "figure " + name for p in R.problems)
        if absent and not failed:
            R.p(absent)
        else:
            R.parts.append("<p class='missing'>figure not built: %s (see the list at the top)</p>" % name)
        return
    R.figure(f, "fig_" + name, width)
    R.look(look)


# ============================================================================ the page
def summary():
    R.h2("0. Summary")
    cov, xs, ac = S.get("cov"), S.get("xs", {}), S.get("ac", {})
    tests = S.get("tests")
    sent = ("Five <i>Drosera</i> species were split into the two ancestral subgenomes A and B by comparing each gene "
            "copy with the two homeologous copies of <i>Dionaea muscipula</i> (Dio_A, Dio_B), with <i>Nepenthes "
            "gracilis</i> as outgroup. Of %s <i>Drosera</i> genes, %s%% carry a measured Δ; segments and blocks built "
            "from those carry a label to %s%%. " % (fi(S.get("dros_genes")), ff(S.get("pct_measured_all")),
                                                     ff(S.get("pct_labelled_all"))))
    if "A" in xs:
        good = all(xs[s]["p"] < 1e-3 and xs[s]["f"] > 0.5 and xs[s]["f"] > xs["ctl"][s][2] for s in ("A", "B"))
        sent += ("%s Drosera-to-Drosera dS, which Δ never uses, groups A with A in %s%% of %s tests and B with B in "
                 "%s%% of %s (chance 50%%; with labels permuted %s%% and %s%%). " % (
                     "The labels pass a test independent of how they were made:" if good else
                     "The labels do not clearly pass the independent test:",
                     ff(100 * xs["A"]["f"]), fi(xs["A"]["n"]), ff(100 * xs["B"]["f"]), fi(xs["B"]["n"]),
                     ff(100 * xs["ctl"]["A"][0]), ff(100 * xs["ctl"]["B"][0])))
    if ac:
        hs = ac["band"]["species"]["hi"][0]
        sent += ("Δ is %s between neighbouring measured genes on a chromosome (lag-1 correlation %s; 97.5%% of "
                 "shuffles %s). " % ("correlated" if ac["obs"][0] > hs else "not detectably correlated",
                                     ff(ac["obs"][0], 3), ff(hs, 3)))
    if S.get("cn_q"):
        sent += ("95%% of loci carry at most %s copies (DR22 predicts %s). " % (
            "; ".join("%s %s A and %s B" % (sp(g), ff(q["A95"], 0), ff(q["B95"], 0)) for g, q in S["cn_q"].items()),
            pred_text()))
    if S.get("best"):
        parts = []
        for f in chrono_families():
            g = S["best"][f]
            what = fam_label(f)
            parts.append("%s: %s (%s)" % (what, "converged" if g["citable"] else "not citable",
                                          "%s, pooled ESS %s" % (g["names"], fi(g["pooled_min"])) if g["citable"]
                                          else g["why_not"]))
        for x in stalled():
            parts.append("the %s runs (%s-tip tree%s%s) have produced no samples" % (
                x["fam"], fi(x["ntip"]), ", correctly rooted" if x["rooted"] else "",
                ", which would date the A/B split" if x["has_AB"] else ""))
        sent += "Dating: %s. " % "; ".join(parts)
    R.p(sent)
    rows = [["Δ measured, % of genes"] + [ff(cov.pct_measured[g]) if cov is not None and g in cov.index else "–"
                                          for g in DROS],
            ["label (measured or propagated), % of genes"] + [ff(cov.pct_labelled[g]) if cov is not None and
                                                               g in cov.index else "–" for g in DROS]]
    if "seg_counts" in S:
        sc = S["seg_counts"]
        rows.append(["segments called A / B / ambiguous / too small"] + [
            " / ".join(fi(sc.loc[g].get(l, 0)) if g in sc.index else "–" for l in ("A", "B", "ambiguous", "unresolved"))
            for g in DROS])
    if S.get("unit_sd_by_sp"):
        rows.append(["noise of Δ per gene (sd within units)"] + [ff(S["unit_sd_by_sp"].get(g), 2) for g in DROS])
    if S.get("ac", {}).get("per"):
        pr = S["ac"]["per"]
        rows.append(["lag-1 correlation of Δ (97.5% of shuffles within species)"] + [
            "%s (%s)" % (ff(pr[g]["obs"], 3), ff(pr[g]["hi_species"], 3)) if g in pr else "–" for g in DROS])
    if S.get("cn_q"):
        rows.append(["A / B copies per locus, 95th percentile (predicted)"] + [
            "%s / %s (%s / %s)" % (ff(S["cn_q"][g]["A95"], 0), ff(S["cn_q"][g]["B95"], 0),
                                   fi(S["cn_q"][g]["predA"]), fi(S["cn_q"][g]["predB"])) if g in S["cn_q"] else "–"
            for g in DROS])
    R.table([""] + [isp(g) for g in DROS], rows)
    # status per item
    st = []
    if "A" in xs:
        ok_lab = all(xs[s]["p"] < 1e-3 and xs[s]["f"] > 0.5 and xs[s]["f"] > xs["ctl"][s][2] for s in ("A", "B"))
        nf = len(tests["failed"]) if tests else 0
        st.append(["A/B labels (DR02, DR02c)",
                   ("supported" if not nf else "supported; %d assertion%s failing" % (nf, "s" if nf > 1 else ""))
                   if ok_lab else "not shown",
                   "cross-species test (4b) %s chance and the permuted control for both subgenomes; DRtest_phasing.R "
                   "%s" % ("above" if ok_lab else "not clearly above",
                           "%d of %d assertions pass%s" % (len(tests["passed"]), len(tests["rows"]),
                                                           " (failing: %s)" % "; ".join(tests["failed"]) if nf else "")
                           if tests else "not run")])
    if S.get("cn_ok"):
        st.append(["constitution (DR21, DR22)", "consistent" if all(S["cn_ok"].values()) else "mixed",
                   "95th-percentile copy numbers equal the predicted ceilings for %d of %d species (%s)" % (
                       sum(S["cn_ok"].values()), len(S["cn_ok"]),
                       ", ".join(sp(g) for g, v in S["cn_ok"].items() if not v) + " differ" if not all(
                           S["cn_ok"].values()) else "all")])
    if "c05f_gap" in S:
        st.append(["regia with Dionaea (DR05f)", S["c05f_verdict"],
                   "regia minus the highest other species: %+.3f (DR05f's rule: > %s clear, ≤ %s none, between "
                   "suggestive)" % (S["c05f_gap"], ff(F05_CLEAR, 2), ff(F05_NONE, 2))])
    for f in chrono_families():
        g = S["best"][f]
        st.append(["dating, %s" % fam_label(f),
                   "citable" if g["citable"] else "not citable: " + g["why_not"],
                   "best replicate group %s: pooled ESS %s on its worst node" % (g["names"], fi(g["pooled_min"]))])
    if S.get("runs"):
        for x in stalled():
            st.append(["dating, the %s runs (%s-tip tree%s)" % (x["fam"], fi(x["ntip"]),
                                                               ", A/B split" if x["has_AB"] else ""),
                       "blocked: no samples",
                       "%d runs, none with samples; rooted on Nepenthes: %s; root calibrated in both ctl and tree: %s"
                       % (len(x["runs"]), "yes" if x["rooted"] else "no", ", ".join(x["double"]) or "none")])
    R.h3("Status")
    R.table(["item", "status", "basis (computed)"], st)


def material_section():
    R.h2("1. Material and terms")
    if "material" in S:
        R.table(["species", "role", "genes (GENESPACE, array representatives)", "sequences with genes",
                 "copies in the locus set", "loci with a copy"], S["material"])
    R.p("<b>Terms, from first principles.</b> A <b>locus</b> is one set of homologous genes across the seven "
        "genomes (one Nepenthes anchor gene or one Dionaea homeolog pair, DR00); a <b>copy</b> (tip) is one gene at "
        "one locus. <b>Dio_A</b> and <b>Dio_B</b> are the two Dionaea copies of a gene, one on each chromosome of a "
        "homeologous pair; <b>A</b> is by definition the side of the pair that kept more genes "
        "(<code>retained_more</code> in fractionation_by_chrpair.csv), <b>B</b> the other. For one Drosera copy, "
        "gene_Drosera, <b>Δ = [d(gene_Drosera, Dio_A) − d(gene_Drosera, Dio_B)] / d(Dio_A, Dio_B)</b>, with d the "
        "synonymous divergence dS (yn00): negative means closer to Dio_A, so A-derived. A <b>region</b> is the "
        "Nepenthes chromosome a Dionaea pair maps to (eight, one to one). A <b>track</b> is one species × region × "
        "chromosome series of measured copies; a <b>segment</b> is a piece of a track after changepoint splitting, "
        "called A if its median Δ is below %s, B if above %s, ambiguous in between, unresolved under %s genes "
        "(DR02). A <b>block</b> is a run of measured copies along a chromosome with one region and one label "
        "(DR02c); <b>propagation</b> gives every gene inside a block's boundaries that label, including genes "
        "without a Δ." % (ff(SEG_LO, 2), ff(SEG_HI, 2), fi(MINSEG)))


def pipeline_section():
    R.h2("2. Pipeline at a glance")
    rows = []
    ls = D.get("locsum")
    if "loci" in S:
        extra = ""
        if ls is not None and "n_dionaea" in ls:
            nep = ls["nep"] if "nep" in ls else ls.get("has_nep")
            extra = "; %s with a Nepenthes copy, %s with exactly two Dionaea copies" % (
                fi(nep.astype(str).str.upper().isin(["TRUE", "T"]).sum()) if nep is not None else "–",
                fi((ls.n_dionaea == 2).sum()))
        rows.append(["1 locus set", "DR00_prep.R", "every Nepenthes-anchored syntenic orthogroup plus every Dionaea "
                     "homeolog pair, all copies", "%s loci, %s copies%s" % (fi(S["loci"]), fi(S["tips"]), extra)])
    if "ks_rows" in S:
        rows.append(["2 alignments, dS", "DR00_build.sh, ks/run_chunk.py (yn00)", "codon alignment per locus; dS "
                     "between every two copies; kept if 0 ≤ dS < %s and ≥ 100 codons (DR02)" % fi(DSMAX),
                     "%s pairs at %s loci; %s kept" % (fi(S["ks_rows"]), fi(S["ks_loci"]), fi(S["ks_kept"]))])
    if "axis" in S:
        rows.append(["3 Dionaea axis", "DR02_label.R §1", "loci where Dionaea kept one copy on each homeolog; "
                     "dropped if d(Dio_A, Dio_B) < %s (gene conversion)" % ff(DCONV, 2),
                     "%s with one copy each side, %s measured, %s dropped, %s used" % (
                         fi(S["axis_one_each"]), fi(S["axis_measured"]), fi(S["axis_conv"]), fi(S["axis"]))])
    if "measured_total" in S:
        rows.append(["4 Δ per copy", "DR02 §2", "Δ for every Drosera copy at an axis locus, minus the Nepenthes "
                     "offset", "%s genes" % fi(S["measured_total"])])
    if "tracks" in S:
        sc = S["seg_counts"].sum()
        rows.append(["5 segments", "DR02 §6", "each track split only where a changepoint beats %s shuffles; "
                     "median Δ gives the call" % fi(NPERM), "%s tracks → %s segments: %s" % (
                         fi(S["tracks"]), fi(S["segments"]), ", ".join("%s %s" % (fi(v), k) for k, v in sc.items()))])
    if "labelled_total" in S:
        rows.append(["6 blocks, propagation", "DR02c_blocks.R", "runs of measured copies with one region and label "
                     "(≥ %s pure, ≥ %s copies); every gene inside takes the label" % (ff(PUR_MIN, 2), fi(NVOTE_MIN)),
                     "%s genes labelled (%s%% of Drosera genes)" % (fi(S["labelled_total"]),
                                                                    ff(S["pct_labelled_all"]))])
    if "cn" in D:
        rows.append(["7 copies per locus", "DR21, DR22", "A and B copies per locus and species",
                     "%s locus × species" % fi(len(D["cn"]))])
    ngt = len(glob.glob("DR/tree/genetrees/*.treefile"))
    rows.append(["8 trees", "DR05, DR09, DR10", "concatenated ML (IQ-TREE), ASTRAL, ASTRAL-Pro on labelled copies",
                 "%s gene trees" % fi(ngt) if ngt else "gene-tree folder not found"])
    if "tracts_n" in S:
        rows.append(["9 riparian", "DRrip6_riparian_FINAL.R", "GENESPACE riparian, blocks coloured by label",
                     "%s tracts; %s braids dropped" % (fi(S["tracts_n"]), fi(S.get("dropped_n")))])
    if S.get("runs"):
        UD = {0: "prior only", 1: "exact likelihood", 2: "approximate likelihood", 3: "preparing the approximation"}
        ud = sorted({r["usedata"] for r in S["runs"]})
        ncal = sorted({len(r["calibs"]) for r in S["runs"] if r.get("tips") is not None})
        rows.append(["10 dating", "MCMCtree (DR/dating)", "Bayesian divergence times: %s; %s calibrated node%s per "
                     "tree" % (", ".join("usedata %d = %s" % (u, UD.get(u, "?")) for u in ud),
                               " or ".join(str(n) for n in ncal) or "–", "" if ncal == [1] else "s"),
                     "%d runs, %d with samples" % (len(S["runs"]), sum(1 for r in S["runs"] if r["mcmc"]))])
    R.table(["step", "code", "what it does", "what comes out"], rows)


def labels_section():
    R.h2("3. How the labels are made")
    R.p("One gene's Δ is noisy, so the call is made per segment: the median Δ of many neighbouring copies. Top row: "
        "per-gene Δ (dashed: the centres DR02 predicts for A and B, %s and %s). Bottom row: segment medians, coloured "
        "by the call." % (ff(PRED_A, 2), ff(PRED_B, 2)))
    show("delta", "the shape of the bottom row, not its colours (colours follow the thresholds %s and %s by "
         "construction). Two humps, one on each side of 0, mean segments separate into two kinds; one hump straddling "
         "0 would mean the thresholds cut through noise." % (ff(SEG_LO, 2), ff(SEG_HI, 2)))
    R.h3("How noisy is one gene's Δ?")
    if "nep" in S:
        chk = ""
        if "off_file" in S:
            same = abs(S["off_file"] - S["nep_median"]) < 1e-9
            chk = (" Re-derived here exactly as DR02 does, the Nepenthes median is %+.4f; DR02 subtracted %+.4f from "
                   "every Drosera Δ (%s)." % (S["nep_median"], S["off_file"], "identical" if same else
                                               "DIFFERENT: the re-derivation does not reproduce DR02"))
        rho = S.get("unit_rho", float("nan"))
        R.p("Two independent measures. Nepenthes sits outside the A/B split, so its true Δ is 0: across %s loci its "
            "Δ has median %+.3f (Wilcoxon p = %s: %s) and sd %s (robust sd %s). Inside one species × region × "
            "chromosome unit with at least %d measured genes, Δ spreads with a pooled sd of %s (median %s over %s "
            "units). If units straddling an A/B boundary inflated this, larger units would spread more; size and "
            "spread have Spearman correlation %s%s. DR02_label.R's design assumed a per-gene sd of about %s.%s" % (
                fi(len(S["nep"])), S["nep_median"], "%.2g" % S["nep_p"],
                "no detectable offset" if S["nep_p"] >= 0.05 else "an offset, which DR02 subtracts",
                ff(S["nep_sd"], 2), ff(S["nep_rsd"], 2), MIN_UNIT, ff(S.get("unit_sd"), 2),
                ff(S.get("unit_sd_median"), 2), fi(len(S.get("units", []))), ff(rho, 2),
                "" if not ok(rho) else (", so they do not" if rho < 0.2 else ", so part of the spread may be "
                                         "boundaries rather than noise"), ff(NOISE_PRED, 2), chk))
    show("noise", "whether the grey histogram is centred on 0 (no offset), and where the dots sit against the dashed "
         "line, the noise DR02 assumed: above it, single genes are noisier than the design assumed.")
    if "sep_sd_gene" in S:
        R.p("What this means for a segment. DR02 predicts the A and B centres %s apart; divided by the per-gene sd "
            "that is %s sd per gene as measured (%s in DR02's design). Averaging n genes multiplies the separation "
            "by √n, so at the median called segment (%s measured genes) the centres are %s sd apart (%s in the "
            "design). DR02's own per-segment score, conf = |median Δ| / (MAD/√n), the distance of the median from 0 "
            "in standard errors, is at least 3 for %s%% of called segments; %s%% when the larger standard error of "
            "a median (1.25 × MAD/√n) is used." % (
                ff(SEP, 2), ff(S["sep_sd_gene"], 2), ff(S["sep_sd_gene_design"], 2), fi(S.get("seg_n_median")),
                ff(S.get("sep_at_median"), 1), ff(S.get("sep_at_median_design"), 1),
                ff(100 * S.get("conf3", float("nan")), 0), ff(100 * S.get("conf3_strict", float("nan")), 0)))
    R.h3("How much is measured, how much carried")
    show("coverage", "the dark part of each bar is the evidence; the grey part is inherited from blocks. On the "
         "right, curves rising early mean most carried labels sit close to a measured gene.")
    cov = S.get("cov")
    if cov is not None:
        rows = []
        for g in DROS:
            if g not in cov.index:
                continue
            sc = S["seg_counts"].loc[g] if g in S.get("seg_counts", pd.DataFrame()).index else {}
            rows.append([isp(g), fi(cov.total[g]), fi(cov.measured[g]), ff(cov.pct_measured[g]), fi(cov.labelled[g]),
                         ff(cov.pct_labelled[g]), " / ".join(fi(sc.get(l, 0)) for l in ("A", "B")),
                         ff(100 * S.get("near", {}).get(g, float("nan"))),
                         ff(100 * S.get("far", {}).get(g, float("nan")))])
        R.table(["species", "genes", "Δ measured", "%", "labelled", "%", "segments A / B",
                 "carried labels within %s kb of a measured gene, %%" % fi(1000 * NEAR_MB),
                 "farther than %s Mb, %%" % ff(FAR_MB, 0)], rows)
    if "purity_median" in S:
        R.p("Vote purity (DR14): within a called segment, the share of measured genes whose own Δ agrees in sign "
            "with the call, genes with |Δ| &lt; %s left out. %sMedian purity %s over %s called segments; %s are below "
            "DR14's marker of %s. %s" % (
                ff(GENE_AMBIG, 2),
                ("With centres at ±%s and a per-gene sd of %s, a correctly called segment is expected to reach only "
                 "%s, so low purity is what per-gene noise produces, not by itself a wrong call. "
                 % (ff(abs(PRED_A), 2), ff(S["unit_sd"], 2), ff(S["purity_expected"], 2)))
                if "purity_expected" in S else "", ff(S["purity_median"], 2), fi(S["purity_n"]),
                fi(S["purity_marginal"]), ff(PUR_MARG, 1),
                ("DR14's binomial test finds the majority significant (p &lt; 0.05) in %s of the %s segments with at "
                 "least 3 votes." % (fi(S["purity_sig"]), fi(S["purity_tested"]))) if "purity_sig" in S else ""))


def validation_section():
    R.h2("4. Do the labels hold up?")
    R.p("Checks that do not reuse the quantity they test. Two earlier checks are deliberately absent because they "
        "are circular: run lengths of propagated labels (every gene in a block inherits the block's call, so long "
        "runs are guaranteed) and scoring windows of genes against a call made from those same genes.")
    R.h3("a. Is Δ correlated along chromosomes?")
    a = S.get("ac")
    if a:
        def reach(b):
            return ("stays above them up to lag %d, a median of %s Mb apart" % (
                b["reach"], ff(a["dmb"][b["reach"] - 1], 1)) if b["reach"] else "is not above them even at lag 1")
        bs, bc = a["band"]["species"], a["band"]["chromosome"]
        R.p("If ancestry comes in blocks, two measured genes close together on a chromosome should have similar Δ. "
            "Two shuffles say how similar chance makes them. Shuffling Δ among all measured genes of a species "
            "removes every positional signal: the observed lag-1 correlation, %s, is against a shuffle range of %s "
            "to %s, and %s. Shuffling only within each chromosome keeps each chromosome's average, so it asks for "
            "blocks shorter than a chromosome: range %s to %s, and the observed correlation %s. This uses no labels "
            "and no segmentation." % (ff(a["obs"][0], 3), ff(bs["lo"][0], 3), ff(bs["hi"][0], 3), reach(bs),
                                      ff(bc["lo"][0], 3), ff(bc["hi"][0], 3), reach(bc)))
    show("autocorr", "the black line above the light band at short lags (positional signal at all) and above the "
         "darker band (blocks shorter than a chromosome), falling into the bands at long lags. A line inside the "
         "light band would mean Δ carries no positional signal.")
    R.h3("b. Do A copies of different species group together?")
    xs = S.get("xs")
    if xs:
        R.p("Drosera-to-Drosera dS was never used to build Δ, which only measures distances to Dionaea. At %s loci "
            "where two species each carry an A- and a B-labelled copy, the test asks whether the first species' A "
            "copy is closer to the other species' A copy than to its B copy (and the same for B). Chance is 50%%. A: "
            "%s%% of %s tests (95%% CI %s–%s, p = %s); B: %s%% of %s (%s–%s, p = %s); with labels permuted within each "
            "species and locus, %s%% and %s%% (%d repeats). %s" % (
                fi(xs["loci"]), ff(100 * xs["A"]["f"]), fi(xs["A"]["n"]), ff(100 * xs["A"]["ci"][0]),
                ff(100 * xs["A"]["ci"][1]), "%.2g" % xs["A"]["p"], ff(100 * xs["B"]["f"]), fi(xs["B"]["n"]),
                ff(100 * xs["B"]["ci"][0]), ff(100 * xs["B"]["ci"][1]), "%.2g" % xs["B"]["p"],
                ff(100 * xs["ctl"]["A"][0]), ff(100 * xs["ctl"]["B"][0]), N_CTRL,
                ("%s copies group more consistently than %s copies (two-proportion test p = %.2g)." % (
                    (("A", "B") if xs["A"]["f"] > xs["B"]["f"] else ("B", "A")) + (xs["AvsB_p"],)))
                if xs["AvsB_p"] < 0.05 else "A and B do not differ (p = %.2g)." % xs["AvsB_p"]))
    show("cross", "green and orange dots right of the dashed line for every pair, grey crosses on it. Pairs whose "
         "dots sit near 50% are where the labels are weakest.")
    R.h3("c. Does Δ agree with the four-point condition?")
    if S.get("mcn"):
        R.p("The four-point condition on (gene_Drosera, Dio_A, Dio_B, Nepenthes) is immune to a rate difference "
            "between the Dio_A and Dio_B lineages; Δ is not. Both are noisy per gene, so the test is symmetry, not "
            "agreement: under noise alone the two kinds of disagreement are equally common (McNemar). %s%% of copies "
            "fall outside the A/B split by the four-point call and are left out." % ff(100 * S.get("outside", np.nan)))
        R.table(["species", "copies", "agreement, %", "Δ says A, four-point B", "Δ says B, four-point A",
                 "McNemar p"], S["mcn"])
    R.h3("d. The project's own assertions")
    t = S.get("tests")
    if t:
        R.p("DR/DRtest_phasing.R checks that the propagated product matches the segment calls (it was written after "
            "four bugs of the same kind). Run during this build with %s: %d passed, %d failed (exit status %d)." % (
                t["rscript"], len(t["passed"]), len(t["failed"]), t["rc"]))
        R.table(["assertion", "result", "detail when failing"], t["rows"])
        if t["failed"]:
            R.p("A failing assertion that pins a number (such as a count of regions) may mean the data changed since "
                "the assertion was written rather than a bug; either way it must be resolved before the labels are "
                "used.")
    elif A.no_tests:
        R.p("Not run in this build (--no_tests).")


def results_section():
    R.h2("5. Constitution, trees, riparian")
    R.h3("How many A and B copies does each species carry?")
    if S.get("cn_q"):
        R.p("Fractionation only removes copies, so across thousands of loci the upper tail of the per-locus count "
            "recovers the original number of copies. Predicted ceilings: %s. Squares mark the prediction; the dashed "
            "line is A = 2 × B. Bubble area is the number of loci (the largest bubble is labelled)." % pred_text())
    show("constitution", "where the outer envelope of bubbles sits: reaching the square and not beyond it supports "
         "the prediction; bubbles past it mean more copies than the constitution allows (or mislabelled copies). "
         "Most loci sit near the origin because fractionation removed copies.")
    if S.get("cn_q"):
        R.table(["species", "loci", "A p95 / p99 / max", "B p95 / p99 / max", "predicted A : B",
                 "A p95 if every unlabelled copy were A", "loci with fewer A than predicted (of these, % with no "
                 "unlabelled copy)"],
                [[isp(g), fi(q["loci"]), "%s / %s / %s" % (ff(q["A95"], 0), ff(q["A99"], 0), fi(q["Amax"])),
                  "%s / %s / %s" % (ff(q["B95"], 0), ff(q["B99"], 0), fi(q["Bmax"])),
                  "%s : %s" % (fi(q["predA"]), fi(q["predB"])), ff(q.get("A95_worst"), 0),
                  "%s (%s%%)" % (fi(q.get("short")), ff(100 * q.get("short_explained", float("nan"))))]
                 for g, q in S["cn_q"].items()])
        R.p("The last two columns are DR22's truncation check. If a missing A copy were present but unlabelled, "
            "counting every unlabelled copy as A would raise the ceiling; at a locus with no unlabelled copy there "
            "is nothing that could be the missing one, so its shortfall is a real loss.")
    R.h3("Species trees")
    R.p("ASTRAL-Pro on multi-copy gene trees (DR09; figure DR10), branches annotated with q1, the share of "
        "gene-tree quartets supporting the branch (1/3 = no signal; red below %s), and branch length in coalescent "
        "units (short = much incomplete lineage sorting)." % ff(Q1_RED, 2))
    R.img(PIPE_FIGS["astralpro"], 80)
    R.look("whether A tips and B tips form two clades, and where regia sits in each.")
    R.p("The 13-tip concatenated tree (DR05h), with real branch lengths (substitutions per site); nodes with low site "
        "concordance are drawn in red by that script.")
    R.img(PIPE_FIGS["concat13"], 75)
    if "c05f" in S:
        R.h3("Does regia group with Dionaea, or is that the labelling?")
        txt = ("DR05f takes each species in turn as the focal species. At every locus where the focal species and "
               "another Drosera species carry a copy with the same label (A or B), the four-point condition on (focal "
               "copy, Dionaea's copy with that label, the other species' copy, Nepenthes) picks which two of the four "
               "are closest relatives, and each locus votes by majority. Labels were made by closeness to Dionaea, so "
               "every species is pulled towards Dionaea a little; that pull is measured by the species known to sit "
               "inside Drosera. A real regia–Dionaea affinity raises regia alone.")
        if "c05f_gap" in S:
            r = S["c05f_rg"]
            lo, hi, top = S["c05f_oth"]
            txt += (" regia pairs with Dionaea at %s%% of %s loci (95%% CI %s–%s%%); the other species at %s–%s%% "
                    "(highest: %s). The gap, %+.1f percentage points, gives DR05f's verdict by its own rule (gap > %s: "
                    "clearly highest; ≤ %s: not distinguishable; between: suggestive): <b>%s</b>." % (
                        ff(100 * r.focal_Dionaea), fi(r.loci), ff(100 * r.lo), ff(100 * r.hi), ff(100 * lo),
                        ff(100 * hi), top, 100 * S["c05f_gap"], ff(F05_CLEAR, 2), ff(F05_NONE, 2), S["c05f_verdict"]))
        R.p(txt)
        show("regia", "regia's bar against the other species' bars (their height is the labelling pull), and whether it "
             "clears the dashed line.", 70)
    R.h3("Riparian")
    R.img(PIPE_FIGS["riparian"], 100)
    R.look("green (A) and orange (B) strips under each chromosome, and braids that land on their own colour.")
    if "tracts_n" in S:
        R.p("%s tracts of at least %s genes are drawn (%s), and %s chromosomes carry both labels. A braid is coloured "
            "only if ≥ %s%% of its block's genes share one label, and %s braids were dropped because an end lands on "
            "the other subgenome (DR/out/DR25_dropped.csv): possibly genuine chimeric blocks, not checked." % (
                fi(S["tracts_n"]), fi(TRACT_MIN), ", ".join("%s %s" % (fi(v), k) for k, v in S["tracts_lab"].items()),
                fi(S["tracts_both"]), ff(100 * BRAID_PUR, 0), fi(S.get("dropped_n"))))


def dating_section():
    R.h2("6. Dating")
    runs = S.get("runs")
    if not runs:
        return
    def cal_ma(c):
        if not c:
            return ""
        lo, hi = c
        return ("%s–%s" % (ff(lo * TIME_UNIT_MA), ff(hi * TIME_UNIT_MA)) if ok(lo) and ok(hi) else
                "over %s" % ff(lo * TIME_UNIT_MA) if ok(lo) else "under %s" % ff(hi * TIME_UNIT_MA))
    ex = next((r["ctl"]["RootAge"] for r in runs if calib(r["ctl"].get("RootAge", ""))), None)
    R.p("MCMCtree runs in DR/dating, one folder each. MCMCtree's time unit here is %g Myr%s. A tree is rooted "
        "correctly if Nepenthes is one of the root's two children. ESS (effective sample size, the number of "
        "independent draws a chain is worth) is computed as coda::effectiveSize does; a replicate group's pooled ESS "
        "is the sum over its chains. Ages are cited only for a correctly rooted tree whose pooled ESS is ≥ %d on "
        "every node; for other runs only their convergence is shown." % (
            TIME_UNIT_MA, ", so RootAge = %s means %s Ma" % (html_text(ex.strip("'\"")), cal_ma(calib(ex))) if ex
            else "", ESS_MIN))
    keys = ["usedata", "clock", "RootAge", "nsample", "sampfreq", "burnin", "finetune", "BDparas", "rgene_gamma",
            "sigma2_gamma", "alpha", "model", "seqfile", "treefile"]
    keys = [k for k in keys if any(k in r["ctl"] for r in runs)]
    same = [k for k in keys if len({r["ctl"].get(k, "") for r in runs}) == 1]
    vary = [k for k in keys if k not in same]
    if same:
        R.p("Settings shared by all %d runs: %s." % (len(runs), "; ".join(
            "%s = %s%s" % (k, html_text(runs[0]["ctl"].get(k, "")),
                           " (%s Ma)" % cal_ma(calib(runs[0]["ctl"][k])) if k == "RootAge" and calib(
                               runs[0]["ctl"].get(k, "")) else "") for k in same)))
    rows = []
    for r in runs:
        c = r["ctl"]
        rows.append([r["name"]] + [c.get(k, "") for k in vary] + [
            "yes" if r["rooted_ok"] else ("NO" if r["rooted_ok"] is not None else "tree not found"),
            "YES" if r["double_root"] else "no", fi(len(r["calibs"])) if r.get("tips") is not None else "–",
            ("yes" if r["inBV"] else "MISSING") if r["usedata"] == 2 else "–",
            fi(r.get("samples", 0)) if r["mcmc"] else "none", fi(r.get("ess_min_nodes")) if r["mcmc"] else "–"])
    R.table(["run"] + vary + ["rooted on Nepenthes", "root calibrated in ctl and tree", "calibrated nodes",
                              "in.BV (needed for usedata 2)", "samples", "min ESS over nodes"], rows)
    show("ess", "bars past the dashed line on the left; dots near 0 on the right. Replicate chains that disagree "
         "on a node have not settled, whatever their ESS.")
    best = S.get("best", {})
    if best:
        R.p("Best replicate group per tree: %s." % "; ".join(
            "%s → %s, pooled ESS %s (%s)" % (f, g["names"], fi(g["pooled_min"]),
                                            "citable" if g["citable"] else "not citable: " + g["why_not"])
            for f, g in ((f, best[f]) for f in chrono_families())))
    show("chrono", "bars (95% HPD intervals) inside the grey calibration boxes on calibrated nodes; uncalibrated "
         "nodes show what the sequences say.", absent="No tree has converged, correctly rooted chains yet, so no "
         "chronogram is drawn.")
    for f in chrono_families():
        g = best[f]
        what = fam_label(f)
        if g["citable"]:
            R.table(["%s: node (t_n)" % what, "posterior mean, Ma", "95% HPD, Ma", "pooled ESS", "calibration, Ma"],
                    [["%s (t_n%d)" % (s["name"], s["num"]), ff(s["mean"]), "%s–%s" % (ff(s["lo"]), ff(s["hi"])),
                      fi(s["ess"]), cal_ma(s["cal"])] for s in g["summary"]])
        elif g["rooted_ok"]:            # converging but not there yet: which nodes lag is worth seeing
            R.table(["%s: node (t_n)" % what, "ages", "pooled ESS (%s)" % g["names"], "calibration, Ma"],
                    [["%s (t_n%d)" % (s["name"], s["num"]), "withheld: " + g["why_not"], fi(s["ess"]),
                      cal_ma(s["cal"])] for s in g["summary"]])
    show("prior", "the posterior (colour) narrower than or shifted from the prior (grey, the same run without "
         "sequence data). Where they coincide, the age comes from the calibration, not from the sequences.",
         absent="Prior-only runs (usedata = 0) are compared with the posterior only for citable trees; none yet."
         if not any(g["citable"] for g in best.values()) else None)
    noch = [r for r in runs if not r["mcmc"] and r["usedata"] != 3]
    if noch:
        R.p("Runs without samples: %s. Logs found in their folders or in DR/dating are in the appendix." % ", ".join(
            r["name"] for r in noch))
    mis = [r["name"] for r in runs if r["rooted_ok"] is False]
    if mis:
        note = ""
        if S.get("dr24"):
            src = S.get("dr24_src") or []
            used = sorted({x for _, rr in src for x in rr} & set(mis))
            note = (". %s exists; %s" % (PIPE_FIGS["dr24"], (
                "it is drawn by %s from %s, so it must not be used either" % (
                    ", ".join(f for f, _ in src), ", ".join(used)) if used else
                "it is drawn by %s, which does not name these runs" % ", ".join(f for f, _ in src) if src else
                "no script under DR/ names it, so check which runs it was drawn from")))
        R.p("<b>Mis-rooted</b> (Nepenthes not a child of the root): %s. Ages from these runs must not be used%s." % (
            ", ".join(mis), note))


def limitations_section():
    R.h2("7. Limitations")
    li = []
    if "pct_measured_all" in S:
        li.append("Only %s%% of Drosera genes carry a measured Δ: Δ needs both Dionaea homeologs at the locus, and "
                  "fractionation removed one at most loci. Trees and the riparian run on labels carried from blocks "
                  "(%s%% of genes)." % (ff(S["pct_measured_all"]), ff(S["pct_labelled_all"])))
    if "unit_sd" in S:
        li.append("Δ's per-gene noise is %s, against %s assumed in DR02_label.R: single-gene labels are not reliable; "
                  "only segment-level calls are." % (ff(S["unit_sd"], 2), ff(NOISE_PRED, 2)))
    li.append("A is defined as the more-retained Dionaea homeolog of each chromosome pair. Fractionation makes each "
              "pair asymmetric, but nothing in it links the dominant member of one pair to the dominant member of "
              "another: phasing across pairs assumes genome-wide subgenome dominance in Dionaea.")
    bad = S.get("q06_bad")
    if bad is not None and len(bad):
        li.append("DR06 quartet counts: in %d of %d rows the three resolutions do not add up to n (%s). DR06 names a "
                  "split by whichever of its two pairs the unrooted gene tree happens to hold as a clade, so one split "
                  "can be counted under two names and the third-ranked name drops out. No ILS-versus-introgression "
                  "symmetry claim can rest on DR06 until it is fixed." % (
                      len(bad), len(D["q06"]), ", ".join(bad.test.astype(str).head(4))))
    if S.get("modelfit"):
        li.append("DR10_SUPP_modelfit.pdf is not an independent check: ASTRAL computes coalescent-unit branch lengths "
                  "from the same quartet frequencies, so predicted and observed concordance agree by construction. "
                  "Its subtitle says they are independent.")
    if S.get("runs"):
        for x in stalled():
            li.append("The %s runs (%s-tip tree) have produced no samples%s. Root calibrated in both ctl and tree "
                      "file: %s." % (x["fam"], fi(x["ntip"]), ", so the age of the A/B split is unknown" if x["has_AB"]
                                     else "", ", ".join(x["double"]) or "in none of these runs"))
        for f in chrono_families():
            g = S["best"][f]
            if not g["citable"]:
                li.append("%s: no citable chains (%s); its ages are not reported." % (
                    fam_label(f)[0].upper() + fam_label(f)[1:], g["why_not"]))
    if S.get("untracked"):
        scripts = [os.path.basename(f) for f in S["untracked"] if not f.endswith(".ctl")]
        li.append("%d code or config files under DR/ are not under version control (listed in the appendix)%s." % (
            len(S["untracked"]), ", among them " + ", ".join(scripts[:6]) if scripts else ""))
    if "dropped_n" in S:
        li.append("%s riparian braids were dropped as possible chimeric blocks; whether they are real was not "
                  "checked." % fi(S["dropped_n"]))
    t = S.get("tests")
    if t and t["failed"]:
        li.append("DRtest_phasing.R fails %d of %d assertions (%s); see 4d." % (
            len(t["failed"]), len(t["rows"]), "; ".join(t["failed"])))
    R.ul(li)


def next_section():
    R.h2("8. Next steps and asks")
    li = []
    if S.get("runs"):
        for x in stalled():
            ref = x["ref"]
            li.append("Find out why the %s runs stop: their last log lines are in the appendix%s%s%s." % (
                x["fam"],
                "; in.BV is missing in %s" % ", ".join(x["no_inbv"]) if x["no_inbv"] else "",
                "; they differ from %s, which produced samples, in %s" % (ref["name"], "; ".join(x["diff"]))
                if ref and x["diff"] else "",
                "; and the root is calibrated in both the ctl (RootAge) and the tree file in %s, the cause suspected "
                "earlier" % ", ".join(x["double"]) if x["double"] else ""))
        for f in ("A", "B"):
            g = S.get("best", {}).get(f)
            if g and not g["citable"] and g["rooted_ok"]:
                li.append("Run longer subgenome %s chains (or more replicates) until pooled ESS ≥ %d on every node; "
                          "the worst node now has %s." % (f, ESS_MIN, fi(g["pooled_min"])))
    t = S.get("tests")
    if t and t["failed"]:
        li.append("Resolve the failing DRtest_phasing.R assertion%s (bug, or a pinned number to update)." % (
            "s" if len(t["failed"]) > 1 else ""))
    if S.get("q06_bad") is not None and len(S["q06_bad"]):
        li.append("Fix DR06's quartet naming (name each split canonically, e.g. by the pair containing the first tip) "
                  "before any symmetry test.")
    if S.get("untracked"):
        li.append("Commit the untracked scripts and run configs.")
    li.append("Ask: agree whether the A/B assignment across Dionaea chromosome pairs (fractionation dominance) is "
              "acceptable as stated, or needs an independent test before publication.")
    R.ul(li)


def appendix():
    R.h2("Appendix")
    t = S.get("tests")
    if t:
        R.details("DRtest_phasing.R output (exit %d)" % t["rc"], text=t["out"])
    runs = S.get("runs", [])
    if runs:
        # one control file in full, then per run only the settings that differ from it
        ref = runs[0]
        R.details("mcmctree.ctl of %s (reference for the comparison below)" % ref["name"], ref["ctlp"])
        diff = []
        for r in runs[1:]:
            d = ["%s = %s (%s: %s)" % (k, r["ctl"].get(k, "(absent)"), ref["name"], ref["ctl"].get(k, "(absent)"))
                 for k in sorted(set(r["ctl"]) | set(ref["ctl"]), key=lambda k: list(ref["ctl"]).index(k)
                                 if k in ref["ctl"] else 99) if r["ctl"].get(k) != ref["ctl"].get(k)]
            diff.append("%s: %s" % (r["name"], "; ".join(d) if d else "identical"))
        R.details("mcmctree.ctl of every other run: settings that differ from %s" % ref["name"], text="\n".join(diff))
        # node numbering, once per distinct tree file
        trees = {}
        for r in runs:
            if r.get("tips") is not None:
                h = hashlib.md5(r["treetext"].encode()).hexdigest()
                trees.setdefault(h, []).append(r)
        for h, rs in trees.items():
            r = rs[0]
            txt = "tips: %s\n" % ", ".join("%d %s" % (i + 1, x) for i, x in enumerate(r["tips"]))
            txt += "\n".join("t_n%d  %s%s" % (nd["num"], nd["name"], "  " + nd["label"] if nd["label"] else "")
                             for nd in r["nodes"])
            txt += "\n\n" + r["treetext"].strip()
            R.details("node numbers of the tree file used by %s" % ", ".join(x["name"] for x in rs),
                      text=txt)
        for r in runs:
            if not r["mcmc"]:
                for lg in r["logs"][:3]:
                    try:
                        tail = "".join(open(lg, errors="replace").readlines()[-25:])
                    except OSError:
                        continue
                    R.have(lg)
                    R.details("%s: last lines of %s" % (r["name"], os.path.basename(lg)), text=tail)
        for lg in S.get("dating_logs", [])[:6]:
            try:
                tail = "".join(open(lg, errors="replace").readlines()[-25:])
            except OSError:
                continue
            R.have(lg)
            R.details("last lines of %s" % lg, text=tail)
    for sb in sorted(glob.glob("DR/dating/*.sbatch")):
        R.details("job script %s" % sb, sb)
    if S.get("untracked") is not None:
        R.details("code and configs under DR/ not under version control (%d)" % len(S["untracked"]),
                  text="\n".join(S["untracked"]) or "none")
    if S.get("modified"):
        R.details("tracked files with uncommitted changes", text="\n".join(S["modified"]))


for name, fn in (("summary", summary), ("material", material_section), ("pipeline", pipeline_section),
                 ("labels", labels_section), ("validation", validation_section), ("results", results_section),
                 ("dating", dating_section), ("limitations", limitations_section), ("next steps", next_section),
                 ("appendix", appendix)):
    R.guarded("section " + name, fn)
R.write()
