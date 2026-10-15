#!/usr/bin/env python3
"""ETA for a running CARP run, from its own log and from a finished run of another species.

usage: carp_eta.py RUNNING_LOG RUNNING_BENCHMARK_DIR REFERENCE_BENCHMARK_DIR [RUN_NAME REF_NAME]
Each remaining step is estimated as the reference species' time for that step x a scale factor,
the median of (running / reference) over the steps both have finished (> 5 min in the reference).
Steps are summed (CARP mostly runs them one after another), so the ETA leans late.
"""
import datetime as dt
import glob
import os
import re
import statistics
import sys

log, bm_run, bm_ref = sys.argv[1:4]
RUN, REF = (sys.argv[4:6] + ["this run", "reference"][len(sys.argv[4:6]):])[:2]


def bench(d):
    out = {}
    for f in glob.glob(os.path.join(d, "*.tsv")):
        rows = [l.split("\t") for l in open(f) if l.strip()]
        if len(rows) >= 2:
            try:
                out[os.path.basename(f)[:-4]] = float(rows[-1][0])
            except ValueError:
                pass
    return out


ref, run = bench(bm_ref), bench(bm_run)
ts_re = re.compile(r"^\[(\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})\]$")
lines = open(log, errors="replace").read().split("\n")
# only this attempt: from the last CARP header on
start = max([i for i, l in enumerate(lines) if l.startswith("=== CARP ") and " for " in l] or [0])
lines = lines[start:]
ts, cur, jobrule, started, finished, failed, stats = None, None, {}, {}, {}, set(), []
in_stats = False
for l in lines:
    m = ts_re.match(l.strip())
    if m:
        ts = dt.datetime.strptime(re.sub(" +", " ", m.group(1)), "%a %b %d %H:%M:%S %Y")
        continue
    if l.startswith("Job stats:") and not stats:
        in_stats = True
        continue
    if in_stats:
        p = l.split()
        if len(p) == 2 and p[1].isdigit() and p[0] != "total":
            stats.append(p[0])
        elif p and p[0] == "total":
            in_stats = False
        continue
    m = re.match(r"^(?:local)?rule (\w+):$", l)
    if m:
        cur = m.group(1)
        continue
    m = re.match(r"^\s+jobid: (\d+)$", l)
    if m and cur:
        jobrule[m.group(1)] = cur
        started.setdefault(cur, ts)
        cur = None
        continue
    m = re.match(r"^Finished job(?:id:)? (\d+)", l)
    if m and m.group(1) in jobrule:
        finished[jobrule[m.group(1)]] = ts
        continue
    m = re.match(r"^Error in rule (\w+):", l)
    if m:
        failed.add(m.group(1))

now = dt.datetime.now()
ratios = {r: run[r] / ref[r] for r in run if r in ref and ref[r] > 300}
vals = sorted(ratios.values())
scale = statistics.median(vals) if vals else 3.0
if len(vals) >= 4:
    lo, _, hi = statistics.quantiles(vals, n=4)      # middle half of the steps' ratios
elif vals:
    lo, hi = vals[0], vals[-1]
else:
    lo, hi = 2.0, 4.0
print("   scale, %s time / %s time, steps both finished (over 5 min for %s):" % (RUN, REF, REF))
for r, x in sorted(ratios.items(), key=lambda kv: -ref[kv[0]]):
    print("     %-44s %s %5.2f h  %s %5.2f h  x%.1f" % (r, REF, ref[r] / 3600, RUN, run[r] / 3600, x))
print("     median x%.1f; the range below uses the middle half of these, x%.1f to x%.1f" % (scale, lo, hi))

todo = [r for r in stats if r != "all" and r not in finished and r not in failed]
running = [r for r in todo if r in started]
waiting = [r for r in todo if r not in started]
rows, tot, tot_lo, tot_hi = [], 0.0, 0.0, 0.0
for r in running + waiting:
    base = ref.get(r, 120.0)
    est, elo, ehi = base * scale, base * lo, base * hi
    if r in started:
        el = (now - started[r]).total_seconds()
        left = max(est - el, 0.1 * est)
        left_lo, left_hi = max(elo - el, 0.1 * elo), max(ehi - el, 0.1 * ehi)
        note = "running %.1f h (%s took %.1f h)" % (el / 3600, REF, base / 3600)
    else:
        left, left_lo, left_hi = est, elo, ehi
        note = "%s took %.2f h" % (REF, base / 3600) if r in ref else "no %s time; 2 min assumed" % REF
    tot += left
    tot_lo += left_lo
    tot_hi += left_hi
    rows.append((left, r, note))
print("   still to do in this attempt (%d of %d steps; estimated %s hours, steps over 5 min):"
      % (len(todo), len([r for r in stats if r != 'all']), RUN))
for left, r, note in rows:
    if left >= 300:
        print("     %-44s ~%4.1f h   %s" % (r, left / 3600, note))
small = sum(1 for left, r, n in rows if left < 300)
print("     + %d shorter steps, %.1f h together" % (small, sum(left for left, r, n in rows if left < 300) / 3600))
eta, eta_lo, eta_hi = (now + dt.timedelta(seconds=x) for x in (tot, tot_lo, tot_hi))
f = "%a %H:%M"
print("   CARP done around %s (range %s to %s). Steps are summed, but CARP runs some side by side, so this leans late.\n"
          "   Publishing takes ~10 min after that, and the CO_smk plots follow within 10 min"
      % (eta.strftime(f), eta_lo.strftime(f), eta_hi.strftime(f)))
if failed:
    print("   failed so far in this attempt (the run continues): %s" % ", ".join(sorted(failed)))
