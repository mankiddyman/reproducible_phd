# ---------------------------------------------------------------------------
# Repeat and rDNA annotation with CARP (Comprehensive Annotation of Repeats
# Pipeline, https://github.com/kavonrtep/CARP), one run per species on the
# frozen two-haplotype {species}_chr.fa.
#
# Why one run on both haplotypes: CARP builds its repeat library and numbers its
# tandem-repeat families (TRC_<n>) from the genome it is given. One run over
# both haplotypes gives one library and one set of family numbers, so TRC_12 on
# chr1_hap1 and TRC_12 on chr2_hap2 are the same satellite family.
#
# CARP is a container that runs its own Snakemake inside: DANTE (protein
# domains), DANTE_LTR / DANTE_TIR / DANTE_LINE (whole elements), TideCluster
# (tandem repeats and satellites, including whole 45S and 5S rDNA arrays) and
# RepeatMasker with the library it builds (TEs, rDNA subunits 18S/5.8S/25S/IGS,
# simple repeats). Its authoritative output is Repeat_Annotation_Unified.gff3:
# non-overlapping, one classification per base.
#
# Resume-safe: CARP keeps its intermediates when a run stops (--keep-incomplete)
# and the retried job continues with --rerun-incomplete. Only the final files
# are declared as outputs, so Snakemake never deletes CARP's work in progress.
#
# Side steps do not block the annotation: CARP's Snakemake runs with --keep-going
# (config carp_extra), and a run that exits non-zero is accepted when the two
# files published from it (Repeat_Annotation_Unified.gff3, summary_statistics.csv)
# are complete. The steps that failed (density tracks, the HTML report) are named
# in carp.done and in the published MANIFEST.
#
# The image is fetched once on the login node by scripts/fetch_carp_sif.sh
# (compute nodes need no internet). Consumers (CO_smk) read the published copy
# in refs/{species}_repeats/, never results/.
# ---------------------------------------------------------------------------

REPEATS_CFG = config.get("repeats", {})
CARP_VERSION = str(REPEATS_CFG.get("carp_version", "1.9.5"))
CARP_SIF = REPEATS_CFG.get(
    "carp_sif",
    "/netscratch/dep_mercier/grp_marques/Aaryan/methods/carp/carp_%s.sif" % CARP_VERSION)
REPEAT_SPECIES = list(REPEATS_CFG.get("species", ["Drosera_paradoxa", "Drosera_binata"]))


def _carp_res(species):
    """Slurm resources for one species' CARP job (config: repeats.resources)."""
    res = dict(REPEATS_CFG.get("resources", {}).get("default", {}))
    res.update(REPEATS_CFG.get("resources", {}).get(species, {}))
    res.setdefault("threads", 48)
    res.setdefault("mem_mb", 360000)
    res.setdefault("partition", "normal")
    res.setdefault("time", "2-00:00:00")
    return res


rule annotate_repeats_carp:
    """CARP on the frozen two-haplotype assembly: TEs, satellites, rDNA arrays."""
    input:
        genome=lambda wc: frozen_chr_fasta(wc.species),
        sif=CARP_SIF,
    output:
        gff="results/{species}/repeats/carp/Repeat_Annotation_Unified.gff3",
        stats="results/{species}/repeats/carp/summary_statistics.csv",
        done="results/{species}/repeats/carp_run/carp.done",
    params:
        outdir="results/{species}/repeats/carp",
        rundir="results/{species}/repeats/carp_run",
        rm_sensitivity=REPEATS_CFG.get("repeatmasker_sensitivity", "default"),
        rm_culling=int(REPEATS_CFG.get("repeatmasker_culling_limit", 2)),
        # concurrent workers in CARP's density-track (BigWig) steps; CARP's default 0 = one per thread
        bigwig_workers=int(REPEATS_CFG.get("bigwig_max_workers", 16)),
        tc_culling=int(REPEATS_CFG.get("tidecluster_reannotate_culling_limit", 0)),
        version=CARP_VERSION,
        # extra arguments for CARP's own Snakemake (config: repeats.resources.<species>.carp_extra),
        # e.g. "--set-threads dante=32" to cap one step's workers
        extra=lambda wc: str(_carp_res(wc.species).get("carp_extra", "")),
    threads: lambda wc: int(_carp_res(wc.species)["threads"])
    resources:
        cpus_per_task=lambda wc: int(_carp_res(wc.species)["threads"]),
        mem_mb=lambda wc: int(_carp_res(wc.species)["mem_mb"]),
        slurm_partition=lambda wc: _carp_res(wc.species)["partition"],
        # no inner quotes here: the plugin pastes this into the sbatch command line as it is
        # (the extra quotes in profiles/slurm/config.yaml are removed when Snakemake reads the profile)
        slurm_extra=lambda wc: "--time=%s --exclude=%s" % (_carp_res(wc.species)["time"],
                                                        _carp_res(wc.species).get("exclude", "hpc15")),
    wildcard_constraints:
        species="|".join(REPEAT_SPECIES),
    retries: 1
    log:
        "logs/annotate_repeats_carp/{species}.log",
    benchmark:
        "benchmarks/annotate_repeats_carp/{species}.tsv"
    shell:
        r"""
        REPO_ROOT=$(pwd)
        LOG="$REPO_ROOT/{log}"
        OUT="$REPO_ROOT/{params.outdir}"
        RUN="$REPO_ROOT/{params.rundir}"
        GENOME=$(readlink -f "$REPO_ROOT/{input.genome}")
        rm -rf "$RUN/tmp"
        mkdir -p "$OUT" "$RUN/tmp" "$(dirname "$LOG")"
        set +eu
        module load singularity-binds 2>/dev/null || module load singularity 2>/dev/null
        set -eu
        RT=$(command -v apptainer || command -v singularity || true)
        if [ -z "$RT" ]; then
            echo "no apptainer or singularity on $(hostname)" >> "$LOG"; exit 1
        fi
        # threads: what Slurm gave the job. Memory: CARP sizes its memory-hungry worker pools from -m.
        # These nodes use strict overcommit (vm.overcommit_memory=2): the kernel refuses new memory once
        # the node's committed total reaches CommitLimit, and processes outside Slurm can hold hundreds
        # of GB of it (698 GB on hpc04 on 2026-10-14, when DANTE's fork of 127 workers failed with ENOMEM).
        # So -m is the smallest of: 90% of the allocation, 80% of what the node can still commit,
        # 90% of what it has available.
        T=${{SLURM_CPUS_PER_TASK:-{threads}}}
        ALLOC_GB=$(( {resources.mem_mb} * 9 / 10 / 1000 ))
        read CL_KB CA_KB MA_KB < <(awk '/^CommitLimit:/ {{cl = $2}} /^Committed_AS:/ {{ca = $2}} /^MemAvailable:/ {{ma = $2}} END {{print cl, ca, ma}}' /proc/meminfo)
        MEM_GB=$ALLOC_GB
        ROOM_GB=$(( (CL_KB - CA_KB) / 1048576 * 8 / 10 ))
        AVAIL_GB=$(( MA_KB / 1048576 * 9 / 10 ))
        if [ "$(cat /proc/sys/vm/overcommit_memory)" = 2 ] && [ $ROOM_GB -lt $MEM_GB ]; then MEM_GB=$ROOM_GB; fi
        if [ $AVAIL_GB -lt $MEM_GB ]; then MEM_GB=$AVAIL_GB; fi
        if [ $MEM_GB -lt 16 ]; then MEM_GB=16; fi
        # open files: the soft limit is 1024 on these nodes, and CARP raises it for only a few of its steps
        ulimit -n "$(ulimit -Hn)" 2> /dev/null || true
        CFG="$RUN/carp_config.yaml"
        cat > "$CFG" <<EOF
genome_fasta: $GENOME
output_dir: $OUT
repeatmasker_sensitivity: {params.rm_sensitivity}
repeatmasker_culling_limit: {params.rm_culling}
tidecluster_reannotate_culling_limit: {params.tc_culling}
tidecluster_detect_rdna: true
max_memory_gb: $MEM_GB
bigwig_max_workers: {params.bigwig_workers}
cleanup_intermediates: minimal
EOF
        {{
          echo "=== CARP {params.version} for {wildcards.species} on $(hostname), $(date)"
          echo "    job ${{SLURM_JOB_ID:-local}}, $T threads, $MEM_GB GB, runtime $RT ($("$RT" --version 2>&1))"
          echo "    image {input.sif}"
          echo "    node memory: overcommit_memory=$(cat /proc/sys/vm/overcommit_memory) ratio=$(cat /proc/sys/vm/overcommit_ratio);" \
               "$(grep -E '^(MemTotal|MemAvailable|CommitLimit|Committed_AS)' /proc/meminfo | tr -s ' ' | tr '\n' ' ')"
          echo "    limits: virtual memory $(ulimit -v), processes $(ulimit -u), open files $(ulimit -n)"
          echo "    memory budget for CARP: $MEM_GB GB (90% of the allocation $ALLOC_GB GB; 80% of the node's" \
               "commit room $(( ROOM_GB * 10 / 8 )) GB; 90% of available $(( AVAIL_GB * 10 / 9 )) GB)"
          echo "    CARP's Snakemake gets: --rerun-incomplete --nolock {params.extra}"
          cat "$CFG"
        }} >> "$LOG"
        # Resuming: a TideCluster step that was interrupted leaves files that TideCluster refuses to
        # overwrite ('pass --force'), so the rerun fails at once. CARP's Snakemake marks the outputs of
        # an interrupted step as incomplete; when a TideCluster step is marked, clear its folder first.
        python3 - "$RUN/.snakemake/incomplete" "$OUT" >> "$LOG" 2>&1 <<'PYEOF'
import base64, os, shutil, sys
inc, out = sys.argv[1], sys.argv[2]
paths = []
for root, dirs, files in os.walk(inc):
    rel = os.path.relpath(root, inc)
    pre = "" if rel == "." else "".join(x.lstrip("@") for x in rel.split(os.sep))
    for f in files:
        b = pre + f
        try:
            paths.append(base64.urlsafe_b64decode(b + "=" * (-len(b) % 4)).decode())
        except Exception:
            pass
for sub in ("TideCluster/default", "TideCluster/short_monomer"):
    d = os.path.join(out, sub)
    if os.path.isdir(d) and any(p.startswith(d + "/") for p in paths):
        shutil.rmtree(d)
        print("    resume: removed the interrupted %s (CARP reruns that step from the start)" % d)
PYEOF
        # bind /netscratch (repo, image, outputs) and the genome's top directory if
        # it lives elsewhere; /tmp inside the container is a folder on netscratch
        # (node-local /scratch filled up during ANNEVO)
        BINDS="-B /netscratch"
        GTOP="/$(echo "$GENOME" | cut -d/ -f2)"
        [ "$GTOP" = /netscratch ] || BINDS="$BINDS -B $GTOP"
        # hpc09's Singularity runs without setuid, so every start unpacks the image
        # into a temporary sandbox (a few minutes): put that on netscratch, not in
        # the node's /tmp, and clear what a killed earlier attempt left behind
        export SINGULARITY_TMPDIR="$RUN/sif_tmp" APPTAINER_TMPDIR="$RUN/sif_tmp"
        rm -rf "$RUN/sif_tmp"
        mkdir -p "$RUN/sif_tmp"
        cd "$RUN"
        set +e
        "$RT" run --cleanenv $BINDS \
            -B "$RUN/tmp:/tmp" -B "$RUN/tmp:/var/tmp" \
            "{input.sif}" \
            -c "$CFG" -t $T -m $MEM_GB \
            -S="--rerun-incomplete --nolock {params.extra}" >> "$LOG" 2>&1
        rc=$?
        set -e
        echo "=== CARP exit $rc at $(date)" >> "$LOG"
        FAILED=none
        if [ $rc -ne 0 ]; then
            # CARP ran with --keep-going, so a failed side step did not stop the rest. Accept the run when
            # the two files published from it are complete: present, not empty, and not marked incomplete
            # by CARP's Snakemake. Otherwise fail (and Snakemake retries once; CARP resumes).
            python3 - "$RUN/.snakemake/incomplete" "$OUT/Repeat_Annotation_Unified.gff3" "$OUT/summary_statistics.csv" >> "$LOG" 2>&1 <<'PYEOF' || exit $rc
import base64, os, sys
inc, need = sys.argv[1], sys.argv[2:]
marked = set()
for root, dirs, files in os.walk(inc):
    rel = os.path.relpath(root, inc)
    pre = "" if rel == "." else "".join(x.lstrip("@") for x in rel.split(os.sep))
    for f in files:
        b = pre + f
        try:
            marked.add(base64.urlsafe_b64decode(b + "=" * (-len(b) % 4)).decode())
        except Exception:
            pass
bad = [p for p in need if not os.path.exists(p) or os.path.getsize(p) == 0 or p in marked]
for p in bad:
    print("    not complete: %s" % p)
sys.exit(1 if bad else 0)
PYEOF
            FAILED=$(grep -o -E "Error in rule [A-Za-z0-9_]+" "$LOG" | awk '{{print $4}}' | sort -u | paste -sd, - || true)
            FAILED=${{FAILED:-unknown}}
            echo "=== CARP: the annotation and its statistics are complete; side steps that failed: $FAILED" >> "$LOG"
        fi
        test -s "$OUT/Repeat_Annotation_Unified.gff3"
        test -s "$OUT/summary_statistics.csv"
        rm -rf "$RUN/tmp" "$RUN/sif_tmp"
        {{
          echo "carp_version: {params.version}"
          echo "image: {input.sif}"
          echo "genome: $GENOME"
          echo "finished: $(date -Iseconds)"
          echo "host: $(hostname)"
          echo "carp_exit: $rc"
          echo "carp_failed_steps: $FAILED"
        }} > "$REPO_ROOT/{output.done}"
        """


rule publish_repeats:
    """Freeze the repeat annotation for consumers (CO_smk) in refs/{species}_repeats/.

    Copies the unified annotation and the files needed to read it, writes
    subsets for rDNA and tandem repeats, and a MANIFEST recording the CARP
    version, the image checksum, the genome it ran on and the repo commit.
    """
    input:
        gff="results/{species}/repeats/carp/Repeat_Annotation_Unified.gff3",
        stats="results/{species}/repeats/carp/summary_statistics.csv",
        done="results/{species}/repeats/carp_run/carp.done",
        genome=lambda wc: frozen_chr_fasta(wc.species),
        sif=CARP_SIF,
    output:
        gff=REFS_ROOT + "/{species}_repeats/Repeat_Annotation_Unified.gff3",
        rdna=REFS_ROOT + "/{species}_repeats/rDNA.gff3",
        tandem=REFS_ROOT + "/{species}_repeats/tandem_repeats.gff3",
        stats=REFS_ROOT + "/{species}_repeats/summary_statistics.csv",
        seqs=REFS_ROOT + "/{species}_repeats/sequences.tsv",
        manifest=REFS_ROOT + "/{species}_repeats/MANIFEST.txt",
    params:
        outdir="results/{species}/repeats/carp",
        rundir="results/{species}/repeats/carp_run",
        version=CARP_VERSION,
    wildcard_constraints:
        species="|".join(REPEAT_SPECIES),
    resources:
        mem_mb=8000,
        cpus_per_task=2,
        slurm_extra="--time=02:00:00 --exclude=hpc15",
    log:
        "logs/publish_repeats/{species}.log",
    run:
        import collections, datetime, hashlib, os, re, shutil, subprocess
        from urllib.parse import unquote

        def top_class(cls):
            """Class_I/LTR/<superfamily>, rDNA/<45S_rDNA|5S_rDNA>, or the first word otherwise."""
            parts = cls.split("/")
            if parts[:2] == ["Class_I", "LTR"] and len(parts) > 2:
                return "/".join(parts[:3])
            if parts[0] in ("Class_I", "rDNA") and len(parts) > 1:
                return "/".join(parts[:2])
            return re.split(r"[(\[]", parts[0])[0] or "Unknown"

        dest = os.path.dirname(output.manifest)
        os.makedirs(dest, exist_ok=True)

        def sha256(path, chunk=1 << 24):
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for block in iter(lambda: fh.read(chunk), b""):
                    h.update(block)
            return h.hexdigest()

        shutil.copyfile(input.gff, output.gff)
        shutil.copyfile(input.stats, output.stats)

        # extra files worth having next to the annotation (copied if CARP wrote them)
        extras = {
            "TideCluster/default/TideCluster_rdna.tsv": "tidecluster_rdna.tsv",
            "TideCluster/default/TideCluster_report/data/trc_table.tsv": "trc_table_default.tsv",
            "TideCluster/short_monomer/TideCluster_report/data/trc_table.tsv": "trc_table_short.tsv",
            "TideCluster/default/TideCluster_consensus_dimer_library.fasta": "trc_consensus_dimers.fasta",
            "Repeat_Annotation_Unified.te_derived_trc.csv": "te_derived_trc.csv",
            "run_provenance.json": "carp_run_provenance.json",
            "carp_manifest.json": "carp_manifest.json",
        }
        copied = []
        for src_rel, dst_name in extras.items():
            src = os.path.join(params.outdir, src_rel)
            if os.path.exists(src):
                shutil.copyfile(src, os.path.join(dest, dst_name))
                copied.append(dst_name)
        cfg = os.path.join(params.rundir, "carp_config.yaml")
        if os.path.exists(cfg):
            shutil.copyfile(cfg, os.path.join(dest, "carp_config.yaml"))
            copied.append("carp_config.yaml")

        # rDNA and tandem-repeat subsets, and a count of bp by top-level class
        bp = collections.Counter()
        n_feat = collections.Counter()
        seqids = set()
        with open(input.gff) as src, open(output.rdna, "w") as rd, open(output.tandem, "w") as td:
            for line in src:
                if line.startswith("#"):
                    rd.write(line)
                    td.write(line)
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                attrs = {k: unquote(v) for k, v in
                         (kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)}
                cls = attrs.get("classification", "Unknown")
                seqids.add(f[0])
                if not attrs.get("ID", "").startswith("UA_L2_"):
                    top = top_class(cls)
                    bp[top] += int(f[4]) - int(f[3]) + 1
                    n_feat[top] += 1
                if cls.startswith("rDNA"):
                    rd.write(line)
                if cls.startswith(("Satellite", "rDNA/45S_rDNA", "rDNA/5S_rDNA")) and \
                        attrs.get("source_tool", "").startswith(("TideCluster", "TideHunter")):
                    td.write(line)

        # sequences the annotation was made on (name, length), from the genome index
        fai = input.genome + ".fai"
        if not os.path.exists(fai):
            subprocess.run("samtools faidx %s" % input.genome, shell=True, check=True)
        total = 0
        with open(fai) as fh, open(output.seqs, "w") as out:
            out.write("seqid\tlength\n")
            for line in fh:
                name, length = line.split("\t")[:2]
                out.write("%s\t%s\n" % (name, length))
                total += int(length)

        try:
            commit = subprocess.run("git rev-parse --short HEAD", shell=True,
                                    capture_output=True, text=True).stdout.strip()
            dirty = subprocess.run("git status --porcelain -- workflow config", shell=True,
                                   capture_output=True, text=True).stdout.strip()
        except Exception:
            commit, dirty = "unknown", ""

        with open(output.manifest, "w") as fh:
            fh.write("species:          %s\n" % wildcards.species)
            fh.write("what:             repeat and rDNA annotation (CARP) of both haplotypes\n")
            fh.write("published:        %s\n" % datetime.datetime.now().isoformat(timespec="seconds"))
            fh.write("carp_version:     %s\n" % params.version)
            fh.write("carp_image:       %s\n" % os.path.abspath(input.sif))
            fh.write("carp_image_sha256: %s\n" % sha256(input.sif))
            fh.write("source_fasta:     %s\n" % os.path.abspath(input.genome))
            fh.write("frozen_stage:     %s\n" % annotation_df.loc[wildcards.species, "frozen_stage"])
            fh.write("source_results:   %s\n" % os.path.abspath(params.outdir))
            fh.write("repo_commit:      %s%s\n" % (commit, "  (DIRTY WORKING TREE)" if dirty else ""))
            carp_done = dict(l.rstrip("\n").split(": ", 1) for l in open(input.done) if ": " in l)
            fh.write("carp_exit:        %s\n" % carp_done.get("carp_exit", "0"))
            fh.write("carp_failed_steps: %s (side steps; the files here do not depend on them)\n"
                     % carp_done.get("carp_failed_steps", "none"))
            fh.write("sequences:        %d (%.3f Gb); annotated: %d\n" % (
                sum(1 for _ in open(output.seqs)) - 1, total / 1e9, len(seqids)))
            fh.write("files:            Repeat_Annotation_Unified.gff3, rDNA.gff3, tandem_repeats.gff3, "
                     "summary_statistics.csv, sequences.tsv%s\n" % (
                         (", " + ", ".join(copied)) if copied else ""))
            fh.write("\nbp by class (top-level features of the unified annotation; Simple_repeat and\n"
                     "Low_complexity may overlap other classes, CARP's summary_statistics.csv is authoritative):\n")
            for k, v in sorted(bp.items(), key=lambda kv: -kv[1]):
                fh.write("  %-34s %14d bp  %6.2f%%  %9d features\n" % (k, v, 100.0 * v / max(total, 1), n_feat[k]))
            if dirty:
                fh.write("\nWARNING: published from a dirty working tree (workflow/ or config/);\n")
                fh.write("the commit above does not fully describe what produced this annotation.\n")

        with open(log[0], "w") as fh:
            fh.write("%s: %d sequences, %.3f Gb, rDNA %d bp, copied %s\n" % (
                wildcards.species, len(seqids), total / 1e9,
                sum(v for k, v in bp.items() if k.startswith("rDNA")), ", ".join(copied)))
