# ---------------------------------------------------------------------------
# Publishing references for downstream repos (CO_smk, etc.)
#
# reproducible_phd PRODUCES assemblies and annotations; other repos CONSUME
# them. Rather than have a consumer reach into this repo's results/ tree (which
# couples it to our internal layout), we publish a frozen, versioned artifact
# with a manifest recording exactly what produced it. Consumers point at
# refs/ via absolute paths and never at results/.
# ---------------------------------------------------------------------------

REFS_ROOT = config.get("refs_root",
                       "/netscratch/dep_mercier/grp_marques/Aaryan/refs")


def _hap_contigs(fai_path, hap):
    """Contig names for one haplotype, in numeric chromosome order.

    Assemblies are phased dual-haplotype: chr1_hap1, chr1_hap2, ... A single
    haplotype must be extracted before the reference can be used for variant
    calling -- against the dual-hap assembly every site looks HOMOZYGOUS,
    because the two haplotypes are separate contigs rather than two alleles at
    one locus.
    """
    names = []
    with open(fai_path) as fh:
        for line in fh:
            c = line.split("\t")[0]
            if c.endswith(f"_{hap}"):
                names.append(c)

    def key(c):
        stem = c[:-(len(hap) + 1)]
        stem = stem[3:] if stem.startswith("chr") else stem
        try:
            return (0, int(stem), "")
        except ValueError:
            return (1, 0, stem)

    return sorted(names, key=key)


rule publish_haplotype_reference:
    """Freeze a single-haplotype reference (FASTA + GFF3 + manifest) into refs/.

    Downstream repos consume refs/{species}_{hap}/ and never results/. The
    manifest records the source paths and the git commit, so a published
    reference can always be traced back to what produced it.
    """
    input:
        fasta=lambda wc: frozen_chr_fasta(wc.species),
        gff="results/{species}/annotation/final/{species}.final.gff3",
    output:
        fasta=REFS_ROOT + "/{species}_{hap}/genome.fa",
        fai=REFS_ROOT + "/{species}_{hap}/genome.fa.fai",
        gff=REFS_ROOT + "/{species}_{hap}/annotation.gff3",
        names=REFS_ROOT + "/{species}_{hap}/contigs.txt",
        manifest=REFS_ROOT + "/{species}_{hap}/MANIFEST.txt",
    wildcard_constraints:
        hap="hap[12]",
    log:
        "logs/publish_haplotype_reference/{species}.{hap}.log",
    resources:
        mem_mb=8000,
        runtime=120,
    run:
        import os, subprocess, datetime, shutil

        fai_src = input.fasta + ".fai"
        if not os.path.exists(fai_src):
            subprocess.run(f"samtools faidx {input.fasta}", shell=True, check=True)

        contigs = _hap_contigs(fai_src, wildcards.hap)
        if not contigs:
            raise ValueError(
                f"no contigs ending in _{wildcards.hap} in {fai_src}. "
                f"Is this assembly phased dual-haplotype?")

        os.makedirs(os.path.dirname(output.fasta), exist_ok=True)
        with open(output.names, "w") as fh:
            fh.write("\n".join(contigs) + "\n")

        # FASTA: extract the haplotype's contigs, then index
        subprocess.run(
            f"samtools faidx {input.fasta} {' '.join(contigs)} > {output.fasta}",
            shell=True, check=True)
        subprocess.run(f"samtools faidx {output.fasta}", shell=True, check=True)

        # GFF3: keep only features on those contigs. An annotation covering
        # contigs absent from the index gives STAR junctions it cannot place.
        keep = set(contigs)
        n_feat = n_gene = 0
        with open(input.gff) as src, open(output.gff, "w") as dst:
            for line in src:
                if line.startswith("#"):
                    dst.write(line)
                    continue
                f = line.split("\t", 3)
                if len(f) > 2 and f[0] in keep:
                    dst.write(line)
                    n_feat += 1
                    if f[2] == "gene":
                        n_gene += 1

        bp = 0
        with open(output.fasta + ".fai") as fh:
            for line in fh:
                bp += int(line.split("\t")[1])

        try:
            commit = subprocess.run("git rev-parse --short HEAD", shell=True,
                                    capture_output=True, text=True).stdout.strip()
            dirty = subprocess.run("git status --porcelain", shell=True,
                                   capture_output=True, text=True).stdout.strip()
        except Exception:
            commit, dirty = "unknown", ""

        with open(output.manifest, "w") as fh:
            fh.write(f"species:          {wildcards.species}\n")
            fh.write(f"haplotype:        {wildcards.hap}\n")
            fh.write(f"published:        {datetime.datetime.now().isoformat(timespec='seconds')}\n")
            fh.write(f"source_fasta:     {os.path.abspath(input.fasta)}\n")
            fh.write(f"source_gff:       {os.path.abspath(input.gff)}\n")
            fh.write(f"frozen_stage:     {annotation_df.loc[wildcards.species, 'frozen_stage']}\n")
            fh.write(f"repo_commit:      {commit}{'  (DIRTY WORKING TREE)' if dirty else ''}\n")
            fh.write(f"contigs:          {len(contigs)}\n")
            fh.write(f"total_bp:         {bp}\n")
            fh.write(f"gff_features:     {n_feat}\n")
            fh.write(f"gff_genes:        {n_gene}\n")
            if dirty:
                fh.write("\nWARNING: published from a dirty working tree; the commit\n")
                fh.write("above does not fully describe what produced this reference.\n")

        with open(log[0], "w") as fh:
            fh.write(f"{wildcards.species} {wildcards.hap}: {len(contigs)} contigs, "
                     f"{bp/1e9:.3f} Gb, {n_gene} genes\n")
        print(f"published {wildcards.species}_{wildcards.hap}: "
              f"{len(contigs)} contigs, {bp/1e9:.3f} Gb, {n_gene} genes")
