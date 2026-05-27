# HapHiC install (for this pipeline)

HapHiC's official PyPI/bioconda packages are unmaintained per their README.
We install from upstream git, into a custom conda env with the alignment tools
needed for our scaffolding rules.

## Steps

```bash
# 1. Clone HapHiC (from repo root)
mkdir -p methods
cd methods
git clone https://github.com/zengxiaofei/HapHiC.git
cd HapHiC

# Pin commit for reproducibility (commit hash is in envs/haphic.git_commit)
git checkout $(cat ../../envs/haphic.git_commit)
cd ../..

# 2. Create the conda env (alignment tools + HapHiC deps in one env)
#    samtools comes in old via bioconda channel priority issues; we update
#    it post-create to a modern version.
micromamba env create -f envs/haphic.yaml \
  -p /path/to/your/micromamba_envs/haphic -y

# 3. Update samtools to a modern version (required for rule shell commands)
micromamba activate /path/to/your/micromamba_envs/haphic
micromamba update -y samtools

# 4. Verify
methods/HapHiC/haphic check       # should exit 0, all checks Successful
samtools --version | head -1      # should be 1.20+
which bwa samblaster              # both should resolve
```

## Pinned versions

- HapHiC commit: see `envs/haphic.git_commit`
- env definition: `envs/haphic.yaml`
- Verified working: samtools 1.23.1, bwa 0.7.19, samblaster 0.1.26, haphic 1.0.7
