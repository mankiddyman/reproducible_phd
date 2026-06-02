#!/bin/bash
# Telomere-flipping curation loop helper (NOT a snakemake rule).
#
#   teloflip.sh <species> [stage]
#
# 1. Finds the newest out_JBAT.review.assembly under
#    results/<species>/scaffolding/[<stage>/]  (stage default: search all)
# 2. Derives sibling out_JBAT.liftover.agp + the work/ ref fasta
# 3. juicer post -> teloflip.FINAL.fa (overwritten each iteration)
# 4. Reads n_scaffolds = chr_number_2n from config/species.csv
# 5. module load genespace; Rscript plot_telomeres.R (x11 over ssh -X)
#
# Loop: curate/flip in Juicebox -> export out_JBAT.review.assembly
#       -> teloflip.sh <species> -> inspect -> repeat.
set -euo pipefail

SPECIES="${1:?usage: teloflip.sh <species> [stage]}"
STAGE="${2:-}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"

HAPHIC_JUICER="$REPO_ROOT/methods/HapHiC/utils/juicer"
PLOT_R="$REPO_ROOT/workflow/scripts/telomere/plot_telomeres.R"

# --- locate newest review.assembly ---
if [ -n "$STAGE" ]; then
    SEARCH="results/$SPECIES/scaffolding/$STAGE"
else
    SEARCH="results/$SPECIES/scaffolding"
fi
REVIEW=$(find "$SEARCH" -name "out_JBAT.review.assembly" -printf '%T@ %p\n' 2>/dev/null \
         | sort -rn | head -1 | cut -d' ' -f2-)
if [ -z "$REVIEW" ]; then
    echo "ERROR: no out_JBAT.review.assembly under $SEARCH" >&2
    exit 1
fi
BUILD_DIR="$(dirname "$REVIEW")"
PASS_DIR="$(dirname "$BUILD_DIR")"
LIFTOVER="$BUILD_DIR/out_JBAT.liftover.agp"
REF_FA=$(ls "$PASS_DIR"/work/*.fa 2>/dev/null | head -1)

echo "Picked review : $REVIEW"
echo "Liftover      : $LIFTOVER"
echo "Ref fasta     : $REF_FA"
[ -f "$LIFTOVER" ] || { echo "ERROR: liftover not found: $LIFTOVER" >&2; exit 1; }
[ -n "$REF_FA" ] && [ -f "$REF_FA" ] || { echo "ERROR: no ref fasta in $PASS_DIR/work/" >&2; exit 1; }

# --- n_scaffolds from species.csv (chr_number_2n) ---
N_SCAFF=$(python3 -c "
import pandas as pd, sys
d = pd.read_csv('config/species.csv').fillna('').set_index('species_id')
v = str(d.loc['$SPECIES','chr_number_2n']).strip()
print(int(float(v)) if v else sys.exit('chr_number_2n empty for $SPECIES'))
")
echo "n_scaffolds   : $N_SCAFF (chr_number_2n)"

# --- juicer post -> teloflip.FINAL.fa (overwrite working copy) ---
OUT_PREFIX="$BUILD_DIR/teloflip"
echo "=== juicer post ==="
export PATH=/usr/bin:$PATH
"$HAPHIC_JUICER" post -o "$OUT_PREFIX" "$REVIEW" "$LIFTOVER" "$REF_FA"
FINAL_FA="${OUT_PREFIX}.FINAL.fa"
[ -f "$FINAL_FA" ] || { echo "ERROR: juicer post did not produce $FINAL_FA" >&2; exit 1; }
echo "Working FINAL : $FINAL_FA"

# --- plot ---
echo "=== telomere plot (press Enter in this terminal when done viewing) ==="
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load genespace
Rscript "$PLOT_R" "$FINAL_FA" "$N_SCAFF"
