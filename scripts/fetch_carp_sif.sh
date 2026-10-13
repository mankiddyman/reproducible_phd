#!/bin/bash
# Fetch the CARP container image (a Singularity/Apptainer SIF that the CARP authors
# publish to GHCR as an ORAS artefact) into methods/carp/. Run on the LOGIN node:
# compute nodes need no internet. Does nothing if the image is already there.
# Usage: bash scripts/fetch_carp_sif.sh [version]          (default 1.9.5)
VER="${1:-1.9.5}"
DEST=/netscratch/dep_mercier/grp_marques/Aaryan/methods/carp
SIF="$DEST/carp_${VER}.sif"
REPO=kavonrtep/carp/sif
mkdir -p "$DEST/cache" || exit 1
if [ -s "$SIF" ]; then
    echo "image already there: $SIF ($(du -h "$SIF" | cut -f1))"
    exit 0
fi
module load singularity-binds 2>/dev/null || module load singularity 2>/dev/null || true
export APPTAINER_CACHEDIR="$DEST/cache" SINGULARITY_CACHEDIR="$DEST/cache"
rm -f "$SIF.part"
for rt in apptainer singularity; do
    command -v "$rt" >/dev/null 2>&1 || continue
    echo "pulling oras://ghcr.io/$REPO:$VER with $rt ($("$rt" --version 2>&1))"
    if "$rt" pull --force "$SIF.part" "oras://ghcr.io/$REPO:$VER"; then
        mv "$SIF.part" "$SIF"
        echo "image: $SIF ($(du -h "$SIF" | cut -f1))"
        exit 0
    fi
    echo "$rt pull failed; trying the next way"
done
# Fallback: the registry's HTTP API with curl (anonymous pull token); the file is
# checked against the layer's sha256 digest before it is kept.
echo "fetching with curl from the GHCR registry API"
TOKEN=$(curl -fsSL "https://ghcr.io/token?service=ghcr.io&scope=repository:$REPO:pull" \
        | python3 -c 'import sys, json; print(json.load(sys.stdin)["token"])') \
    || { echo "could not get a registry token"; exit 1; }
DIGEST=$(curl -fsSL -H "Authorization: Bearer $TOKEN" \
        -H "Accept: application/vnd.oci.image.manifest.v1+json" \
        "https://ghcr.io/v2/$REPO/manifests/$VER" \
        | python3 -c 'import sys, json; m = json.load(sys.stdin); print(max(m["layers"], key=lambda l: l["size"])["digest"])') \
    || { echo "no manifest for $REPO:$VER"; exit 1; }
echo "layer $DIGEST"
curl -fL --retry 3 -H "Authorization: Bearer $TOKEN" -o "$SIF.part" \
    "https://ghcr.io/v2/$REPO/blobs/$DIGEST" || { echo "download failed"; exit 1; }
GOT="sha256:$(sha256sum "$SIF.part" | cut -d' ' -f1)"
if [ "$GOT" != "$DIGEST" ]; then
    echo "checksum mismatch: got $GOT, expected $DIGEST"
    exit 1
fi
mv "$SIF.part" "$SIF"
echo "image: $SIF ($(du -h "$SIF" | cut -f1))"
