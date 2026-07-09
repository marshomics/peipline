#!/usr/bin/env bash
# ============================================================================
# STEP 1 of 2 — run this on a machine WITH internet.
#
# Builds the six conda environments from conda-forge and bioconda only, then
# packs each into a relocatable tarball with conda-pack. Nothing from the
# `defaults` / `anaconda` / `main` / `r` channels is used or contacted.
#
#   ./setup/build_envs.sh [OUTDIR]        # default: ./dist
#
# The tarballs are glibc- and architecture-specific. Build them on a machine
# whose OS is no newer than the cluster's, and on the same architecture
# (linux-64 for almost every SGE cluster). If your laptop is macOS or arm64,
# build inside a container instead:
#
#   docker run --rm -v "$PWD":/w -w /w --platform linux/amd64 \
#       condaforge/miniforge3:latest ./setup/build_envs.sh /w/dist
#
# Then copy dist/ to the cluster and run setup/install_envs.sh there.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."

OUT="${1:-$PWD/dist}"
BUILD="$OUT/build"
PACK="$OUT/envs"
mkdir -p "$BUILD" "$PACK"

ENVS=(hmmer py phylo network prodigal r)

# ---------------------------------------------------------------------------
# Pick a solver. mamba/micromamba are much faster; conda works.
if command -v micromamba &> /dev/null; then
    SOLVER="micromamba"
elif command -v mamba &> /dev/null; then
    SOLVER="mamba"
elif command -v conda &> /dev/null; then
    SOLVER="conda"
else
    echo "ERROR: need micromamba, mamba or conda on PATH." >&2
    echo "  Miniforge is the conda-forge-only installer, and is what you want:" >&2
    echo "  https://github.com/conda-forge/miniforge" >&2
    exit 1
fi
echo "[build_envs] solver: $SOLVER"

# --override-channels is the load-bearing flag. Without it, conda appends
# whatever is in your ~/.condarc, which on a stock Anaconda install is the
# `defaults` channel. `nodefaults` in the env YAML covers the same ground; both
# are set, because either one alone can be undone by a stray .condarc.
CHAN=(--override-channels -c conda-forge -c bioconda)
export CONDA_CHANNEL_PRIORITY=strict

for e in "${ENVS[@]}"; do
    PREFIX="$BUILD/$e"
    if [ -d "$PREFIX" ]; then
        echo "[build_envs] $e already built, skipping"
    else
        echo "[build_envs] === solving $e ==="
        case "$SOLVER" in
          micromamba) micromamba create -y -p "$PREFIX" "${CHAN[@]}" -f "envs/$e.yaml" ;;
          mamba)      mamba env create -p "$PREFIX" -f "envs/$e.yaml" ;;
          conda)      conda env create -p "$PREFIX" -f "envs/$e.yaml" ;;
        esac
    fi
done

# conda-pack itself, in its own env so it cannot perturb the six
if ! command -v conda-pack &> /dev/null; then
    echo "[build_envs] installing conda-pack"
    case "$SOLVER" in
      micromamba) micromamba create -y -p "$BUILD/_pack" "${CHAN[@]}" conda-pack
                  export PATH="$BUILD/_pack/bin:$PATH" ;;
      *)          conda install -y "${CHAN[@]}" -n base conda-pack ;;
    esac
fi

for e in "${ENVS[@]}"; do
    TGZ="$PACK/$e.tar.gz"
    [ -f "$TGZ" ] && { echo "[build_envs] $e already packed"; continue; }
    echo "[build_envs] === packing $e ==="
    conda-pack -p "$BUILD/$e" -o "$TGZ" --n-threads 4
done

# ---------------------------------------------------------------------------
# Snakemake and its SGE executor plugin, as wheels, for offline pip install.
# Snakemake runs on the submit host, not on the compute nodes, so it does not
# need to be a conda env -- but it does need to be installed without internet.
echo "[build_envs] === downloading snakemake wheels ==="
mkdir -p "$OUT/pip"
python3 -m pip download --dest "$OUT/pip" \
    "snakemake>=8.20,<9" snakemake-executor-plugin-cluster-generic pyyaml

# provenance: exactly what got solved
mkdir -p "$OUT/locks"
for e in "${ENVS[@]}"; do
    case "$SOLVER" in
      micromamba) micromamba list -p "$BUILD/$e" --json > "$OUT/locks/$e.json" ;;
      *)          conda list -p "$BUILD/$e" --explicit > "$OUT/locks/$e.txt" ;;
    esac
done

cat > "$OUT/MANIFEST.txt" <<EOF
built:    $(date -u +%Y-%m-%dT%H:%M:%SZ)
host:     $(uname -srm)
solver:   $SOLVER
channels: conda-forge, bioconda (override-channels, strict priority, nodefaults)
envs:     ${ENVS[*]}
EOF

echo
echo "[build_envs] done. Copy this to the cluster:"
du -sh "$OUT"
echo "  rsync -a $OUT/ user@cluster:/path/to/c71_dist/"
echo "  then on the cluster: ./setup/install_envs.sh /path/to/c71_dist"
