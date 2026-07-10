#!/usr/bin/env bash
# ============================================================================
# STEP 2 of 2 — run this on the cluster. No internet required.
#
#   ./setup/install_envs.sh /path/to/c71_dist [ENVS_ROOT]
#
# Unpacks the conda-pack tarballs, rewrites their absolute paths with
# conda-unpack, installs snakemake from the vendored wheels, and verifies that
# every binary the pipeline calls is actually present and runnable.
#
# ENVS_ROOT defaults to <repo>/envs_installed. Put it somewhere the compute
# nodes can read: a shared filesystem, not /tmp on the submit host.
# Then set `envs_root:` in config.yaml to that path.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
REPO="$PWD"

DIST="${1:?usage: install_envs.sh /path/to/c71_dist [ENVS_ROOT]}"
ENVS_ROOT="${2:-$REPO/envs_installed}"
# `selection` was missing here while the Snakefile REQUIRED it, so an
# offline cluster had no supported way to obtain it: the run aborted
# telling the user to run this script, which did not install it.
ENVS=(hmmer py phylo network prodigal r selection)

[ -d "$DIST/envs" ] || { echo "ERROR: $DIST/envs not found. Did build_envs.sh run?" >&2; exit 1; }
mkdir -p "$ENVS_ROOT"

for e in "${ENVS[@]}"; do
    TGZ="$DIST/envs/$e.tar.gz"
    DST="$ENVS_ROOT/$e"
    [ -f "$TGZ" ] || { echo "ERROR: missing $TGZ" >&2; exit 1; }
    if [ -x "$DST/bin/conda-unpack" ] || [ -f "$DST/.unpacked" ]; then
        echo "[install_envs] $e already installed"
        continue
    fi
    echo "[install_envs] === $e ==="
    rm -rf "$DST"; mkdir -p "$DST"
    tar -xzf "$TGZ" -C "$DST"
    # rewrites the build-machine prefix baked into scripts, shebangs and rpaths
    ( set +u; source "$DST/bin/activate"; conda-unpack; )
    touch "$DST/.unpacked"
done

# --- snakemake on the submit host, from vendored wheels ---------------------
if [ -d "$DIST/pip" ]; then
    echo "[install_envs] === snakemake (offline pip) ==="
    python3 -m pip install --no-index --find-links "$DIST/pip" \
        "snakemake>=8.20,<9" snakemake-executor-plugin-cluster-generic \
        || echo "[install_envs] WARNING: pip install failed; install snakemake yourself"
fi

# --- verify -----------------------------------------------------------------
echo
echo "[install_envs] === verifying ==="
declare -A NEED=(
  [hmmer]="hmmsearch hmmalign hmmconvert"
  [prodigal]="prodigal"
  [phylo]="iqtree2 trimal seqkit mmseqs python"
  [network]="diamond mmseqs seqkit python"
  [py]="python"
  [r]="Rscript"
  [selection]="hyphy iqtree2 python"
)

# Binaries are not enough. module_trees.py runs in `phylo` and imports numpy,
# pandas, Bio and matplotlib; selection.py runs in `selection` and imports
# numpy/pandas/scipy/statsmodels. A binary-only check passed happily while
# `import numpy` was about to fail on a compute node, hours into a run.
declare -A PYIMPORTS=(
  [py]="numpy pandas scipy sklearn matplotlib Bio logomaker statsmodels yaml"
  [phylo]="numpy pandas Bio matplotlib yaml"
  [hmmer]="numpy pandas yaml"
  [network]="numpy pandas scipy sklearn igraph matplotlib yaml"
  [prodigal]="yaml"
  [selection]="numpy pandas scipy statsmodels matplotlib yaml"
)
FAIL=0
for e in "${ENVS[@]}"; do
    for tool in ${NEED[$e]}; do
        if ( set +u; source "$ENVS_ROOT/$e/bin/activate"; command -v "$tool" &>/dev/null ); then
            printf '  [ok]   %-9s %s\n' "$e" "$tool"
        else
            printf '  [FAIL] %-9s %s\n' "$e" "$tool"; FAIL=1
        fi
    done
done

echo "  --- python imports, per environment ---"
for e in "${ENVS[@]}"; do
    mods="${PYIMPORTS[$e]:-}"
    [ -z "$mods" ] && continue
    ( set +u; source "$ENVS_ROOT/$e/bin/activate"
      MODS="$mods" ENVNAME="$e" python -c "
import importlib.util, os, sys
mods = os.environ['MODS'].split()
env = os.environ['ENVNAME']
bad = [m for m in mods if importlib.util.find_spec(m) is None]
if bad:
    print(f'  [FAIL] {env:<9} missing imports: {bad}')
    sys.exit(1)
print(f'  [ok]   {env:<9} imports')
" ) || FAIL=1
done

echo "  --- R packages ---"
( set +u; source "$ENVS_ROOT/r/bin/activate"
  Rscript -e 'p <- c("ape","phylolm","caper","data.table","yaml")
              m <- p[!sapply(p, requireNamespace, quietly=TRUE)]
              if (length(m)) { cat("  [FAIL] missing:", m, "\n"); quit(status=1) }
              cat("  [ok]   R packages\n")' ) || FAIL=1

echo
if [ "$FAIL" -ne 0 ]; then
    echo "[install_envs] FAILED. Do not run the pipeline." >&2
    exit 1
fi

echo "[install_envs] all environments OK."
echo
echo "Now set this in config.yaml:"
echo
echo "    envs_root: $ENVS_ROOT"
echo
echo "and run:  ./run.sh preflight"
