#!/usr/bin/env bash
# ===========================================================================
# Build the PM reference databases for the optional structure_search stack:
#   - an HH-suite profile DB  (hhsuite_pm_db)  from the PM-exclusive OG alignments
#   - a Foldseek structure DB (foldseek_pm_db) from predicted reference structures
#
# There is no public "pm_refs" database: both are built from the six PM-exclusive
# reference proteins (Mur alpha/beta/gamma/delta, the MraY-like GT, the CPS). This
# script builds whichever half you give it inputs for -- do the sequence-only
# HH-suite DB first (no GPU, no structures), the Foldseek DB once you have
# reference structures.
#
# The output prefixes are what go in config.yaml under
#   specificity.structure_search.{hhsuite_pm_db, foldseek_pm_db}
# and they are exactly the files structure_search.py / preflight check for
# (<prefix>_hhm.ffdata for HH-suite, <prefix>.dbtype for Foldseek).
#
# Usage (either half, or both):
#   setup/build_structure_dbs.sh \
#       [--alignments DIR --hhsuite-out PREFIX] \
#       [--structures DIR  --foldseek-out PREFIX] \
#       [--envs-root ROOT]
#
#   --alignments   directory of the Lupo PM-exclusive OG MSAs, e.g.
#                  Lupo-.../ompapa/alignments (files named *OG0001163*.fasta ...).
#   --hhsuite-out  output prefix, e.g. /ebio/.../hhsuite/pm_refs
#   --structures   directory of reference PDB/mmCIF structures (predict the six OG
#                  representatives with ESMFold/ColabFold; name them <OGID>.pdb so
#                  Foldseek hits are interpretable).
#   --foldseek-out output prefix, e.g. /ebio/.../foldseek/pm_refs
#   --envs-root    optional pre-built conda envs root (offline cluster). hhsuite +
#                  foldseek are activated from the `structure` env; hhsuite tools
#                  (reformat.pl, hhmake, cstranslate, hhsuitedb.py) must be present.
#
# Offline: nothing here reaches the network.
# ===========================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ALN="" HH_OUT="" STRUCT="" FS_OUT="" ENVS_ROOT=""
while [ $# -gt 0 ]; do
    case "$1" in
        --alignments)   ALN="$2"; shift 2;;
        --hhsuite-out)  HH_OUT="$2"; shift 2;;
        --structures)   STRUCT="$2"; shift 2;;
        --foldseek-out) FS_OUT="$2"; shift 2;;
        --envs-root)    ENVS_ROOT="$2"; shift 2;;
        *) echo "unknown arg: $1" >&2; exit 2;;
    esac
done
if [ -z "$HH_OUT$FS_OUT" ]; then
    echo "Give at least one of: (--alignments + --hhsuite-out) or (--structures + --foldseek-out)" >&2
    exit 2
fi

act() { [ -n "$ENVS_ROOT" ] && { set +u; source "$ENVS_ROOT/structure/bin/activate"; set -u; }; return 0; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "[build_structure_dbs] required tool not on PATH: $1" >&2; exit 3; }; }

# The six PM-exclusive OGs, from pmur_reference so this stays in sync.
OGS=$(python3 - "$HERE" <<'PY'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
import pmur_reference as R
print(" ".join(R.PM_EXCLUSIVE.keys()))
PY
)
echo "[build_structure_dbs] PM-exclusive OGs: $OGS"
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT

# ---- HH-suite profile DB ---------------------------------------------------
if [ -n "$HH_OUT" ]; then
    [ -n "$ALN" ] || { echo "--hhsuite-out needs --alignments" >&2; exit 2; }
    act; need reformat.pl; need hhmake
    mkdir -p "$(dirname "$HH_OUT")" "$work/a3m"
    n=0
    for og in $OGS; do
        f="$(find "$ALN" -maxdepth 1 -type f \( -iname "*${og}*.fasta" -o -iname "*${og}*.fas" -o -iname "*${og}*.a3m" \) | head -n1 || true)"
        if [ -z "$f" ]; then echo "[build_structure_dbs] MISSING alignment for $og in $ALN" >&2; continue; fi
        case "$f" in
            *.a3m) cp "$f" "$work/a3m/${og}.a3m";;
            *)     reformat.pl fas a3m "$f" "$work/a3m/${og}.a3m" >/dev/null 2>&1;;
        esac
        n=$((n+1)); echo "  a3m: $og  <- $(basename "$f")"
    done
    [ "$n" -ge 1 ] || { echo "[build_structure_dbs] no alignments found; nothing to build" >&2; exit 4; }
    # hhsuitedb.py is the canonical builder: it runs hhmake + cstranslate and emits
    # <prefix>_{a3m,hhm,cs219}.ff{data,index} -- exactly what stage_ready checks.
    if command -v hhsuitedb.py >/dev/null 2>&1; then
        hhsuitedb.py --ia3m "$work/a3m/*.a3m" -o "$HH_OUT"
    else
        echo "[build_structure_dbs] hhsuitedb.py not found; building the ffindex DB by hand" >&2
        need ffindex_build; need cstranslate
        ( cd "$work/a3m" && for a in *.a3m; do hhmake -i "$a" -o "${a%.a3m}.hhm" >/dev/null 2>&1; done )
        ffindex_build -s "${HH_OUT}_a3m.ffdata" "${HH_OUT}_a3m.ffindex" "$work/a3m"/*.a3m
        ffindex_build -s "${HH_OUT}_hhm.ffdata" "${HH_OUT}_hhm.ffindex" "$work/a3m"/*.hhm
        # cs219 needs the context libraries shipped with hh-suite ($HHLIB/data)
        cs_lib="${HHLIB:-}/data/cs219.lib"; ctx_lib="${HHLIB:-}/data/context_data.lib"
        cstranslate -f -x 0.3 -c 4 -I a3m -i "${HH_OUT}_a3m" -o "${HH_OUT}_cs219" \
            ${HHLIB:+-A "$cs_lib" -D "$ctx_lib"}
    fi
    if [ -s "${HH_OUT}_hhm.ffdata" ]; then
        echo "[build_structure_dbs] OK  HH-suite DB -> ${HH_OUT}  (set hhsuite_pm_db: ${HH_OUT})"
    else
        echo "[build_structure_dbs] FAIL: ${HH_OUT}_hhm.ffdata was not produced" >&2; exit 5
    fi
fi

# ---- Foldseek structure DB -------------------------------------------------
if [ -n "$FS_OUT" ]; then
    [ -n "$STRUCT" ] || { echo "--foldseek-out needs --structures" >&2; exit 2; }
    act; need foldseek
    ns=$(find "$STRUCT" -maxdepth 1 -type f \( -iname "*.pdb" -o -iname "*.cif" -o -iname "*.mmcif" -o -iname "*.pdb.gz" -o -iname "*.cif.gz" \) | wc -l)
    [ "$ns" -ge 1 ] || { echo "[build_structure_dbs] no structures in $STRUCT (need .pdb/.cif)" >&2; exit 4; }
    echo "[build_structure_dbs] $ns reference structures -> foldseek createdb"
    mkdir -p "$(dirname "$FS_OUT")"
    foldseek createdb "$STRUCT" "$FS_OUT"
    foldseek createindex "$FS_OUT" "$work/fs_tmp" >/dev/null 2>&1 || true   # optional speedup
    if [ -s "${FS_OUT}.dbtype" ]; then
        echo "[build_structure_dbs] OK  Foldseek DB -> ${FS_OUT}  (set foldseek_pm_db: ${FS_OUT})"
    else
        echo "[build_structure_dbs] FAIL: ${FS_OUT}.dbtype was not produced" >&2; exit 5
    fi
fi

echo
echo "[build_structure_dbs] done. In config.yaml under specificity.structure_search set:"
[ -n "$HH_OUT" ] && echo "    hhsuite_pm_db:  $HH_OUT"
[ -n "$FS_OUT" ] && echo "    foldseek_pm_db: $FS_OUT"
echo "  then set enabled: true and re-run: python scripts/preflight.py"
