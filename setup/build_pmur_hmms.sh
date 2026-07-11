#!/usr/bin/env bash
# ===========================================================================
# Build the PM-exclusive marker HMMs for the out-of-order pseudomurein screen.
#
# Taxonomy-first design (Lupo et al. 2025): the primary PM call is GTDB order,
# with no markers. These six HMMs exist for ONE job -- the out-of-order discovery
# block that hunts PM in lineages OUTSIDE Methanobacteriales/Methanopyrales. The
# block only fires when >=1 muramyl ligase + the MraY-like GT + the CPS are all
# present, so build all six.
#
# You supply one protein FASTA per orthogroup (the member sequences from Lupo
# 2025's OrthoFinder output, Orthogroup_Sequences/OG000....fa). This aligns each
# with MAFFT and runs hmmbuild, naming the output <OGID>.hmm so cellwall_genotype's
# classify_markers() picks up the role automatically.
#
# Usage:
#   setup/build_pmur_hmms.sh  INPUT_FAA_DIR  [OUTPUT_HMM_DIR]  [ENVS_ROOT]
#
#   INPUT_FAA_DIR   directory containing the six OG FASTAs. Each file's name must
#                   contain its OG id, e.g. OG0001163.fa / OG0001163_MraY.faa.
#   OUTPUT_HMM_DIR  where to write the .hmm files. Default: the pmur_hmm_dir in
#                   config.yaml.
#   ENVS_ROOT       optional: pre-built conda envs root (offline cluster). If set,
#                   the hmmer and phylo envs are activated for hmmbuild and mafft.
#                   If unset, mafft and hmmbuild must already be on PATH.
#
# Offline note: nothing here reaches the network. It needs mafft (phylo env) and
# hmmbuild (hmmer env).
# ===========================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN_DIR="${1:?INPUT_FAA_DIR required (one FASTA per OG, OG id in the filename)}"
OUT_DIR="${2:-}"
ENVS_ROOT="${3:-}"

# The six PM-exclusive OGs and their role, read from pmur_reference so this stays
# in sync with the catalogue. Also resolves pmur_hmm_dir from config as the default.
read -r DEFAULT_OUT OG_LIST < <(python3 - "$HERE" <<'PY'
import sys, os, yaml
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
import pmur_reference as R
cfg = yaml.safe_load(open(os.path.join(sys.argv[1], "config.yaml")))
default_out = cfg["specificity"]["pmur_hmm_dir"]
ogs = list(R.PM_EXCLUSIVE.keys())          # OG0001148 ... OG0001014
print(default_out, ",".join(ogs))
PY
)
OUT_DIR="${OUT_DIR:-$DEFAULT_OUT}"
IFS=',' read -r -a OGS <<< "$OG_LIST"

act() { [ -n "$ENVS_ROOT" ] && { set +u; source "$ENVS_ROOT/$1/bin/activate"; set -u; }; }

echo "[build_pmur] output dir: $OUT_DIR"
echo "[build_pmur] required OGs: ${OGS[*]}"
mkdir -p "$OUT_DIR"
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT

built=0
for og in "${OGS[@]}"; do
    # locate an input FASTA whose name contains this OG id
    faa="$(find "$IN_DIR" -maxdepth 1 -type f \
             \( -iname "*${og}*.fa" -o -iname "*${og}*.faa" -o -iname "*${og}*.fasta" \) \
             | head -n1 || true)"
    if [ -z "$faa" ]; then
        echo "[build_pmur] MISSING input for $og in $IN_DIR (looked for *${og}*.{fa,faa,fasta})" >&2
        continue
    fi
    n_seq="$(grep -c '^>' "$faa" || echo 0)"
    aln="$work/${og}.aln"
    out="$OUT_DIR/${og}.hmm"
    echo "[build_pmur] $og: $n_seq sequences from $(basename "$faa")"
    ( act phylo; mafft --auto --anysymbol "$faa" > "$aln" 2> "$work/${og}.mafft.log" )
    ( act hmmer; hmmbuild --amino -n "$og" -o "$work/${og}.hmmbuild.log" "$out" "$aln" )
    built=$((built + 1))
done

echo
echo "[build_pmur] built $built/${#OGS[@]} HMMs into $OUT_DIR"

# Show the role classification the pipeline will infer, and whether the block is
# complete -- so you catch a naming/role problem here, not six hours into a run.
python3 - "$HERE" "$OUT_DIR" <<'PY'
import sys, os, glob
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
from cellwall_genotype import classify_markers
hmms = sorted(glob.glob(os.path.join(sys.argv[2], "*.hmm")))
stems = [os.path.splitext(os.path.basename(h))[0] for h in hmms]
roles = classify_markers(stems)
print("\n[build_pmur] role classification:")
for m in stems:
    print(f"    {m:20s} -> {roles[m]}")
lig = [m for m, r in roles.items() if r == "muramyl_ligase"]
mray = [m for m, r in roles.items() if r == "mray_like"]
cps = [m for m, r in roles.items() if r == "cps"]
trap = [m for m, r in roles.items() if r == "trap"]
block_ok = bool(lig) and bool(mray) and bool(cps)
print(f"\n[build_pmur] ligases={lig} mray={mray} cps={cps}")
if trap:
    print(f"[build_pmur] WARNING: TRAP families present (remove them): {trap}")
if block_ok:
    print("[build_pmur] OK: the PM-exclusive block is complete; out-of-order "
          "discovery is enabled.")
else:
    print("[build_pmur] INCOMPLETE block: out-of-order discovery will fall back to "
          "marker count (weaker). Need >=1 muramyl ligase + MraY-like + CPS.")
    sys.exit(1)
PY

echo
echo "[build_pmur] done. Re-run preflight to confirm: python scripts/preflight.py"
