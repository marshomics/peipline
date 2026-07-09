#!/usr/bin/env bash
# ============================================================================
# One batch: concatenate its proteomes, build the reversed decoys, and run
# hmmsearch for every profile against both -- all on node-local scratch.
#
#   search_batch.sh BATCH_TSV BATCH_ID PROFILES_DIR HMM_OUT_DIR MAP_OUT \
#                   THREADS CONFIG DO_DECOY
#
# Why this is one job instead of four rules.
#
# The batch FASTA is pure intermediate. At 1,000 genomes per batch it is roughly
# 1 GB, and the reversed decoy is another 1 GB. Across ~365 batches that is
# ~700 GB written to the shared filesystem and then read straight back by
# hmmsearch, for data that no other job will ever look at. On a cluster where
# the proteomes already live on shared storage, that write-then-read is the
# dominant cost of the entire screen: hmmsearch itself only needs a few tens of
# core-hours for 3x10^11 residues.
#
# Fusing the four steps lets the intermediate live in $TMPDIR on the execution
# host, which SGE creates per job and removes on exit. Only the domtblout files
# and the per-sample map -- kilobytes, not gigabytes -- reach shared storage.
#
# The cost is granularity: re-running hmmsearch now re-reads the proteomes. That
# is the right trade, because re-reading the proteomes was always the expensive
# half.
# ============================================================================
set -euo pipefail

BATCH_TSV="${1:?batch tsv}"
BATCH_ID="${2:?batch id}"
PROFILES_DIR="${3:?profiles dir}"
HMM_OUT="${4:?hmm out dir}"
MAP_OUT="${5:?map out}"
THREADS="${6:-4}"
CONFIG="${7:?config}"
DO_DECOY="${8:-1}"

SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVS_ROOT="${ENVS_ROOT:-}"
ALLOW_MISSING="${ALLOW_MISSING:---allow-missing}"

act() {   # activate a pre-built env, or do nothing when snakemake manages conda
    [ -n "$ENVS_ROOT" ] || return 0
    set +u; source "$ENVS_ROOT/$1/bin/activate"; set -u
}

# Node-local scratch. SGE sets TMPDIR per job and wipes it on exit; the trap is
# for the local/interactive case and for signals SGE does not clean up after.
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/c71.${BATCH_ID}.XXXXXXXX")"
trap 'rm -rf "$SCRATCH"' EXIT INT TERM

mkdir -p "$HMM_OUT/decoy" "$(dirname "$MAP_OUT")"

echo "[search_batch] batch $BATCH_ID  scratch=$SCRATCH  threads=$THREADS" >&2
df -h "$SCRATCH" 2>/dev/null | tail -1 >&2 || true

# --- 1. concatenate + rename, straight to scratch ---------------------------
(
  act py
  python "$SCRIPTS/batch_faa.py" \
      --batch "$BATCH_TSV" --batch-id "$BATCH_ID" \
      --out-faa "$SCRATCH/batch.faa" --out-map "$MAP_OUT" $ALLOW_MISSING

  if [ "$DO_DECOY" = "1" ] && [ -s "$SCRATCH/batch.faa" ]; then
      python "$SCRIPTS/make_decoys.py" \
          --in-faa "$SCRATCH/batch.faa" --out-faa "$SCRATCH/decoy.faa"
  fi

  # profile name -> threshold flags, resolved from config once
  python - "$CONFIG" > "$SCRATCH/profiles.tsv" <<'PY'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1]))
for name, d in sorted(cfg["profiles"].items()):
    t = str(d["threshold"])
    flags = f"--{t}" if t.startswith("cut_") else f"-T {t} --domT {t} --incT {t} --incdomT {t}"
    print(f"{name}\t{flags}")
PY
)

# --- 2. hmmsearch, target and decoy -----------------------------------------
act hmmer

run_one() {   # run_one <hmm> <seqdb> <domtblout> <tblout> <flags...>
    local hmm="$1" db="$2" dom="$3" tbl="$4"; shift 4
    if [ ! -s "$db" ]; then
        printf '#\n' > "$dom"; printf '#\n' > "$tbl"
        return 0
    fi
    hmmsearch --cpu "$THREADS" --noali "$@" \
        --domtblout "$dom" --tblout "$tbl" -o /dev/null "$hmm" "$db"
}

while IFS=$'\t' read -r NAME FLAGS; do
    [ -n "$NAME" ] || continue
    HMM="$PROFILES_DIR/$NAME.hmm"
    [ -s "$HMM" ] || { echo "[search_batch] missing $HMM" >&2; exit 1; }

    # shellcheck disable=SC2086
    run_one "$HMM" "$SCRATCH/batch.faa" \
            "$HMM_OUT/batch_${BATCH_ID}.${NAME}.domtblout" \
            "$HMM_OUT/batch_${BATCH_ID}.${NAME}.tblout" $FLAGS

    if [ "$DO_DECOY" = "1" ]; then
        # shellcheck disable=SC2086
        run_one "$HMM" "$SCRATCH/decoy.faa" \
                "$HMM_OUT/decoy/batch_${BATCH_ID}.${NAME}.domtblout" \
                "$HMM_OUT/decoy/batch_${BATCH_ID}.${NAME}.tblout" $FLAGS
    fi
done < "$SCRATCH/profiles.tsv"

echo "[search_batch] batch $BATCH_ID done" >&2
