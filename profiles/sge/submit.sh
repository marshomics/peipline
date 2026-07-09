#!/usr/bin/env bash
# ============================================================================
# qsub wrapper for snakemake's cluster-generic executor.
#
#   submit.sh RULE THREADS MEM_MB H_RT QUEUE <jobscript>
#
# Three things this exists to get right, none of which a bare `qsub ...` line
# in the profile can express.
#
# 1. h_vmem is a CONSUMABLE resource on this cluster, which in SGE means it is
#    requested PER SLOT. `-pe parallel 16 -l h_vmem=64G` asks for 1024 GB, not
#    64 GB. Every such job either waits forever or is killed by the reaper for
#    the wrong reason. We take the job's total memory and divide by slots.
#
# 2. standard.q has h_rt=24:00:00 (soft 23:55). A job requesting more is never
#    schedulable there. Rules that genuinely need days (iqtree, hmmalign) are
#    routed to long.q; anything mis-routed gets promoted with a warning rather
#    than silently pending.
#
# 3. `-pe parallel 1` is legal but pointlessly narrows the eligible hosts on
#    some SGE builds. Single-threaded jobs get no -pe at all.
#
# Probe facts this encodes (sge_probe_report, 2026-07-06):
#   SGE 8.1.9; PE `parallel` = $pe_slots (single host, correct for our threaded
#   jobs); PE `openmpi` = $round_robin (multi-host, wrong for every rule here).
#   h_vmem: requestable YES, consumable YES.  s_vmem/mem_free: not consumable.
#   Queue h_rt: standard.q 24h, long.q 672h, cryo-em.q 672h, test.q unlimited.
#
# Env knobs:
#   SGE_EXTRA   extra qsub args, appended verbatim
#   SGE_MEM_FREE=1  also request -l mf=<per-slot>, to avoid landing on a node
#                   that is already swapping (mem_free is not consumable, so
#                   this is advisory scheduling only)
# ============================================================================
set -euo pipefail

RULE="${1:?rule}"; THREADS="${2:?threads}"; MEM_MB="${3:?mem_mb}"
H_RT="${4:?h_rt}"; QUEUE="${5:?queue}"
shift 5
JOBSCRIPT="${*: -1}"

mkdir -p logs/sge

# --- slots ------------------------------------------------------------------
SLOTS="$THREADS"
[[ "$SLOTS" =~ ^[0-9]+$ ]] || SLOTS=1
[ "$SLOTS" -lt 1 ] && SLOTS=1

# Smallest exec host has 32 cores (node529). Asking for more than that on a
# $pe_slots PE means the job can only ever run on the handful of 64/128/256
# core hosts, and will queue behind them.
MAX_SLOTS="${SGE_MAX_SLOTS:-32}"
if [ "$SLOTS" -gt "$MAX_SLOTS" ]; then
    echo "submit.sh: $RULE asked for $SLOTS slots; capping at $MAX_SLOTS (smallest node)" >&2
    SLOTS="$MAX_SLOTS"
fi

# --- memory: total -> per slot ----------------------------------------------
[[ "$MEM_MB" =~ ^[0-9]+$ ]] || MEM_MB=8000
PER_SLOT_MB=$(( (MEM_MB + SLOTS - 1) / SLOTS ))
[ "$PER_SLOT_MB" -lt 512 ] && PER_SLOT_MB=512

# --- walltime vs queue ------------------------------------------------------
to_sec() {   # HH:MM:SS -> seconds
    local h m s
    IFS=: read -r h m s <<< "$1"
    echo $(( 10#${h:-0} * 3600 + 10#${m:-0} * 60 + 10#${s:-0} ))
}
REQ=$(to_sec "$H_RT")

case "$QUEUE" in
    standard.q) CAP=$(to_sec "23:55:00") ;;   # s_rt, not h_rt: SIGUSR1 lands first
    long.q|cryo-em.q) CAP=$(to_sec "670:00:00") ;;
    test.q) CAP=$(( 10**9 )) ;;
    *) CAP=$(( 10**9 )) ;;
esac

if [ "$REQ" -gt "$CAP" ]; then
    if [ "$QUEUE" = "standard.q" ]; then
        echo "submit.sh: $RULE wants $H_RT but standard.q caps at 23:55:00; promoting to long.q" >&2
        QUEUE="${SGE_LONG_QUEUE:-long.q}"
    else
        echo "submit.sh: $RULE wants $H_RT which exceeds $QUEUE; submitting anyway" >&2
    fi
fi

# --- build ------------------------------------------------------------------
ARGS=(-cwd -V -terse
      -N "c71.${RULE}"
      -q "$QUEUE"
      -l h_rt="$H_RT"
      -l h_vmem="${PER_SLOT_MB}M"
      -o logs/sge/ -e logs/sge/)

# `parallel` is the $pe_slots PE: all slots on one host. Never `openmpi`
# ($round_robin), which would scatter the slots and break every threaded tool
# in this pipeline.
[ "$SLOTS" -gt 1 ] && ARGS+=(-pe "${SGE_PE:-parallel}" "$SLOTS")

[ "${SGE_MEM_FREE:-0}" = "1" ] && ARGS+=(-l mf="${PER_SLOT_MB}M")

# shellcheck disable=SC2206
[ -n "${SGE_EXTRA:-}" ] && ARGS+=(${SGE_EXTRA})

if [ "${SGE_DRY_RUN:-0}" = "1" ]; then
    echo "qsub ${ARGS[*]} $JOBSCRIPT" >&2
    echo "999999"
    exit 0
fi

exec qsub "${ARGS[@]}" "$JOBSCRIPT"
