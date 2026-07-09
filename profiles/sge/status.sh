#!/usr/bin/env bash
# SGE/UGE job-status probe for snakemake's cluster-generic executor.
# Must print exactly one of: running | success | failed
set -u
JOBID="$1"

# 1. Still in the queue/running?
if qstat -j "$JOBID" &> /dev/null; then
    echo "running"
    exit 0
fi

# 2. Finished -- ask the accounting database for the exit status.
#    qacct can lag behind qstat by a few seconds on busy clusters, so retry.
for _ in 1 2 3 4 5 6; do
    OUT=$(qacct -j "$JOBID" 2> /dev/null) && break
    sleep 5
done

if [ -z "${OUT:-}" ]; then
    # Never appeared in accounting. Treat as running rather than failing the
    # DAG; snakemake's latency-wait will catch a genuinely dead job.
    echo "running"
    exit 0
fi

EXIT=$(echo "$OUT" | awk '$1=="exit_status"{print $2; exit}')
FAILED=$(echo "$OUT" | awk '$1=="failed"{print $2; exit}')

if [ "${EXIT:-1}" = "0" ] && [ "${FAILED:-1}" = "0" ]; then
    echo "success"
else
    echo "failed"
fi
