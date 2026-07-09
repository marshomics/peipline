#!/usr/bin/env bash
# Driver for the C71 pipeline on an SGE/UGE cluster.
#
#   ./run.sh preflight # validate every input path, column, tip label and HMM
#   ./run.sh dry       # DAG only, touches nothing
#   ./run.sh local     # everything on the current node (testing / small inputs)
#   ./run.sh cluster   # submit to SGE (default)
#
# Run preflight first. It takes a minute and catches the mistakes that otherwise
# surface six hours into a 700-job array.
set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:-cluster}"
mkdir -p logs/sge
chmod +x profiles/sge/status.sh setup/*.sh 2> /dev/null || true

# If envs_root is set, the environments are pre-built and snakemake must NOT try
# to create or solve any of its own -- there is no internet on the compute nodes.
ENVS_ROOT=$(python3 -c "import yaml,sys; print(yaml.safe_load(open('config.yaml')).get('envs_root') or '')" 2>/dev/null || echo "")
if [ -n "$ENVS_ROOT" ]; then
    DEPLOY=()
    echo "[run.sh] pre-built environments: $ENVS_ROOT (snakemake conda disabled)"
else
    DEPLOY=(--software-deployment-method conda --conda-frontend mamba)
    echo "[run.sh] no envs_root set: snakemake will create conda envs (needs internet)"
fi

case "$MODE" in
  preflight)
    shift || true
    python3 scripts/preflight.py --config config.yaml "$@"
    ;;
  dry)
    snakemake -n -p --rerun-incomplete "${DEPLOY[@]}"
    ;;
  local)
    snakemake --cores "${CORES:-32}" "${DEPLOY[@]}" --rerun-incomplete --keep-going -p
    ;;
  cluster)
    # Run snakemake itself inside screen/tmux, or submit it as a long single job.
    snakemake --profile profiles/sge "${DEPLOY[@]}"
    ;;
  unlock)
    snakemake --unlock
    ;;
  dag)
    snakemake --rulegraph | dot -Tsvg > rulegraph.svg
    echo "wrote rulegraph.svg"
    ;;
  *)
    echo "usage: $0 {preflight|dry|local|cluster|unlock|dag}" >&2; exit 1
    ;;
esac
