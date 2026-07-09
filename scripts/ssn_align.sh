#!/usr/bin/env bash
# All-vs-all for the sequence similarity network.
#
# Nodes are 100%-identity representatives, then MMseqs2 representatives if the
# set is still above ssn.max_nodes. An SSN with 10^5 nodes is unreadable and the
# all-vs-all is O(n^2); EFI-EST itself dereplicates for the same reason.
#
# Edges come from DIAMOND blastp with --very-sensitive. EFI-EST's edge weight is
# -log10(E-value), computed in ssn.py from column 11 of this m8.
set -euo pipefail

FAA="$1"; M8="$2"; NODES="$3"; CFG="$4"; TMP="$5"; THREADS="$6"
mkdir -p "$TMP" "$(dirname "$M8")"

python - "$CFG" > "$TMP/params.sh" <<'PY'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))["ssn"]
print(f'EVALUE={c["evalue_max"]}')
print(f'SENS="{c["diamond_sensitivity"]}"')
print(f'MAXN={c["max_nodes"]}')
r = yaml.safe_load(open(sys.argv[1]))["redundancy"]
print(f'MINID={r["cluster_min_seq_id"]}')
print(f'COV={r["cluster_coverage"]}')
PY
source "$TMP/params.sh"

seqkit rmdup -s -o "$TMP/derep.faa" "$FAA" 2> "$TMP/rmdup.log"
N=$(grep -c '^>' "$TMP/derep.faa" || true)
echo "[ssn_align] $N unique sequences after 100% dereplication" >&2

if [ "$N" -gt "$MAXN" ]; then
    echo "[ssn_align] above max_nodes=$MAXN; clustering at $MINID identity" >&2
    mmseqs easy-cluster "$TMP/derep.faa" "$TMP/clu" "$TMP/mm" \
        --min-seq-id "$MINID" -c "$COV" --cov-mode 0 --threads "$THREADS" -v 1
    cp "$TMP/clu_rep_seq.fasta" "$NODES"
else
    cp "$TMP/derep.faa" "$NODES"
fi
echo "[ssn_align] $(grep -c '^>' "$NODES") network nodes" >&2

diamond makedb --in "$NODES" --db "$TMP/nodes.dmnd" --threads "$THREADS" --quiet
diamond blastp \
    --query "$NODES" --db "$TMP/nodes.dmnd" \
    --out "$M8" --outfmt 6 qseqid sseqid pident length mismatch gapopen \
        qstart qend sstart send evalue bitscore qlen slen \
    --evalue "$EVALUE" $SENS \
    --max-target-seqs 0 --threads "$THREADS" --quiet

echo "[ssn_align] $(wc -l < "$M8") raw alignments -> $M8" >&2
