#!/usr/bin/env bash
# Generic Linux runner (any server / GitHub Actions). Resumable: re-run the SAME command after any stop.
# Usage: bash scripts/run_linux.sh [conditions=core8|all|C029,C030,...] [rep_from=1] [rep_to=200] [workers=$(nproc)]
#        env: OUT (results_core8)  R_PINNED_LIB (pinned package library)  TIME_BUDGET_MIN (Inf)
#             SHARD / NSHARDS (split work across machines or CI jobs; default 1/1)
set -euo pipefail
cd "$(dirname "$0")/../sim"
C=${1:-core8}; FROM=${2:-1}; TO=${3:-200}; W=${4:-$(nproc)}
[ "$C" = core8 ] && C=C029,C030,C031,C032,C047,C050,C053,C056
OUT=${OUT:-results_core8}; SHARD=${SHARD:-1}; NSHARDS=${NSHARDS:-1}; BUDGET=${TIME_BUDGET_MIN:-Inf}
[ -n "${R_PINNED_LIB:-}" ] && export R_LIBS_USER="$R_PINNED_LIB"
mkdir -p "$OUT/logs"; rm -f "$OUT"/C*/rep_*.rds.tmp_* 2>/dev/null || true
NW=$((W * NSHARDS))
for i in $(seq 1 "$W"); do
  WK=$(( (SHARD - 1) * W + i ))
  IRTREE_COND="$C" IRTREE_REP_FROM="$FROM" IRTREE_REP_TO="$TO" IRTREE_OUT="$OUT" \
  IRTREE_WORKER="$WK" IRTREE_NWORKERS="$NW" IRTREE_SEED_SCHEME=main IRTREE_ORACLE_DRAWS=1 IRTREE_REQUIRE_FROZEN=1 \
  IRTREE_TIME_BUDGET_MIN="$BUDGET" IRTREE_EXPECT_REP_MIN=11 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  Rscript main_sim.R >> "$OUT/logs/worker_${WK}.log" 2>&1 &
done
wait
echo "saved: $(find "$OUT" -name 'rep_*.rds' | wc -l)"
