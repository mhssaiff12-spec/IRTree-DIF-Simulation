#!/usr/bin/env bash
# Launches W independent, resumable workers for the main simulation.
# Usage: bash run_parallel.sh <W> <conditions: all | C029,C030,...> <rep_from> <rep_to> [out_dir]
# Each worker takes every W-th (condition, replication) task; finished tasks (rep_XXXX.rds) are skipped,
# so re-running the SAME command after an interruption resumes exactly where it stopped.
W=${1:-2}; CONDS=${2:-all}; FROM=${3:-1}; TO=${4:-100}; OUT=${5:-results}
mkdir -p "$OUT/logs"
for i in $(seq 1 "$W"); do
  IRTREE_COND="$CONDS" IRTREE_REP_FROM="$FROM" IRTREE_REP_TO="$TO" IRTREE_OUT="$OUT" \
  IRTREE_WORKER="$i" IRTREE_NWORKERS="$W" IRTREE_SEED_SCHEME=main IRTREE_ORACLE_DRAWS=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  nohup Rscript main_sim.R > "$OUT/logs/worker_$i.log" 2>&1 &
done
wait
echo "all workers finished"; ls "$OUT"/C*/rep_*.rds 2>/dev/null | wc -l
