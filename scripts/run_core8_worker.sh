#!/usr/bin/env bash
# One resumable worker for the core-8 run (200 reps, new seeds, frozen code). Usage: bash run_core8_worker.sh <worker 1|2>
cd /home/claude/main_sim
IRTREE_COND=C029,C030,C031,C032,C047,C050,C053,C056 IRTREE_REP_FROM=1 IRTREE_REP_TO=200 IRTREE_OUT=results_core8 \
IRTREE_SEED_SCHEME=main IRTREE_ORACLE_DRAWS=1 IRTREE_WORKER=$1 IRTREE_NWORKERS=2 IRTREE_REQUIRE_FROZEN=1 \
IRTREE_TIME_BUDGET_MIN=112 IRTREE_EXPECT_REP_MIN=11 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
Rscript main_sim.R >> results_core8_worker$1.log 2>&1
echo "worker $1 exit $?; saved: $(ls results_core8/C*/rep_*.rds 2>/dev/null | wc -l)/1600"
