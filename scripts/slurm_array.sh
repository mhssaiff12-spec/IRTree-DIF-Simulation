#!/bin/bash
# SLURM alternative: one array task per condition (56 tasks), each running a replication range.
#SBATCH --job-name=irtree_dif
#SBATCH --array=1-56
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=72:00:00
COND=$(printf "C%03d" "$SLURM_ARRAY_TASK_ID")
export IRTREE_COND=$COND IRTREE_REP_FROM=${REP_FROM:-1} IRTREE_REP_TO=${REP_TO:-100} IRTREE_OUT=${OUT:-results}
export IRTREE_SEED_SCHEME=main IRTREE_ORACLE_DRAWS=1 OMP_NUM_THREADS=1
Rscript main_sim.R
