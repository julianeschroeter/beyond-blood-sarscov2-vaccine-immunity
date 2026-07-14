#!/bin/bash

#SBATCH --output=/opt/mesh/raasay/chris/aim-cv/output/output.%N_%A_%a.out
#SBATCH --job-name=AIM-CV
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --nodelist=mingulay,tiree,eorsa,harris,canna,eriskay,pabbay,taransay,eigg
#SBATCH --array=1-400 # number of tasks
#SBATCH --nice=10000

cd /opt/mesh/raasay/chris/aim-cv

Rscript aim_cv_hpc.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT
