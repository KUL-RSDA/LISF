#!/bin/bash -l
#SBATCH -t 00:59:59
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH -o compile_log.txt
#SBATCH -e compile_out.txt

cd $SLURM_SUBMIT_DIR
source ../load_modules

./compile -j $SLURM_CPUS_PER_TASK

