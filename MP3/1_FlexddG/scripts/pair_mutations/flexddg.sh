#!/bin/bash

#SBATCH --job-name="flex_ddg_LCB3"
#SBATCH --ntasks=1                              # MPI process number
#SBATCH --cpus-per-task=20                      # CPU number
#SBATCH --partition=cpu

python run_mutation.py
