#!/bin/bash

#SBATCH --job-name="flex_ddg"
#SBATCH --ntasks=1                              # MPI pricess number
#SBATCH --cpus-per-task=20                      # CPU number
#SBATCH --partition=cpu

python run_saturation.py
