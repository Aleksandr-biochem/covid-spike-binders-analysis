#!/bin/bash

#SBATCH --job-name="flex_ddg_S29"
#SBATCH --ntasks=1                              # Количество MPI процессов
#SBATCH --cpus-per-task=20                      # Требуемое кол-во CPU
#SBATCH --partition=cpu

python run_saturation.py
