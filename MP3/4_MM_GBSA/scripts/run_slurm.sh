#!/bin/bash

#SBATCH --job-name="mm_gbsa"
#SBATCH --ntasks=20                             # MPI process number
#SBATCH --partition=cpu
#SBATCH --nodelist=bionmr-mom-001

# state amber source
# note! if MPI was used with task preparations then the amber build with MPI is required
source /path/to/amber/amber.sh
# state python environment source
source ~/path/to/python/venv/bin/activate


typeset -i i END
# Folder containing all trajectories
TRJ_FOLDER="../../2_MD_Amber"

# an array of mutants
declare -a mutants=("wt+mp3" "alpha+mp3" "delta+mp3" "delta_plus+mp3" "omicron+mp3" "delta_p+mp3_d37r" "delta_p+mp3_d37r_t10w")

cd "${TRJ_FOLDER}"

for complex in ${mutants[@]}; do 
	# folders are specifyed according to names provided in prep_MMPBSA_jobs.py
	cd "${complex}/7_mmpbsa/igb2_salt150"
	mpirun -np 10 MMPBSA.py.MPI -O -i ../script150_igb2.in -o complex_MMPBSA.dat -sp ../../0_prepare/box.prmtop -cp ../com_igb2.top -rp ../rec_igb2.top -lp ../lig_igb2.top -y ../mmpbsa.nc

	cd "../igb8_salt150"
	mpirun -np 10 MMPBSA.py.MPI -O -i ../script150_igb8.in -o complex_MMPBSA.dat -sp ../../0_prepare/box.prmtop -cp ../com_igb8.top -rp ../rec_igb8.top -lp ../lig_igb8.top -y ../mmpbsa.nc
	cd "../../.."
done
