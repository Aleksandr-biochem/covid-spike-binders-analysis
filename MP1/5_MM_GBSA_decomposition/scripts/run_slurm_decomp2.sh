#!/bin/bash

#SBATCH --job-name="mm_gbsa_wt_o"
#SBATCH --ntasks=20           # Количество MPI процессов
#SBATCH --partition=cpu
#SBATCH -w, --nodelist=bionmr-mom-002,bionmr-mom-003

# We need amber 
source /home/oleg/amber/amber.sh
# We need python enviroment
source /home/xenia/cov2/venv/bin/activate

typeset -i i END
TRJ_FOLDER="../../2_MD_Amber/sample/" # Trajectory folder (i.e. where 1_build resides)
export STRAINS=$(echo wt+mp1 alpha+mp1 delta+mp1 delta_plus+mp1 omicron+mp1)

for complex in $STRAINS; do 
	cd "${TRJ_FOLDER}/${complex}/8_mmpbsa_decomposition/igb2_salt150"
	mpirun -np 10 --oversubscribe MMPBSA.py.MPI -O -i igb2_salt150_decomp.in -o complex_MMPBSA.dat -sp ../../0_prepare/box.prmtop -cp ../com_igb2.top -rp ../rec_igb2.top -lp ../lig_igb2.top -y ../mmpbsa.nc
	
	cd "${TRJ_FOLDER}/${complex}/8_mmpbsa_decomposition/igb8_salt150"
	mpirun -np 10 --oversubscribe MMPBSA.py.MPI -O -i igb8_salt150_decomp.in -o complex_MMPBSA.dat -sp ../../0_prepare/box.prmtop -cp ../com_igb8.top -rp ../rec_igb8.top -lp ../lig_igb8.top -y ../mmpbsa.nc
done;
