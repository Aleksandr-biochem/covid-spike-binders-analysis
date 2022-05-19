# 1. FlexddG
Requirements:
- pyxmolpp2 
- Rosetta
- slurm

[FlexddG](https://pubs.acs.org/doi/10.1021/acs.jpcb.7b11367) is a method within the Rosetta macromolecular modeling suite, that samples conformational diversity using “backrub” to generate an ensemble of models and then applies torsion minimization, side chain repacking, and averaging across this ensemble to estimate interface ΔΔG values.  
The present analysis is based on FlexddG [tutorial](https://github.com/Kortemme-Lab/flex_ddG_tutorial) and consists of saturation mutational scan and pair mutations scan.

## 1.1 Saturation analysis

During this step all residues from specifyed range in MP1 are mutated to each possible amino acid and respective ΔΔG is estimated. Scripts for this step are located in 'scripts/saturation_analysis'.

1. Navigate to 'scripts/saturation_analysis'
```sh
cd scripts/saturation_analysis
```
2. Scipts will operate upon 'template' folder. The sample is provided for MP1/delta+ variant complex saturation mutational scan:
2.1  Chek the content of 'input/pdb' directory. File `chains_to_move.txt` contains the name of the MP1 chain (as in .pdb file), which will be mutated. Structure `delta_plus+mp1_minimized.pdb` is a structure prepared at the step №0. Files `pdb2rosetta.resmap.json` and `rosetta2pdb.resmap.json` contains information to link rosetta and pdb structure numbering. These are applicable for all MP1-RBD variants since the number of residues is fixed
2.2 Enter path to rosetta_scripts.mpi in `run_saturation.py`. You can also change default FlexddG settings (maximum cpus, backrub steps, covergence threshols etc.)
2.3 This project uses slurm to run computation tasks. You can change job parameters in `flexddg.sh`
3. In `create_flex_ddg_run_dir.py` check 'path_to_pdb' variable, which specifyes path to minimized.pdb. Execute the script. It will create subfolders for each residue in MP3:
```sh
python create_flex_ddg_run_dir.py
```
4. Run analysis in all subdirectories:
```sh
sbatch run_subdirs.sh
```
5. When FlexddG run finishes, execute `analyze_flexddg.py`. It will create 'analysis_output' folders in all residue subdirectories with .csv tables containing results summary:
```sh
python run_flexddG.py
```
6. Tables produced with `run_flexddG.py` contain information with many parameters. In order to filter information and plot the results use`flexddG_sum.py`, which is located above 'saturation_analysis' dir. This script works as command-lone tool, help is available: 
```sh
python flexddG_sum.py -h
```
Example command:
```sh
python3 flexddG_sum.py -i LCB1 -s deltaplus -f fa_talaris2014 -o delta_plus_flexddg_sum.csv
```
7. When results are extracted you can plot heatmap with `heatmap.py`. Provide path to .csv data inside script. You can try script on sample data provided in 'tables' dir (e.x. delta_plus_flexddg_sum.csv). This script works as a command-line tool. The example command:
```sh
python3 heatmap.py -i ../tables/delta_plus_flexddg_sum.csv -s delta+ -o delta_plus_flexddG.png
```
*The plot is saved in the current directory* 

# Obtained results:

All resutls are available in `result_plots` directory.

After performing saturation mutational scan, heatmaps with FlexddG scores for each mutation were plotted. The Y axis represents primary MP1 amino acids and on the X axis mutated residues are placed. Each score stands for the difference in dinding energy between primary and mutated amino acid in MP1. Mutations with scores below -1 are considered to be stabilizing and highlighted with blue frames. Substitutions to the same residue are highlighted in grey frames. Below are results for delta+ complex. It is obvious, that practically none of mutations tends to stabilize the complex. The majority of them lead to the decrease of binding energy, in other words - destabilize the structure. 

<p align="center">
  <img src="plots/delta_plus_flexddG.png" width="700">
</p>

In case of MP1/delta+ a fovourable E3W mutation is found. This mutation leads to the itroduction of space-sufficient hydrophobic aromatic group, which could disrupt the MP1 structure itself and induce complex destabilization. Thus none of mutations can be called favourable. 
