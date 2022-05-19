# 1. FlexddG
Requirements:
- pyxmolpp2 
- Rosetta
- slurm

[FlexddG](https://pubs.acs.org/doi/10.1021/acs.jpcb.7b11367) is a method within the Rosetta macromolecular modeling suite, that samples conformational diversity using “backrub” to generate an ensemble of models and then applies torsion minimization, side chain repacking, and averaging across this ensemble to estimate interface ΔΔG values.  
The present analysis is based on FlexddG [tutorial](https://github.com/Kortemme-Lab/flex_ddG_tutorial) and consists of saturation mutational scan and pair mutations scan.

## 1.1 Saturation analysis

During this step all residues from specifyed range in MP3 are mutated to each possible aminoacid and respective ΔΔG is estimated. Scripts for this step are located at 'scripts/saturation_analysis'.

1. Navigate to 'scripts/saturation_analysis'
```sh
cd scripts/saturation_analysis
```
2. Scipts will operate upon 'template' folder. The sample is provided for MP3/delta+ variant complex saturation mutational scan:
2.1  Chek the content of 'input/pdb' directory. File `chains_to_move.txt` contains the name of the MP3 chain (as in .pdb file), which will be mutated. Structure `delta_plus+mp3_minimized.pdb` is a structure prepared at the step №0. Files `pdb2rosetta.resmap.json` and `rosetta2pdb.resmap.json` contains information to link rosetta and pdb structure numbering. These are applicable for all MP3-RBD variants since the number of residues is fixed
2.2 Enter path to rosetta_scripts.mpi in `run_saturation.py`. You can also change default FlexddG settings (maximum cpus, backrub steps, covergence threshols etc.)
2.3 This project uses slurm to run computation tasks. You can change job parameters in `flexddg.sh`
3. In `create_flex_ddg_run_dir.py` check 'path_to_pdb' variable, which specifyes path to minimized.pdb. Execute the script. It will create subfolders for each residue in MP3:
```sh
python create_flex_ddg_run_dir.py
```
4. Run analysis in all subdirectories:
```sh
sbatch run_sh_in_subdirs.sh
```
5. When FlexddG run finishes, execute `analyze_flexddg.py`. It will create 'analysis_output' folders in all residue subdirectories with .csv tables containing results summary:
```sh
python analyze_flexddg.py
```
6. Tables produced with `analyze_flexddg.py` contain information with many parameters. In order to filter information for plotting use `summarize_ddG.py`, which is located above 'saturation_analysis' dir. This script works as command-lone tool, help is available: 
```sh
python summarize_ddG.py -h
```
Example command:
```sh
python summarize_ddG.py -i ./saturation_analysis -s ddg -b 35000 -n 35 -f fa_talaris2014-gam -o mutation_data_delta_plus_MP3
```
7. When results are extracted you can plot heatmap with `plot_saturation_mutation_results.py`. Provide path to .csv data inside script. You can try script on sample data provided in 'tables' dir (e.x. mutation_data_delta_plus_MP3.csv). This script works as a command-line tool. The example command:
```sh
python plot_saturation_mutation_results.py -i ../tables/mutation_data_delta_plus_MP3.csv -mutant delta+ -o ../plots
```
8. You can analyse discrete mutations and assemble a barplot of contributions to rosetta scores with `plot_score_contributions.py`. You can use sample data in 'tables' dir (output_struct_scores_delta_plus_37Asp.csv). This script works as a command-line tool. Note: this script requires path to .csv at 'analysis_output' subdirectories from Flexddg run. Use data with the suffix 'struct_scores_results.csv'. The sample output_struct_scores_delta_plus_37Asp.csv is provided in 'tables'. The example command:
```sh
python plot_score_contributions.py -i ../tables/output_struct_scores_delta_plus_37Asp.csv -mutant delta+ -start_res D -pos 37 -mut_res R -b 35000
```

**Obtained results:**
After performing saturation mutational scan, heatmaps with FlexddG scores for each mutation can be plotted. Mutations with scores below -1 are considered stabilizing and are highlighted in green. Substitutions to the same residue are highlighted in greay dashed boxes. Below are results for delta+ and omicron complexes. It is evident, that for omicron complex much less stabilizing mutations are found, which is consistent with MM-GBSA results shiwing stabilization of MP3/omicron complex (therefore not many stabilizing mutations are possible).

<p align="center">
  <img src="/results_plots/heatmap_delta+.png" width="700">
</p>

<p align="center">
  <img src="/results_plots/heatmap_omicron.png" width="700">
</p>

In case of MP3/delta+ a fovourable D37R mutation is found. This utation to polar residue may lead to formation of new polar interactions. Patricular contributions to FlexddG score for this mutation can be visualized as bar plot. In this case, polar Δhbond contributions seem not to be major, however a detailed investigation with MD and MM_GBSA decomposition is needed to check the effect of this mutation.

<p align="center">
  <img src="/results_plots/score_contributions_MP3_D37R.png" width="500">
</p>


## 1.2 Pair mutations analysis

In order to investigate pairs of mutation in MP3 execute scripts in 'pair_mutations'. The sample is also provided for delta+ variant complex. Data at 'inputs' dir is organased in same manner as for saturation mutational scan. 
  
Examine `run_mutation.py`. It differs from such script in 'saturation_analysis' in a way which `residue_combinations_to_mutate` variable is specifyed. Here a list of pair mutations to perform is provided as follows: [Chain, tuple(PDB residue numbers), tuple(new residues in single letter notation)].

0. Navigate to:
```sh
cd scripts/pair_mutations
``` 
1. After setting desired parameters and rosetta path in run_mutation.py, execute:
```sh
sbatch flexddg.sh
```
2. Once analysis is filished:
```sh
python analyze_flex_ddG.py
```
3. Draw a barplot of acquired ΔΔG values with `plot_pair_mutation_results.py`. This script accepts path to -struct_scores_results.csv and filters it prior to plotting. You can use sample data from 'tables' dir (delta_plus_MP3_pair_mutations.csv). The example command of the command-line tools script:
```sh
python plot_pair_mutation_results.py -i ../tables/delta_plus_MP3_pair_mutations.csv -p pair_mutations/inputs/pdb_deltap/delta_p+lcb3_minimized.pdb -mutant delta+ -s D37R -b 35000 -n 35 -f fa_talaris2014-gam
```
Note that this script operates upon assumtion of 1 mutation being fixed across all pairs.

**Obtained results:**
In order to find MP3 variants with two mutations a pair mutational scan is conducted. Pairs consist of D37R substitution and each stabilixing mutation found for MP3/delta+ complex earlier. Visualixation of FlexddG scores as box plot shows that all combinations are stabilizing for delta+ complex. the best result corresponds to T10W-D37R combination as it has one of the lowest scores and low std. All combinations of mutations are a destabilizing for omicron complex and are rather illustrative here.

<p align="center">
  <img src="/results_plots/delta+_pair_mutations.png" width="600">
</p>

<p align="center">
  <img src="/results_plots/omicron_pair_mutations.png" width="600">
</p>
