# 4. MM-GBSA binding free energy estimation
Requirements:
- `Amber20` (build with MPI)
- `slurm` (20.11.8)

Within this step MM-GBSA calculations are performed to estimate binding free energy of MP3 complexes with different RDB variants. 

- MM-GBSA jobs are prepared at the MD folders (at `2_MD_Amber`) and saved as `state.dill` objects at '7_mmpbsa' subdirs
- Two models are used in scripts: igb2 and igb8. To our experience, their results are quite consistent, so employng igb8 could be sufficient. However, a comparison of two models can be helful.
- Salt concentration is set to 0.15M
- `script150_igb2.in` and `script150_igb8.in` scripts contain parameters for MM-GBSA runs. These can be modifyed with additional parameters, if neccessary
- By default 1000 frames from trajectory are used with 1 ns stride starting from 500 ns of simulation. This span correspond to an equilibrated complex configuration

1.  Execute `prep_MMPBSA_jobs.py` after scpecifying the list of complexes for analysis in `mutants` variable. Note that this script assumes working with MPI, therefore mpi-amber assembly is needed:
```sh
python prep_MMPBSA_jobs.py
```
2. Run jobs (don't forget to check `mutants`, and specify paths to amber and python env inside the script):
```sh
sbatch run_slurm.sh
```
3. Amber produces summary in .dat format, however for a detailed analysis `_MMPBSA_info` file should be parsed using Amber's API. Execute `analyze_gbsa.py` after specifying `mutants` for analysis. This will create .csv tables with data on complex, ligand and receptor total energies. The column 'delta' will contain binding free energies found as: `G(complex) - G(ligand) - G(receptor)`. Precomputed tables for MP3 are available in `tables`.
```sh
python analyze_gbsa.py
```
4. In order to plot MM-GBSA binding free energies execute `plot_gbsa.py` providing the path to data folder, igb model to summarize and `-s` state to plot (dG or ddG computed as difference from wild-type complex ΔG binding). You can also provide starting frame as `-f` option. By default it is 500 (se explanation in the beggining):
```sh
python plot_gbsa.py -h # to get help
python plot_gbsa.py -i ../tables -igb igb8 -f 500 -s dG
```
**Note:** If you add new complex type to the analysis, modify `mutants` and `labels` variables inside the script.

**Obtained results:**  
MM-GBSA results can be plotted as either ΔG or ΔΔG in relation to MP3/wild-type complex. The latter option allows more convinient comparison of two models. Binding free energy estimations shows that MP3 affinity changes differently across RBD variants (see ΔΔG bar charts below). Depending on the model chosen for MM-GBSA calculation (igb2 or 8) predictions for alpha and delta variants vary. However, a significant destabilization is observed for the complex with delta+ RDB (ΔΔG > 0) and stabilization is evident for omicron binding (ΔΔG < 0). 

These plots also show that mutant MP3 variants have improved binding to delta+ RBD compared to initial ligand. However, addition of second mutation (T10W) adds no significant improvement. This result is investigated further with decomposition analysis.

<p align="center">
  <img src="results_plots/ddG_mmgbsa_igb2.png" width="600">
</p>

<p align="center">
  <img src="results_plots/ddG_mmgbsa_igb8.png" width="600">
</p>
