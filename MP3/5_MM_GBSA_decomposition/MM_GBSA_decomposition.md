# 4. MM-GBSA binding free energy estimation
Requirements:
- Amber20
- mpi

Within this step MM-GBSA decomposition calculations are performed to investigate contributions of distinct residue pairs in MP3/RBD binding. 

- MM-GBSA_decomposition jobs are prepared at the MD folders and saved as state.dill objects at '8_mmpbsa_decomp' subdirs
- Two models are used in scripts: igb2 and igb8. To our experience, their results are quite consistent, so employng igb8 could be sufficient
- Salt concentration is set to 0.15M
- `script150_igb2_decomp.in` and `script150_igb8_decomp.in` are template scripts containing parameters for MM-GBSA runs. These are modifyed with `decomp` parameters during the job preparation. These can be modifyed with additional parameters if neccessary
- By default 1000 frames from trajectory are used with 1 ns stride starting from 500 ns of simulation. This span corresponds to an equilibrated complex configuration

0. A list of interface residues is required to run decomposition. This list is found using `analyze_interface_residues.py`. **Note:** this script operates on 7_mmpbsa/mmpbsa.nc trajectory assembled in step №4 for MM-GBSA. If you want to skip MM-GBSA calculations, at least execute `prep_MMPBSA_jobs.py` to create thess .nc files.
The script will create 'interface_residues' dir, write .txt files with interface residues and their trajectory occurancies in % and print spans of interfce residues ids for decomposition .in scripts:
```sh
python analyze_interface_residues.py
```
1.  Execute `prep_MMPBSA_decomp_jobs.py` after scpecifying the list of complexes for analysis in `mutants` variable. The list of interface residues to be analyzed is provided inside the script as `roi` variable. Note that this script assumes working with MPI, therefore mpi-amber assembly is needed:
```sh
python prep_MMPBSA_decomp_jobs.py
```
2. Run jobs by executing (don't forget to specify mutants):
```sh
sbatch run_slurm_decomposition.sh
```
3. Amber produces summary in .dat format, however for a detailed analysis _MMPBSA_info should be parsed using Amber's API. Execute `analyze_mmgbsa_decomp.py` after specifying `mutants` for analysis. Precomputed tables for MP3 are available in 'tables'.
Note that processing can take a while, because _MMPBSA_info is time-consuming. However, after first script run the parsed data is dumped into picke, so the next parsing will be more rapid. By default frames are analyzed starting from 500. It can be changed in `first_frame` variable inside the script.
```sh
python analyze_mmgbsa_decomp.py
```
4. In order to plot decomposition results as differential contact maps use `plot_decomp_result.py` providing the path to data folder and igb model to summarize:
```sh
python plot_decomp_result.py -i ../tables -igb igb2
```
The plot is saved to 'plots' dir. If you add new complex type to the analysis, modify `mutants` variable inside the script.

