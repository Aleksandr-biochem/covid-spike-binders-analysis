# 2. MD in Amber
Requirements:
- Amber20
- amber_runner
- Propka

0. Before executing scripts, make sure that Amber is added to paths:
```sh
source path/to/amber20
```

1. Go to scripts directory to launch the simulations.
```sh
cd scripts
```

2. The script `init_cov2_1.py` prepares MD jobs by creating a directory for each complex and saving job parameters as state.dill. The script can initiate multiple jobs at once, specify `mutants` variable in the script with the names of MP-RDB complexes.
```sh
python3 init_cov2_1.py
```

**Some default simulation parameters**:
- 1500 ns simulation
- ff14SB force field
- tip4pew water model
- pH 7.5 (protonation state generated with Propka)
- temperature at production 298K

2. `run_cov2.py` will launch the simulations for the complexes specifyed in `mutants` variable inside the script:
```sh
python3 run_cov2.py
```
3. When jrajectory writing is finished it can be converted to .pdb format with 1 ns stride with a scpecifyed stride using `get_pdb_trajectory.py` (complexes are specifyed in `mutants` variable inside the script):
```sh
python3 get_pdb_trajectory.py
```
Note that there are `trj_length` and `stride` parameters in this script (default 1499 and 1000 respectively). Change them if you work with trajectory longer or shorter than 1500 ns or if you want a different stride.
