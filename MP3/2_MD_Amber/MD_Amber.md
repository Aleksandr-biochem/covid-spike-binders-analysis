# 2. MD in Amber
Requirements:
- Amber20
- [amber-runner](https://github.com/sizmailov/amber-runner) (0.0.8)
- [pyxmolpp2](https://github.com/sizmailov/pyxmolpp2) (1.6.0)
- [Propka](https://github.com/jensengroup/propka) (3.1)
- [md-utils](https://github.com/OOLebedenko/md-utils)

0. Before executing scripts, make sure that Amber20 is added to paths:
```sh
source path/to/amber20.sh
```
1. The script `init_cov2.py` prepares MD jobs by creating a directory for each complex and saving job parameters as `state.dill`. The script can initiate multiple jobs at once. If neccesary, edit `mutants` variable in the script with the names of MP/RDB complexes to simulate.
```sh
python init_cov2.py
```

**Some default simulation parameters**:
- 1500 ns simulation
- ff14SB force field
- tip4pew water model
- pH 7.5 (protonation state generated with `Propka`)
- temperature at production 298K

2. `run_cov2.py` will launch the simulations for the complexes specified in `mutants` variable inside the script. Before execution specify paths to python env and Amber20 inside the script.
```sh
python run_cov2.py
```
3. When trajectory writing is finished it can be converted to .pdb format with a scpecifyed stride (default 1 ns or 1000 frames stride, can be changed insode stript in `stride` variable) using `get_pdb_trajectory.py` (complexes are specifyed in `mutants` variable inside the script):
```sh
python get_pdb_trajectory.py
```
Note that there is also `trj_length` parameter in this script (default 1499, takes `stride` into account). Change it if you work with trajectory longer or shorter than 1500 ns or if you want a different stride.
