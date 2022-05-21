import sys 
# Add path to md-utils repository
sys.path.append('path_to_md-utils')
from md_utils.frame.select import get_sec_str_residues_predicate
from pyxmolpp2.pipe import AssembleQuaternaryStructure, Align
from pyxmolpp2 import PdbFile, Trajectory, AmberNetCDF, mName, aName
import os
from tqdm import tqdm

# chose mutant trajectories to process
mutants = ["wt+mp1", "alpha+mp1", "delta+mp1", "delta_plus+mp1", "omicron+mp1"]

for mutant in mutants:
    print(f'Processing {mutant}')

    if mutant in ["omicron+mp1", "wt+mp1"]:
        path_to_trj_dir = f"../sample/{mutant}_2"
    else:
        path_to_trj_dir = f"../sample/{mutant}" # path to trajectory folder
    path_to_trj_ref = os.path.join(path_to_trj_dir, "0_prepare/protein_named.pdb")

    # set trj parameters
    trj_start = 1
    trj_length = 1499
    stride = 1000

    # read trajectory
    trj_ref = PdbFile(path_to_trj_ref).frames()[0]
    traj = Trajectory(trj_ref)
    for ind in tqdm(range(trj_start, trj_length + 1), desc="traj_reading"):
        fname = "{pattern}.{filetype}".format(pattern="run%05d", filetype="nc")
        traj.extend(AmberNetCDF(os.path.join(os.path.join(path_to_trj_dir, "6_run"), fname % (ind))))

    # set predicate sec.str CA atoms
    all_sec_str_ca_predicate = get_sec_str_residues_predicate(frame=trj_ref, molnames=["A", "B"]) & (aName == "CA")

    # save trajectory in pdb
    frame_atoms = None
    with open(f"../{mutant}.pdb", "w") as fout:
        #  run through trajectory
        for frame in tqdm(traj[::stride]
                          | AssembleQuaternaryStructure(of=(mName.is_in("A", "B")),
                                                        by=all_sec_str_ca_predicate,
                                                        reference=trj_ref)
                          | Align(by=all_sec_str_ca_predicate, reference=trj_ref)
                          ):

            if frame_atoms is None:
                frame_atoms = frame.atoms
            frame.to_pdb(fout)
            fout.write("ENDMDL\n")
