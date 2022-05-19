import sys
sys.path.append('../../../md_utils')
from md_utils.frame.select import get_sec_str_residues_predicate
from pyxmolpp2 import PdbFile, Trajectory, AmberNetCDF, calc_rmsd, mName, aName, rId
from pyxmolpp2.pipe import AssembleQuaternaryStructure
from tqdm import tqdm
import os
import numpy as np
import pandas as pd

# specify mutants
mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

os.makedirs('../rmsd_output', exist_ok=True)

for mutant in mutants:
    print(f'Processing {mutant}')

    path_to_trj_dir = f"../../2_MD_Amber/{mutant}/"
    path_to_trj_ref = os.path.join(path_to_trj_dir, "0_prepare/protein_named.pdb")

    # set trj length parameters
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
    mp_sec_str_ca_predicate = get_sec_str_residues_predicate(frame=trj_ref, molnames=["A"]) & (aName == "CA")
    rbd_sec_str_ca_predicate = get_sec_str_residues_predicate(frame=trj_ref, molnames=["B"]) & (aName == "CA")
    all_sec_str_ca_predicate = get_sec_str_residues_predicate(frame=trj_ref, molnames=["A", "B"]) & (aName == "CA")

    #  reference atoms for alignment
    ref_mp_sec_str_ca = trj_ref.atoms.filter(mp_sec_str_ca_predicate)
    ref_mp_interacting_helices_sec_str_ca = trj_ref.atoms.filter(mp_sec_str_ca_predicate & (rId.is_in(set([_ for _ in range(1, 42)]))))
    ref_mp_distanced_helix_sec_str_ca = trj_ref.atoms.filter(mp_sec_str_ca_predicate & (rId.is_in(set([_ for _ in range(42, 65)]))))
    ref_rbd_sec_str_ca = trj_ref.atoms.filter(rbd_sec_str_ca_predicate)
    ref_all_sec_str_ca = trj_ref.atoms.filter(all_sec_str_ca_predicate)

    #  create data variables
    rmsd_mp_sec_str_ca = np.zeros(int(traj.size / stride))
    rmsd_mp_interacting_helices_sec_str_ca = np.zeros(int(traj.size / stride))
    rmsd_mp_distanced_helix_sec_str_ca = np.zeros(int(traj.size / stride))
    rmsd_rbd_sec_str_ca = np.zeros(int(traj.size / stride))
    rmsd_all_sec_str_ca = np.zeros(int(traj.size / stride))

    #  run through trajectory
    for frame in tqdm(traj[::stride] | AssembleQuaternaryStructure(of=(mName.is_in("A", "B")),
                                                                   by=all_sec_str_ca_predicate,
                                                                   reference=trj_ref)):
        if frame.index == 0:
            # create variables for selections from atoms of frame
            frame_mp_sec_str_ca = frame.atoms.filter(mp_sec_str_ca_predicate)
            frame_mp_interacting_helices_sec_str_ca = frame.atoms.filter(mp_sec_str_ca_predicate & (rId.is_in(set([_ for _ in range(1, 42)])))) # for MP3
            frame_mp_distanced_helix_sec_str_ca = frame.atoms.filter(mp_sec_str_ca_predicate & (rId.is_in(set([_ for _ in range(42, 65)])))) # for MP3
            frame_rbd_sec_str_ca = frame.atoms.filter(rbd_sec_str_ca_predicate)
            frame_all_sec_str_ca = frame.atoms.filter(all_sec_str_ca_predicate)

        # rmsd for sec.str. CA atoms of MP
        alignment = frame_mp_sec_str_ca.alignment_to(ref_mp_sec_str_ca)
        crd = frame_mp_sec_str_ca.coords.values.copy()
        crd = crd @ alignment.matrix3d().T + alignment.vector3d().values
        rmsd_mp_sec_str_ca[int(frame.index / stride)] = calc_rmsd(ref_mp_sec_str_ca.coords.values,
                                                                   crd
                                                                   )

        # rmsd for sec.str. CA atoms of MP3 helices interacting with RBD
        alignment = frame_mp_interacting_helices_sec_str_ca.alignment_to(ref_mp_interacting_helices_sec_str_ca)
        crd = frame_mp_interacting_helices_sec_str_ca.coords.values.copy()
        crd = crd @ alignment.matrix3d().T + alignment.vector3d().values
        rmsd_mp_interacting_helices_sec_str_ca[int(frame.index / stride)] = calc_rmsd(ref_mp_interacting_helices_sec_str_ca.coords.values,
                                                                   crd
                                                                   )

        # rmsd for sec.str. CA atoms of mp3 helix distanced from RBD
        alignment = frame_mp_distanced_helix_sec_str_ca.alignment_to(ref_mp_distanced_helix_sec_str_ca)
        crd = frame_mp_distanced_helix_sec_str_ca.coords.values.copy()
        crd = crd @ alignment.matrix3d().T + alignment.vector3d().values
        rmsd_mp_distanced_helix_sec_str_ca[int(frame.index / stride)] = calc_rmsd(ref_mp_distanced_helix_sec_str_ca.coords.values,
                                                                   crd
                                                                   )

        # rmsd for sec.str. CA atoms of Ig-like domain
        alignment = frame_rbd_sec_str_ca.alignment_to(ref_rbd_sec_str_ca)
        crd = frame_rbd_sec_str_ca.coords.values.copy()
        crd = crd @ alignment.matrix3d().T + alignment.vector3d().values
        rmsd_rbd_sec_str_ca[int(frame.index / stride)] = calc_rmsd(ref_rbd_sec_str_ca.coords.values,

                                                                   crd)

        # rmsd for all sec.str. CA atoms
        alignment = frame_all_sec_str_ca.alignment_to(ref_all_sec_str_ca)
        crd = frame_all_sec_str_ca.coords.values.copy()
        crd = crd @ alignment.matrix3d().T + alignment.vector3d().values
        rmsd_all_sec_str_ca[int(frame.index / stride)] = calc_rmsd(ref_all_sec_str_ca.coords.values,
                                                                   crd)

    # write RMSD to file
    col_names = [r"time_ns",
                 r"sec.str. C$\rm\alpha$ MP $\rm\AA$",
                 r"sec.str. C$\rm\alpha$ RDB $\rm\AA$",
                 r"sec.str. C$\rm\alpha$ complex $\rm\AA$",
                 r"sec.str. C$\rm\alpha$ MP interface helices $\rm\AA$",
                 r"sec.str. C$\rm\alpha$ MP distanced helix $\rm\AA$"]

    time = np.linspace(trj_start, trj_length, int(traj.size / stride))
    pd.DataFrame(np.vstack((time, rmsd_mp_sec_str_ca, rmsd_rbd_sec_str_ca, rmsd_all_sec_str_ca,
                            rmsd_mp_interacting_helices_sec_str_ca,rmsd_mp_distanced_helix_sec_str_ca)).T,
                 columns=col_names).to_csv(f"../rmsd_output/rmsd_ca_ss_{mutant}.csv", index=False)
