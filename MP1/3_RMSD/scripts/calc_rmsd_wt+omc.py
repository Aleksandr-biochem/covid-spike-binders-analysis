import os
import sys
sys.path.append('../../handling/LCB1/md-utils/md_utils')
from tqdm import tqdm
import argparse

import numpy as np
import pandas as pd

from frame.select import get_sec_str_residues_predicate
from pyxmolpp2 import PdbFile, Trajectory, AmberNetCDF, calc_rmsd, mName, aName
from pyxmolpp2.pipe import AssembleQuaternaryStructure

def rmsd_calc_wt_omcr(mutant="", trj_path=".", outpth="."):
  """
  Calculate RMSD between Ca atoms of residues on interface
  for wt and omicron variants.

  param: mutant - for which mutant RMSD is calculated (e.g. alpha+lcb1, omicron+lcb1)
  param: trj_path - Path to directory with MD analysis of trajectories
  param: outpth - Path to output directory where result is stored

  """

  path_to_trj_pdb = os.path.join(trj_path, f"{mutant}+mp1_2/")
  path_to_trj_dir = os.path.join(trj_path, f"{mutant}+mp1/")
  path_to_trj_ref = os.path.join(path_to_trj_pdb, "0_prepare/protein_named.pdb")

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
  lcb_sec_str_ca_predicate = get_sec_str_residues_predicate(frame=trj_ref, molnames=["A"]) & (aName == "CA")
  rbd_sec_str_ca_predicate = get_sec_str_residues_predicate(frame=trj_ref, molnames=["B"]) & (aName == "CA")
  all_sec_str_ca_predicate = get_sec_str_residues_predicate(frame=trj_ref, molnames=["A", "B"]) & (aName == "CA")

  #  reference atoms for alignment
  ref_lcb_sec_str_ca = trj_ref.atoms.filter(lcb_sec_str_ca_predicate)
  ref_rbd_sec_str_ca = trj_ref.atoms.filter(rbd_sec_str_ca_predicate)
  ref_all_sec_str_ca = trj_ref.atoms.filter(all_sec_str_ca_predicate)

  #  create data variables
  rmsd_lcb_sec_str_ca = np.zeros(int(traj.size / stride))
  rmsd_rbd_sec_str_ca = np.zeros(int(traj.size / stride))
  rmsd_all_sec_str_ca = np.zeros(int(traj.size / stride))

  #  run through trajectory
  for frame in tqdm(traj[::stride] | AssembleQuaternaryStructure(of=(mName.is_in("A", "B")),
                                                                 by=all_sec_str_ca_predicate,
                                                                 reference=trj_ref)):
      if frame.index == 0:
          # create variables for selections from atoms of frame
          frame_lcb_sec_str_ca = frame.atoms.filter(lcb_sec_str_ca_predicate)
          frame_rbd_sec_str_ca = frame.atoms.filter(rbd_sec_str_ca_predicate)
          frame_all_sec_str_ca = frame.atoms.filter(all_sec_str_ca_predicate)

      # rmsd for sec.str. CA atoms of lcb
      alignment = frame_lcb_sec_str_ca.alignment_to(ref_lcb_sec_str_ca)
      crd = frame_lcb_sec_str_ca.coords.values.copy()
      crd = crd @ alignment.matrix3d().T + alignment.vector3d().values
      rmsd_lcb_sec_str_ca[int(frame.index / stride)] = calc_rmsd(ref_lcb_sec_str_ca.coords.values,
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
  #rmsd_fnout = f"../../2_MD_Amber/sample/RMSD/{mutant}+mp1_rmsd_ca_ss.csv"
  rmsd_fnout = os.path.join(outpth, f"{mutant}+mp1_rmsd_ca_ss.csv")
  col_names = [r"time_ns",
               r"sec.str. C$\rm\alpha$ lcb [$\rm\AA$]",
               r"sec.str. C$\rm\alpha$ rbd [$\rm\AA$]",
               r"sec.str. C$\rm\alpha$ [$\rm\AA$]"]

  time = np.linspace(trj_start, trj_length, int(traj.size / stride))
  pd.DataFrame(np.vstack((time, rmsd_lcb_sec_str_ca, rmsd_rbd_sec_str_ca, rmsd_all_sec_str_ca)).T,
               columns=col_names).to_csv(rmsd_fnout,
                                         index=False)

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Calculate rmsd of wt and omicron complexes per each ns of simulated trajectories and return in a csv format:')
    
    parser.add_argument('-i', nargs=1, help='Path to folder containing subfolders with MD simulations results', required=True, default="../../2_MD_Amber/sample/")

    parser.add_argument('-m', nargs=1, help='Strain to analyze', required=True)
    
    parser.add_argument('-o', nargs=1, help='Path to output folder', required=True, default="./../tables/")
    
    args = parser.parse_args()

    #os.makedirs("../tables/RMSD", exist_ok=True)

    rmsd_calc_wt_omcr(mutant=args.m[0], trj_path=args.i[0], outpth=args.o[0])
