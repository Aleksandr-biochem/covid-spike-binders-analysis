import os
import pickle
from tqdm import tqdm
import argparse

import pandas as pd
import numpy as np

from MMPBSA_mods import API as MMPBSA_API
from pyxmolpp2 import PdbFile, mName, rId

def pickle_mmgbsa_data(path_to_info_file, path_to_output_dir="./", cmplx=None, run=None):
    wd = os.getcwd()
    os.chdir(os.path.dirname(path_to_info_file))
    data = MMPBSA_API.load_mmpbsa_info("_MMPBSA_info")
    os.chdir(wd)
    os.makedirs(path_to_output_dir, exist_ok=True)
    with open(os.path.join(path_to_output_dir, f'MMPBSA_data_{cmplx}+mp1_{run}.pickle'), 'wb') as f:
        pickle.dump(data, f)

def load_pickle_data(pcike_file="MMPBSA_data.pickle"):
    with open(pcike_file, 'rb') as f:
        data = pickle.load(f)
    return data

def sum_gbsa_decomp(trj_path=".", outpth="summary_gbsa_decomp", answ="y"):
    complexes = ["wt", "alpha", "delta", "delta_plus", "omicron"]
  # define gb runs to be analyzed 
    runs = ['igb2_salt150', 'igb8_salt150']

    for cmplx in complexes:
        print(f"Start PBSA analysis for {cmplx}")
        for run in runs:
            print(f"Analysing run for {run} model")
            decomp_path = os.path.join(trj_path, f"{cmplx}+mp1/8_mmpbsa_decomposition/{run}")
            #try:
            #	os.path.exists(decomp_path)
            #except FileNotFoundError:
            #	print("Directory is not found. No decomposition folder.")

            # read reference
            if cmplx in ["wt", "omicron"]:
            	path_to_reference = f"../../2_MD_Amber/sample/{cmplx}+mp1_2/0_prepare/protein_named.pdb"
            else:
            	path_to_reference = f"../../2_MD_Amber/sample/{cmplx}+mp1/0_prepare/protein_named.pdb"
            
            reference = PdbFile(path_to_reference).frames()[0]
            lcb = reference.molecules.filter(mName == "A")
            rbd = reference.molecules.filter(mName == "B")
            rids_lcb = [residue.id.serial for residue in lcb.residues]
            rids_rbd = [residue.id.serial for residue in rbd.residues]

            ### set first rid in case end-to-end numbering (relevant for mmgbsa analysis)
            lcb_start_rid = rids_lcb[0]
            rbd_start_rid = rids_lcb[-1] + 1

            # get data as a nested dictionary
            print('Checking compressed data existence....')
            
            if os.path.exists(os.path.join(decomp_path, f"MMPBSA_data_{cmplx}+mp1_{run}.pickle")):
            	print("Data is compressed. Loading already pickled data...")
            	if answ == 'y':
            		pickle_mmgbsa_data(path_to_info_file=f"{decomp_path}/_MMPBSA_info", path_to_output_dir=decomp_path, cmplx=cmplx, run=run)
            	#else:data = load_pickle_data(pcike_file=f"MMPBSA_data_{cmplx}_{run}.pickle")
            	#data = load_pickle_data(pcike_file=f"MMPBSA_data_{cmplx}_{run}.pickle")
            
            else:
            	print("Data is not compressed. Compressing...")
            	pickle_mmgbsa_data(path_to_info_file=f"{decomp_path}/_MMPBSA_info", path_to_output_dir=decomp_path, cmplx=cmplx, run=run)
            	print("Compression ended")

            print("Loading data...")
            data = load_pickle_data(pcike_file=f"{decomp_path}/MMPBSA_data_{cmplx}+mp1_{run}.pickle")

            print('Finished loading data')

            receptor_residue_pairs = np.array([*data['decomp']["gb"]['receptor']["TDC"]])
            ligand_residue_pairs = np.array([*data['decomp']["gb"]['ligand']["TDC"]])
            complex_residue_pairs = np.array([*data['decomp']["gb"]['complex']["TDC"]])

            mask_intra_residues = np.isin(complex_residue_pairs, receptor_residue_pairs) | np.isin(complex_residue_pairs, ligand_residue_pairs)

            # filter intra-ligand and intra-receptor interactions
            complex_residue_pairs = complex_residue_pairs[~(mask_intra_residues)]

            mask_pairs = [True if int(residue_pairs.split("-")[0]) <= rbd_start_rid else False for residue_pairs in complex_residue_pairs]
            complex_residue_pairs = complex_residue_pairs[mask_pairs]

            out_decomp_df = pd.DataFrame()
            temp = {}
            for residue_pair in tqdm(complex_residue_pairs, desc=f"{run} model {cmplx}+mp1"):
            	rid_first, rid_second = residue_pair.split('-')
            	total_energy = data['decomp']["gb"]['complex']["TDC"][f"{rid_first}-{rid_second}"]['tot'].copy() + data['decomp']["gb"]['complex']["TDC"][f"{rid_second}-{rid_first}"]['tot'].copy()
            	rid_lcb = rids_lcb[int(rid_first) - lcb_start_rid]
            	rid_rbd = rids_rbd[int(rid_second) - rbd_start_rid]

            	rname_lcb = lcb.residues.filter(rId == rid_lcb)[0].name
            	rname_rbd = rbd.residues.filter(rId == rid_rbd)[0].name

            	total_energy_avg, total_energy_std = total_energy[500:].mean(), total_energy[500:].std()

            	if abs(total_energy_avg) > 1.0:
                      temp = {"rId_lcb": rid_lcb, "rName_lcb": rname_lcb, "rId_rbd": rid_rbd, "rName_rbd": rname_rbd, "total_energy_avg": total_energy_avg, "total_energy_std": total_energy_std}
                      out_decomp_df = pd.concat([out_decomp_df, pd.DataFrame(temp, index=[0])])
            # assign brief model name
            model = run.split('_')[0]
            if not os.path.exists(outpth):
              os.mkdir(outpth)
            out_decomp_df.to_csv(f"{outpth}/{cmplx}+mp1_{model}_sum_decomp_GBSA.TDC.significant.tsv", sep='\t', index = False)
            print(f"Decomposition data is saved in {outpth}!")

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Summarize MM-GBSA decomposition data into tsv format')
    
    parser.add_argument('-i', nargs=1, help='Path to folder containing MM-GBSA decomposition info file', required=True, default="../../2_MD_Amber/sample")

    parser.add_argument('-y', nargs=1, help='Overwrite or not already existing pickles: y/n', required=True, default="y")

    parser.add_argument('-o', nargs=1, help='Path to output folder', required=True, default="../tables/")
    
    args = parser.parse_args()

    sum_gbsa_decomp(trj_path=args.i[0], outpth=args.o[0], answ=args.y[0])

