from MMPBSA_mods import API as MMPBSA_API
from collections import OrderedDict
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import os
import os.path
from pyxmolpp2 import PdbFile, mName, rId

# a function to convert data in _MMPBSA_info to pickle for more rapid parsing
def pickle_mmgbsa_data(path_to_info_file, path_to_output_dir=".", outname="MMPBSA_data.pickle"):
    wd = os.getcwd()
    os.chdir(os.path.dirname(path_to_info_file))
    data = MMPBSA_API.load_mmpbsa_info("_MMPBSA_info")
    os.chdir(wd)
    os.makedirs(path_to_output_dir, exist_ok=True)
    with open(os.path.join(path_to_output_dir, outname), 'wb') as f:
        pickle.dump(data, f)

# function to load pickle data
def load_pickle_data(pcike_file="MMPBSA_data.pickle"):
    with open(pcike_file, 'rb') as f:
        data = pickle.load(f)
    return data

# define folder-mutants names to be analyzed 
mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

# define gb decomp runs to be analyzed
runs = ["igb8_salt150_decomp", "igb2_salt150_decomp"]

# only frames starting from first_frmae are analyzed to omit the structure equibration span
first_frame = 500

os.makedirs('../decomposition_output', exist_ok=True)

# analyze each mutant and each run
for compl in mutants:
     print(f"Analysing {compl}")
     for run in runs:
        print(f"Analysing run {run}")

        # define path to info file
        trj_dir = f"../../2_MD_Amber/{compl}"
        decomposition_dir = f"{trj_dir}/8_mmpbsa_decomp/{run}"

        # dump data into pickle if no pickle is found
        if not os.path.isfile(f"{decomposition_dir}/MMPBSA_data_{compl}_{run}.pickle"):
            print("Dumping data to pickle")
            pickle_mmgbsa_data(path_to_info_file=f"{decomposition_dir}/_MMPBSA_info",
                               path_to_output_dir=decomposition_dir,
                               outname=f"MMPBSA_data_{compl}_{run}.pickle")

        # create output dataframe
        out_decomp_df = pd.DataFrame()

        ### read reference structure
        path_to_reference = f"{trj_dir}/0_prepare/protein_named.pdb"
        reference = PdbFile(path_to_reference).frames()[0]
        mp = reference.molecules.filter(mName == "A")
        rbd = reference.molecules.filter(mName == "B")
        rids_mp = [residue.id.serial for residue in mp.residues]
        rids_rbd = [residue.id.serial for residue in rbd.residues]

        ### set first rid in case end-to-end numbering (relevant for mmgbsa analysis)
        mp_start_rid = rids_mp[0]
        rbd_start_rid = rids_mp[-1] + 1

        ### load decomp_data
        decomp_data = load_pickle_data(pcike_file=f"{decomposition_dir}/MMPBSA_data_{compl}_{run}.pickle")

        ### read total energy data
        receptor_residue_pairs = np.array([*decomp_data['decomp']["gb"]['receptor']["TDC"]])
        ligand_residue_pairs = np.array([*decomp_data['decomp']["gb"]['ligand']["TDC"]])
        complex_residue_pairs = np.array([*decomp_data['decomp']["gb"]['complex']["TDC"]])

        ### get interacting residues mp+rbd
        mask_intra_residues = np.isin(complex_residue_pairs, receptor_residue_pairs) | np.isin(complex_residue_pairs, ligand_residue_pairs)
        complex_residue_pairs = complex_residue_pairs[~(mask_intra_residues)]
        mask_pairs = [True if int(residue_pairs.split("-")[0]) <= rbd_start_rid else False for residue_pairs in complex_residue_pairs]
        complex_residue_pairs = complex_residue_pairs[mask_pairs]

        ### retaining only those pairs, whose trajectory-average energy exceeded threshold
        energy_threshold = 1 # kcal/mol
        temp = {}
        for residue_pair in tqdm(complex_residue_pairs, desc=f"{run} {compl}"):
        # for residue_pair in complex_residue_pairs:
            rid_first, rid_second = residue_pair.split("-")

            ### decomp_matrix is pseido-symmetry. sum of the per-residue contributions will equal the total binding free energy
            ### Gohlke, H., C. Kiel, and D. A. Case. 2003 https://pubmed.ncbi.nlm.nih.gov/12850155/
            total_energy = decomp_data['decomp']["gb"]['complex']["TDC"][f"{rid_first}-{rid_second}"]['tot'].copy() + \
                           decomp_data['decomp']["gb"]['complex']["TDC"][f"{rid_second}-{rid_first}"]['tot'].copy()

            rid_mp = rids_mp[int(rid_first) - mp_start_rid]
            rid_rbd = rids_rbd[int(rid_second) - rbd_start_rid]

            rname_mp = mp.residues.filter(rId == rid_mp)[0].name
            rname_rbd = rbd.residues.filter(rId == rid_rbd)[0].name

            total_energy_avg, total_energy_std = total_energy[first_frame:].mean(), total_energy[first_frame:].std()

            ### retaining only those pairs, whose trajectory-average energy exceeded 1 kcal/mol 
            if abs(total_energy_avg) > 1.0:
                temp = {"rId_mp": rid_mp, "rName_mp": rname_mp, "rId_rbd": rid_rbd, "rName_rbd": rname_rbd, "total_energy_avg": total_energy_avg, "total_energy_std": total_energy_std}
                out_decomp_df = pd.concat([out_decomp_df, pd.DataFrame(temp, index=[0])])

        # save data
        out_decomp_df.to_csv(f"../decomposition_output/{compl}_{run}_decomp.csv", index=False)





