import pandas as pd
import numpy as np

import os
import argparse
import pickle
from tqdm import tqdm
import re

from MMPBSA_mods import API as MMPBSA_API
from pyxmolpp2 import PdbFile, mName, rId

def load_pickle_data(pcike_file="MMPBSA_data.pickle"):
	with open(pcike_file, 'rb') as f:
		data = pickle.load(f)
	return data

def polar_contribution_data(path_to_trj=".", outpth=".", mp_residue="", rbd_residue=""):
	"""
    Extract incofrmation of polar contribution
     into total dG of MM-GBSA decomposition results.
	The information is extracted for the pair of residues of interes.

    param: path_to_trj - path to the folder with MD trj, including mm_gbsa_decomposition subfolder
    param: outpth - path to the directory where ouput table will be generated
    param: mp_residue - MP residue of interest, e.g. SER29
    param: rbd_residue - RBD residue of interest, e.g. ASP417

    """

	mutants = ["alpha", "delta", "delta_plus", "omicron"]

	for mut in mutants: 
		decomp_path = os.path.join(path_to_trj, f"{mut}+mp1/8_mmpbsa_decomposition/igb8_salt150")
		data = load_pickle_data(pcike_file=os.path.join(decomp_path, f"MMPBSA_data_{mut}+mp1_igb8_salt150.pickle"))
		if mut == "omicron":
			path_to_reference = f"../../2_MD_Amber/sample/{mut}+mp1_2/0_prepare/protein_named.pdb"
		else:
			path_to_reference = f"../../2_MD_Amber/sample/{mut}+mp1/0_prepare/protein_named.pdb"
		
		reference = PdbFile(path_to_reference).frames()[0]
		mp = reference.molecules.filter(mName == "A")
		rbd = reference.molecules.filter(mName == "B")
		rids_mp = [residue.id.serial for residue in mp.residues]
		rids_rbd = [residue.id.serial for residue in rbd.residues]

		### set first rid in case end-to-end numbering (relevant for mmgbsa analysis)
		mp_start_rid = rids_mp[0]
		rbd_start_rid = rids_mp[-1] + 1

		receptor_residue_pairs = np.array([*data['decomp']["gb"]['receptor']["TDC"]])
		ligand_residue_pairs = np.array([*data['decomp']["gb"]['ligand']["TDC"]])
		complex_residue_pairs = np.array([*data['decomp']["gb"]['complex']["TDC"]])

		mask_intra_residues = np.isin(complex_residue_pairs, receptor_residue_pairs) | np.isin(complex_residue_pairs, ligand_residue_pairs)

		# filter intra-ligand and intra-receptor interactions
		complex_residue_pairs = complex_residue_pairs[~(mask_intra_residues)]

		mask_pairs = [True if int(residue_pairs.split("-")[0]) <= rbd_start_rid else False for residue_pairs in complex_residue_pairs]
		complex_residue_pairs = complex_residue_pairs[mask_pairs]

		for residue_pair in tqdm(complex_residue_pairs, desc=f"igb8_model {mut}"):
			rid_first, rid_second = residue_pair.split('-')
			total_energy = data['decomp']["gb"]['complex']["TDC"][f"{rid_first}-{rid_second}"]['tot'].copy() + data['decomp']["gb"]['complex']["TDC"][f"{rid_second}-{rid_first}"]['tot'].copy()
			polar_energy = data['decomp']["gb"]['complex']["TDC"][f"{rid_first}-{rid_second}"]['pol'].copy() + data['decomp']["gb"]['complex']["TDC"][f"{rid_second}-{rid_first}"]['pol'].copy()
			rid_mp = rids_mp[int(rid_first) - mp_start_rid]
			rid_rbd = rids_rbd[int(rid_second) - rbd_start_rid]

			rname_mp = mp.residues.filter(rId == rid_mp)[0].name
			rname_rbd = rbd.residues.filter(rId == rid_rbd)[0].name

			total_energy_avg = total_energy[500:].mean()

			out_polar_contr = pd.DataFrame(columns=['time', 'total_energy', 'polar_energy'])
			#abs(total_energy_avg) > 1.0 and 
			if (rid_mp == int(mp_residue[3:])) and (rid_rbd == int(rbd_residue[3:])):
				out_polar_contr = pd.DataFrame({
		           	'time': range(500,1500),
		           	'total_energy': total_energy[500:],
		           	'polar_energy': polar_energy[500:]
		           	})
				break
			else:
				continue
		outd_path = outpth
		out_polar_contr.to_csv(os.path.join(outd_path, f"{mut}+mp1_polar_{mp_residue}-{rbd_residue}.tsv"), sep='\t', index=False)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Summarize MM-GBSA decomposition data for residues pf interest in tsv format')
    
    parser.add_argument('-i', nargs=1, help='Path to MD trajectory, where decomposition subfolder is', required=True, default="../../2_MD_Amber/sample")

    parser.add_argument('-r1', nargs=1, help='MP residue in the format AA№, where AA - three letter amino acid code', type=str, required=True)

    parser.add_argument('-r2', nargs=1, help='RBD residue in the format AA№, where AA - three letter amino acid code', type=str, required=True)

    parser.add_argument('-o', nargs=1, help='Path to output folder', required=True, default="../tables/")
    
    args = parser.parse_args()
    mp = str(args.r1[0])
    rbd = str(args.r2[0])
    polar_contribution_data(path_to_trj=args.i[0], outpth=args.o[0], mp_residue=mp, rbd_residue=rbd)

