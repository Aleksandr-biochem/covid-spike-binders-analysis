import sys
import os
from tqdm import tqdm
from collections import defaultdict
from pyxmolpp2.pipe import AssembleQuaternaryStructure, Align
from pyxmolpp2 import PdbFile, Trajectory, AmberNetCDF, mName, AtomPredicate

sys.path.append('../../../md_utils')
from md_utils.frame.extract import extract_residues_on_interface
from md_utils.frame.select import get_sec_str_residues_predicate


# define folder-mutants names to be analyzed 
# note that only wt residues are enough for analysis
mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

# specify distance cutoff for residue extraction
d_cutoff = 6

# set trj parameters
trj_start = 1
trj_length = 1499
stride = 1000

os.makedirs('../interface_residues', exist_ok=True)

for cmplx in mutants:
	print(f'Analysing {cmplx}')
	# define path to trajectory folder
	path_to_trj_dir = f"../../2_MD_Amber/{cmplx}/"

	# get initial and renumbered reference structures from trajectory folder
	path_to_trj_ref = os.path.join(path_to_trj_dir, "0_prepare/box.pdb")
	path_to_xray_ref = os.path.join(path_to_trj_dir, "0_prepare/protein_named.pdb")
	xray_ref = PdbFile(path_to_xray_ref).frames()[0]
	trj_ref = PdbFile(path_to_trj_ref).frames()[0]

	chain_list = [] # generate predicates to connect between initial and renumbered structures
	for m_dest, m_src in zip(trj_ref.molecules, xray_ref.molecules):
		chain_list.append(m_src.name)
		m_dest.name = m_src.name
	chain_set = set(chain_list)
	predicates = {} # allow to save additional inf per residue
	for chain in chain_list: # per molecule
		# avoid hydrogen atoms
		predicates[chain] = AtomPredicate(lambda a: not a.name.startswith("H")) & mName.is_in(chain)

	# in this case the trajectory assempled for mm-gbsa is used
	traj = Trajectory(trj_ref)
	traj.extend(AmberNetCDF(f"../../2_MD_Amber/{cmplx}/7_mmpbsa/mmpbsa.nc"))
	
	# dictinoaries to collect residues with treir occurancies
	interface_res_A, interface_res_B = defaultdict(int), defaultdict(int)

	# analyze mmpbsa.nc
	residues_of_interest = []
	index = traj[0].index
	n_frames = 0
	for frame in tqdm(traj, desc="Extraction of residues on interface: "):
		if frame.index == index:
			# create variables for selections from atoms of frame
			frame_ats = frame.molecules.atoms
			frame_by_chain = {}
			for chain in chain_list:
				frame_by_chain[chain] = frame_ats.filter(predicates[chain]).residues

		# get interface residues
		interface_residues_partner_A, interface_residues_partner_B = extract_residues_on_interface( 
		partner_A=frame_by_chain[chain_list[0]], 
		partner_B=frame_by_chain[chain_list[1]],
		cutoff=d_cutoff)

		# count if residue is found in frame A mol
		for res in interface_residues_partner_A: 
			interface_res_A[int(str(res.id))] += 1

		# count if residue is found in frame B mol
		for res in interface_residues_partner_B: 
			interface_res_B[int(str(res.id))] += 1

		n_frames += 1


	# write output residues with occurancies
	with open(f"../interface_residues/{cmplx}_interface_res_distance_{d_cutoff}.txt", "w") as f:
		f.write("res_id\toccurancy\n")
		for res in interface_res_A:
			occurancy = interface_res_A[res] * 100 / n_frames
			f.write(f"{res}\t{occurancy}\n")
		for res in interface_res_B:
			occurancy = interface_res_B[res] * 100 / n_frames
			f.write(f"{res}\t{occurancy}\n")

	# list all found residues and generate intervals for decomposition input
	residues_for_decompostion = [int(i) for i in interface_res_A.keys()]
	residues_for_decompostion.extend(list(interface_res_B.keys()))
	residues_for_decompostion.sort()

	# find intersections and return intervals
	intersected_rIDs = []
	curr_id, prev_id = residues_for_decompostion[0], residues_for_decompostion[0]
	for rID in residues_for_decompostion[1:]:
		if (prev_id + 1) < rID:
			if curr_id != prev_id:
				intersected_rIDs.append(f"{curr_id}-{prev_id}")
			else:
				intersected_rIDs.append(str(curr_id))
			curr_id = rID

		prev_id = rID

	if curr_id != prev_id:
		intersected_rIDs.append(f"{curr_id}-{prev_id}")
	else:
		intersected_rIDs.append(str(curr_id))

	intersected_rIDs = ','.join(intersected_rIDs)
	print(intersected_rIDs) # print intervals for decomposition input
