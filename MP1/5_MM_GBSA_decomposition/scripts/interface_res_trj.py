import sys
import os
from tqdm import tqdm
from collections import defaultdict
from pyxmolpp2.pipe import AssembleQuaternaryStructure, Align
from pyxmolpp2 import PdbFile, Trajectory, AmberNetCDF, mName, AtomPredicate

sys.path.append('../../../handling/LCB1/md-utils/md_utils')
from frame.extract import extract_residues_on_interface
from frame.select import get_sec_str_residues_predicate

# only wt residues are enough for analysis
complexes = ["wt", "alpha", "delta", "delta_plus","omicron"]

# specify distance threshold for residue extraction
d_cutoff = 6

# set trj parameters
trj_start = 1
trj_length = 1499
stride = 1000

intersected_rIDs_wt = None

for cmplx in complexes:
	print(f'Start {cmplx}+mp1 analysis')
	path_to_trj_dir = f"../../2_MD_Amber/sample/{cmplx}+mp1/"

	# get initial and renumbered reference structures
	path_to_trj_ref = os.path.join(path_to_trj_dir, "0_prepare/box.pdb")
	if cmplx in ["omicron", "wt"]:
		path_to_trj_dir = f"../../2_MD_Amber/sample/{cmplx}+mp1_2/"
		path_to_xray_ref = os.path.join(path_to_trj_dir, "0_prepare/protein_named.pdb")	
	path_to_xray_ref = os.path.join(path_to_trj_dir, "0_prepare/protein_named.pdb")
	xray_ref = PdbFile(path_to_xray_ref).frames()[0]
	trj_ref = PdbFile(path_to_trj_ref).frames()[0]


	chain_list = [] # generate predicates to connect between initial and renumbered structures
	for m_dest, m_src in zip(trj_ref.molecules, xray_ref.molecules):
		chain_list.append(m_src.name)
		m_dest.name = m_src.name
	chain_set = set(chain_list)
	#predicates = [] # list of lists hard to parse
	predicates = {} # allow to save additional inf per residue
	for chain in chain_list: # per molecule
		# avoid hydrogen atoms
		predicates[chain] = AtomPredicate(lambda a: not a.name.startswith("H")) & mName.is_in(chain)

	traj = Trajectory(trj_ref)
	traj.extend(AmberNetCDF(f"../../2_MD_Amber/sample/{cmplx}+mp1/7_mmpbsa/mmpbsa.nc"))
	
	# set predicate sec.str CA atoms -- we can use mmpbsa.nc file, the full frames information
	# all_sec_str_ca_predicate = get_sec_str_residues_predicate(frame=trj_ref, molnames=["A", "B"])
	# traj = Trajectory(xray_ref)
    # for ind in tqdm(range(traj_start, traj_length), desc=f"{cmplx} traj reading"):
    #    traj.extend(AmberNetCDF(f"/home/xenia/cov2/trj/LCB1/{cmplx}/6_run/{ind:05d}.nc"))

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


	# analyze each frame from trajectory assembled from run

	# for frame in tqdm(traj[::stride] | AssembleQuaternaryStructure(of=(mName.is_in("A", "B")),
	# 					                                           by=all_sec_str_ca_predicate,
	# 					                                           reference=trj_ref), desc="Trajectory analysis"):

	# 	interface_residues_partner_A, interface_residues_partner_B = extract_residues_on_interface( 
	# 	partner_A=frame.molecules.filter(mName == "A"), 
	# 	partner_B=frame.molecules.filter(mName == "B"),
	# 	cutoff=distance_threshold)

	# 	for res in interface_residues_partner_A: 
	# 		interface_res_A[int(str(res.id))] += 1

	# 	for res in interface_residues_partner_B: 
	# 		interface_res_B[int(str(res.id))] += 1

	# write output residues with occurancies
	with open(f"../tables/{cmplx}+mp1_distance_{d_cutoff}.txt", "w") as f:
		f.write("res_id\toccurancy\n")
		for res in interface_res_A:
			occurancy = interface_res_A[res] * 100 / n_frames
			f.write(f"{res}\t{occurancy}\n")
		for res in interface_res_B:
			occurancy = interface_res_B[res] * 100 / n_frames
			f.write(f"{res}\t{occurancy}\n")

	
	residues_for_decompostion = [int(i) for i in interface_res_A.keys()]
	residues_for_decompostion.extend(list(interface_res_B.keys()))
	residues_for_decompostion.sort()
	# print(residues_for_decompostion)

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

	intersected_rIDs = ', '.join(intersected_rIDs)
	if cmplx == "wt":
		intersected_rIDs_wt = intersected_rIDs
	print(intersected_rIDs)


	# generate scripts for MM-GBSA with decomposition
	path_to_trj_out = f"../../2_MD_Amber/sample/{cmplx}+mp1/"
	for model in ['2', '8']:
		if not (os.path.exists(os.path.join(path_to_trj_out, f"8_mmpbsa_decomposition/igb{model}_salt150/"))):
			os.makedirs(os.path.join(path_to_trj_out, f"8_mmpbsa_decomposition/igb{model}_salt150/"))

		#os.chdir(os.path.join(path_to_trj_out, f"8_mmpbsa_decomposition/igb{model}_salt150/"))
		with open(os.path.join(path_to_trj_out, f"8_mmpbsa_decomposition/igb{model}_salt150/igb{model}_salt150_decomp.in"), 'w') as f:
			f.write("Per-residue GB decomposition\n" \
				"&general\n" \
				"endframe=99999999, verbose=1,\n" \
				"/\n" \
				"&gb\n" \
				f"igb={model}, saltcon=0.150\n" \
				"/\n" \
				"&decomp\n" \
				"csv_format=1, idecomp=4, dec_verbose=3,\n" \
				f"print_res=\"{intersected_rIDs_wt}\"\n" \
				"/\n")
