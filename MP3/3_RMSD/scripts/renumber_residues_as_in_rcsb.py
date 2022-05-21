from pyxmolpp2 import PdbFile, rId, Frame, ResidueId, mName, rName
import os

# specify mutant complexes for renumbering
mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

# read rcsb reference structure
path_to_rcsb_ref = "../../0_prepare_structures/structures/7jzm.pdb"
rcsb_model = PdbFile(path_to_rcsb_ref).frames()[0]
rcsb_molnames = [mol.name for mol in rcsb_model.molecules.filter(mName.is_in("A", "B"))]

for mutant in mutants:
    # read ref trj pdb
    path_to_trj_ref = f"../../2_MD_Amber/{mutant}/0_prepare/protein.pdb"
    trj_ref = PdbFile(path_to_trj_ref).frames()[0]
    trj_ref_mols = [rcsb_model.molecules.filter(mName == mol_name)[0] for mol_name in rcsb_molnames]

    # specify output
    out_dir = f"../../2_MD_Amber/{mutant}/0_prepare/"
    out_name = "protein_named.pdb"

    # create Frame
    new_frame = Frame()

    for ind_mol, mol_rcsb in enumerate(rcsb_model.molecules.filter(mName.is_in("A", "B"))):

        # create Molecule
        new_frame.add_molecule()
        new_molecule = new_frame.molecules[ind_mol]
        new_molecule.name = mol_rcsb.name

        for ind_res, (res_rcsb, res_trj_ref) in enumerate(zip(mol_rcsb.residues, trj_ref.molecules[ind_mol].residues)):
            # create Residue
            new_molecule.add_residue()
            new_residue = new_molecule.residues[ind_res]
            new_residue.name = res_trj_ref.name
            new_residue.id = res_rcsb.id

            for ind_atom, atom_trj_ref in enumerate(res_trj_ref.atoms):
                # create Atom
                new_residue.add_atom()
                new_atom = new_residue.atoms[ind_atom]
                new_atom.name = atom_trj_ref.name
                new_atom.id = atom_trj_ref.id
                new_atom.r = atom_trj_ref.r

    print(new_frame)
    os.makedirs(out_dir, exist_ok=True)
    new_frame.to_pdb(os.path.join(out_dir, out_name))

    del new_frame
