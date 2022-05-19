from pyxmolpp2 import PdbFile, rId, Frame, ResidueId, mName, rName
import os
import argparse

def renumber_resID():
    """
    Renumber prepared MD .pdb structure according to default numbering.
    Generates new file protein_named.pdb

    For wt and omicron also new directory is created.
    """
    # specify mutant
    mutants = ["wt","alpha", "delta", "delta_plus", "omicron"]

    # read rcsb structure
    path_to_rcsb_ref = "../../0_prepare_structures/structures/7jzu.pdb"
    rcsb_model = PdbFile(path_to_rcsb_ref).frames()[0]
    rcsb_molnames = [mol.name for mol in rcsb_model.molecules.filter(mName.is_in("A", "B"))]

    for mutant in mutants:
        name = mutant.split('+')[0]
        # read ref trj pdb
        path_to_trj_ref = f"../../2_MD_Amber/sample/{name}+mp1/0_prepare/protein.pdb"
        trj_ref = PdbFile(path_to_trj_ref).frames()[0]
        trj_ref_mols = [rcsb_model.molecules.filter(mName == mol_name)[0] for mol_name in rcsb_molnames]

        # specify output
        if mutant in ["omicron+lcb1", "wt+lcb1"]:
            if not os.path.exists(f"../../2_MD_Amber/sample/{name}+mp1_2/0_prepare/"):
                os.mkdir(f"../../2_MD_Amber/sample/{name}+mp1_2/0_prepare/")
           
            out_dir = f"../../2_MD_Amber/sample/{name}+mp1_2/0_prepare/"
        else:
            out_dir = f"../../2_MD_Amber/sample/{name}+mp1/0_prepare/"
        
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
                # new_residue.id = res.id
                new_residue.id = res_rcsb.id

                for ind_atom, atom_trj_ref in enumerate(res_trj_ref.atoms):
                    # create Atom
                    new_residue.add_atom()
                    new_atom = new_residue.atoms[ind_atom]
                    new_atom.name = atom_trj_ref.name
                    new_atom.id = atom_trj_ref.id
                    # new_atom.id = ind_atom + 1
                    new_atom.r = atom_trj_ref.r

        print(new_frame)
        os.makedirs(out_dir, exist_ok=True)
        new_frame.to_pdb(os.path.join(out_dir, out_name))

        del new_frame


if __name__ == "__main__":
    renumber_resID()