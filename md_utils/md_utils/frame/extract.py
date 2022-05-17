from pyxmolpp2 import Frame, mName, Molecule, MoleculeSelection, ResidueSelection, AtomSelection, Atom, aId
from scipy.spatial import cKDTree
from typing import List, Union, Tuple


def extract_one_letter_amino_acid_seq(frame: Frame,
                                      molname: mName) -> list:
    """
    extract one letter amino asids sequence for the specific molecule in frame
    :param frame:
    :param molname: name of molecule in frame
    :return: list of one letter amino acids for the specific molecule in frame
    """
    # dictionary to convert canonical residue names
    three_to_one_leter_name = {"GLY": "G",
                               "ALA": "A",
                               "VAL": "V",
                               "LEU": "L",
                               "ILE": "I",
                               "SER": "S",
                               "THR": "T",
                               "MET": "M",
                               "CYS": "C",
                               "ASN": "N",
                               "GLN": "Q",
                               "ASP": "D",
                               "GLU": "E",
                               "LYS": "K",
                               "ARG": "R",
                               "HIS": "H",
                               "PHE": "F",
                               "TRP": "W",
                               "TYR": "Y",
                               "PRO": "P",
                               }

    return [three_to_one_leter_name[residue.name] for residue in frame.molecules.filter(mName == molname).residues]


def extract_residues_on_interface(partner_A: Union[Molecule, MoleculeSelection],
                                  partner_B: Union[Molecule, MoleculeSelection],
                                  cutoff: float
                                  ) -> List[ResidueSelection]:
    """
    Extract residues on interface of complex between partner_A and partner_B
    :param partner_A: first partner of intermolecular interaction
    :param partner_B: second partner of intermolecular interaction
    :param cutoff: distance cutoff in angstroms
    :return: list of [residues of partner_A, residues of partner_B] - residues on interface
                  of interaction between partner_A and partner_B within the distance cutoff
    """
    # get coordinates and construct binary trees
    a_tree = cKDTree(partner_A.atoms.coords.values)
    b_tree = cKDTree(partner_B.atoms.coords.values)

    # compute distances with a cut-off
    dist = a_tree.sparse_distance_matrix(b_tree, max_distance=cutoff)
    dist = dist.toarray()

    # get atom indexes in cut-off range
    a_index, b_index = dist.nonzero()

    # get corresponding residues
    first_a_id = partner_A.atoms[0].id  # to correct indexes from matrix
    resdue_selection_A = partner_A.atoms.filter(aId.is_in(set(a_index + first_a_id))).residues

    first_b_id = partner_B.atoms[0].id  # to correct indexes from matrix
    resdue_selection_B = partner_B.atoms.filter(aId.is_in(set(b_index + first_b_id))).residues

    return [resdue_selection_A, resdue_selection_B]


def extract_hydrogen_bonds(partner_A: Union[Molecule, MoleculeSelection],
                           partner_B: Union[Molecule, MoleculeSelection],
                           distance_cutoff=3,
                           angle_cutoff=30,
                           ) -> List[Tuple[Union[Atom, AtomSelection]]]:
    """
    extract hydrogen bond is formed between an atom with a hydrogen bonded to it (the donor, D)
            and another atom (the acceptor, A)
    :param partner_A:
    :param partner_B:
    :param distance_cutoff: distance D-A is less than the cut-off distance (default 3.0 Angstroms)
    :param angle_cutoff: angle D-H-A is less than the cut-off angle (default 30 degrees)
    :return:
    """
    ...


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile, mName

    path_to_pdb = "path to your pdb"
    frame = PdbFile(path_to_pdb).frames()[0]

    # example usage of extract_residues_on_interface
    interface_residues_partner_A, interface_residues_partner_B = extract_residues_on_interface(
        partner_A=frame.molecules.filter(mName == "A"),
        partner_B=frame.molecules.filter(mName == "B"),
        cutoff=5)

    # example usage of extract_one_letter_amino_acid_seq
    amino_acid_seq = extract_one_letter_amino_acid_seq(frame, molname="A")
