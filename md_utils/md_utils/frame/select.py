from typing import List
from pyxmolpp2 import AtomSelection, AtomPredicate, mName, aName, rId, PdbFile, Frame
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import functools
import os
import operator


def get_sec_str_residue_ids(frame: Frame,
                            molname: str
                            ) -> List[int]:
    """
    to predict the secondary structure residues use DSSP https://biopython.org/docs/1.76/api/Bio.PDB.DSSP.html
    :param frame:
    :param molname:
    :return: list of residue id, which belongs to secondary-structure region
             helix - H (Alpha helix (4-12), G (3-10 helix) and I (Pi helix)
             or
             strand E (Strand) and B (Isolated beta-bridge residue),
    """
    # predict secondary structure from input pdb
    frame.to_pdb("temp.pdb")
    dssp_tuple = dssp_dict_from_pdb_file("temp.pdb")
    dssp_dict = dssp_tuple[0]

    # get structured residues from target molecule
    sec_str_rids = []

    for key in dssp_dict.keys():
        if (key[0] == molname) and (dssp_dict[key][1] in ['H', 'G', 'I', 'E', 'B']):
            sec_str_rids.append(key[1][1])

    os.remove("temp.pdb")

    return sec_str_rids


def get_sec_str_residues_predicate(frame: Frame,
                                   molnames: list
                                   ) -> AtomPredicate:
    """
    :param frame:
    :param molnames: list of molecules names. for example ["A"] or ["A", "B"]
                    For example, reference structure is complex so that it contains two molecules (A and B)
    :return: CA atoms belongs to secondary structured regions
    """
    sec_str_residues = []
    for molname in molnames:
        sec_str_rids = get_sec_str_residue_ids(frame=frame,
                                               molname=molname)
        sec_str_residues.append((rId.is_in(set(sec_str_rids))) & (mName == molname))

    sec_str_residues_predicate = functools.reduce(operator.or_, sec_str_residues)

    return sec_str_residues_predicate


def select_sec_str_ca_atoms(frame: Frame,
                            molnames: list,
                            ) -> AtomSelection:
    """
    :param frame:
    :param molnames: list of molecules names. for example ["A"] or ["A", "B"]
                    For example, reference structure is complex so that it contains two molecules (A and B)
    :return: CA atoms belongs to secondary structured regions
    """

    sec_str_residues_predicate = get_sec_str_residues_predicate(frame, molnames)
    return frame.atoms.filter(sec_str_residues_predicate & (aName == "CA"))


if __name__ == "__main__":
    path_to_pdb = "7jzu.pdb"
    frame = PdbFile(path_to_pdb).frames()[0]

    sec_str_rids = get_sec_str_residue_ids(frame, molname="A")
    sec_str_CA_atoms = select_sec_str_ca_atoms(frame, molnames=["A"])
