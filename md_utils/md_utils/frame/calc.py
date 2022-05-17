from pyxmolpp2 import Frame, AtomPredicate, calc_rmsd
from select import get_sec_str_residues_predicate


def get_rmsd(reference_frame: Frame,
             probe_frame: Frame,
             by_atoms: AtomPredicate,
             ) -> float:
    """
    Comparison similarity in teo three-dimensional structures by the RMSD after optimal rigid body superposition
    https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
    :param reference_frame: for example X-ray structure or any other case
    :param probe_frame: for example frame from MD trajectory or any other case
    :param by_atoms: rigid body superposition is carried out by this atoms. In case of proteins,
                     this often calculated for the secondary-structure Cα atoms.
    :return: rmsd value by_atoms between reference_frame and probe_frame
    """
    ref_atoms = reference_frame.atoms.filter(by_atoms)  # reference atoms for alignment
    frame_atoms = probe_frame.atoms.filter(by_atoms)  # selection from frame for alignment

    alignment = frame_atoms.alignment_to(ref_atoms)  # alignment
    crd = frame_atoms.coords.values.copy()
    crd = crd @ alignment.matrix3d().T + alignment.vector3d().values  # get coordinates

    return calc_rmsd(ref_atoms.coords.values, crd)


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile

    # test
    frame1 = PdbFile("path1").frames()[0]
    frame2 = PdbFile("path2").frames()[0]
    pred = get_sec_str_residues_predicate(frame=frame1, molnames=["B"])
    print(get_rmsd(frame1, frame2, pred))
