from pyxmolpp2.pipe import TrajectoryProcessor
from pyxmolpp2 import Frame


class CalcRmsd(TrajectoryProcessor):

    def __init__(self, reference: Frame):
        self.reference = reference

    def before_first_iteration(self, frame: Frame):
        # open file for rmsd writing
        ...

    def after_last_iteration(self, frame: Frame):
        # close file
        ...

    def __call__(self, frame: Frame):
        ...
        return frame


if __name__ == "__main__":
    from pyxmolpp2 import Trajectory, PdbFile, TrjtoolDatFile
    from pyxmolpp2.pipe import Run
    import os

    #  setup trajectory parameters
    path_to_trajectory = "path to  your reference"
    path_to_reference = "path to your trajectory"
    trajectory_start = 1
    trajectory_length = 10

    #  load trajectory
    reference = PdbFile(os.path.join(path_to_trajectory, "1_build/ref.pdb"))
    trajectory = Trajectory(reference.frames()[0])
    for ind in range(trajectory_start, trajectory_length + 1):
        fname = "{pattern}.{filetype}".format(pattern="run%05d", filetype="dat")
        trajectory.extend(TrjtoolDatFile(os.path.join(os.path.join(path_to_trajectory, "5_run"), fname % (ind))))

    trajectory | CalcRmsd(reference=reference) | Run()
