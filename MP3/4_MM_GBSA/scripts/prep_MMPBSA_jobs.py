import remote_runner

# the number of mpi to use in analysis run
# note! this requires amber build with MPI
MPI_NUMBER = 20 

class PrepareMMPBSA(remote_runner.Task):
    def __init__(self, path_to_complex_trj_dir,
                 traj_length=1499,
                 salt=150,
                 mpi_number=MPI_NUMBER
                 ):
        from pyxmolpp2 import PdbFile
        from pathlib import Path
        import os

        self.complex_location = path_to_complex_trj_dir

        self.mmpbsa_dir = os.path.join(path_to_complex_trj_dir, "7_mmpbsa")
        Path(self.mmpbsa_dir).mkdir(exist_ok=True)
        super().__init__(wd=Path(self.mmpbsa_dir))

        self.salt = salt
        self.mpi_number = mpi_number

        pdbfile = PdbFile(os.path.join(f"{self.complex_location}", "0_prepare", "protein.pdb")).frames()[0]
        molecules = [m.residues.size for m in pdbfile.molecules]

        print(f"There are {len(molecules)} molecules")

        molecules = molecules[:-1]
        self.ligand_n_reses = sum(molecules)
        self.complex_n_reses = pdbfile.residues.size

        print(f"Ligand is the first {self.ligand_n_reses} residues and {len(molecules)} molecule(s), "
              f"the whole complex is {self.complex_n_reses} residues")

        self.prep_topology()
        self.prep_trajectory(traj_length)

    def run(self):
        from subprocess import call
        import shlex
        call(shlex.split(self.command), cwd=str(self.wd))
        self.save(self.mmpbsa_dir / 'state.dill')

    def prep_topology(self):
        from subprocess import call
        import shlex

        print("Preparing MMPBSA topology files...")

        antemmpbsa = f"ante-MMPBSA.py -p ../0_prepare/box.prmtop -c com_igb8.top -r rec_igb8.top -l lig_igb8.top " \
                     f"-s ':Na+,Cl-,WAT' -m ':1-{self.ligand_n_reses}' --radii=mbondi3"
        call(shlex.split(antemmpbsa), cwd=str(self.wd))

        antemmpbsa = f"ante-MMPBSA.py -p ../0_prepare/box.prmtop -c com_igb2.top -r rec_igb2.top -l lig_igb2.top " \
                     f"-s ':Na+,Cl-,WAT' -m ':1-{self.ligand_n_reses}' --radii=mbondi2"
        call(shlex.split(antemmpbsa), cwd=str(self.wd))

    def prep_trajectory(self, traj_length):
        from subprocess import call

        ncrst_files = ''

        for i in range(0, traj_length + 1):
            ncrst_files += f'trajin {self.complex_location}/6_run/run{i:05d}.ncrst\n'
        ncrst_files = ncrst_files.strip()
        cpptrajin = f"""parm  {self.complex_location}/0_prepare/box.prmtop
                    {ncrst_files}
                    autoimage
                    center :1-{self.complex_n_reses} mass origin
                    image origin center familiar
                    trajout {str(self.wd)}/mmpbsa.nc
                    run
                    """

        with open(f"{str(self.wd)}/mmpbsa_trajectory.in", "w") as f:
            f.write(cpptrajin)

        print("Preparing trajectory file...")
        call(["cpptraj", "-i", f"{str(self.wd)}/mmpbsa_trajectory.in"])

    def change(self, igb_option):
        from pathlib import Path

        self.input_parameters = f"script{self.salt}_igb{igb_option}.in"

        self.wd = Path(f"{str(self.mmpbsa_dir)}/igb{igb_option}_salt{self.salt}")

        self.wd.mkdir(exist_ok=True)
        self.input_parameters = "../../../" + self.input_parameters

        self.command = f"mpirun -np {self.mpi_number} " \
                       f"MMPBSA.py.MPI -O -i {self.input_parameters} " \
                       f"-o complex_MMPBSA.dat -sp ../../0_prepare/box.prmtop " \
                       f"-cp ../com_igb{igb_option}.top " \
                       f"-rp ../rec_igb{igb_option}.top " \
                       f"-lp ../lig_igb{igb_option}.top " \
                       f"-y ../mmpbsa.nc"

        print(self.command)


if __name__ == '__main__':
    import os

    prefix = "../../2_MD_Amber" # specifies the location of MD dirs

    # scpecify mutants for analysis
    # mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]
    mutants = ["wt+mp3", "delta_plus+mp3"]

    for mutant in mutants:

        # NOTE: double check generated .in files
        # change frames for trajectory prep
        task = PrepareMMPBSA(path_to_complex_trj_dir=os.path.join(prefix, mutant),
                             traj_length=1499) # initiate task

        task.change(igb_option=2) # change igb to 2 and save job
        task.save(task.wd / 'state.dill')

        task.change(igb_option=8) # change igb to 8 and save job
        task.save(task.wd / 'state.dill')

        # copy .in script into desired directories
        os.system(f'cp script150_igb2.in {os.path.join(os.path.join(prefix, mutant), "7_mmpbsa")}')  
        os.system(f'cp script150_igb8.in {os.path.join(os.path.join(prefix, mutant), "7_mmpbsa")}')   
