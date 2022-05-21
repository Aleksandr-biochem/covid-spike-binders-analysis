import remote_runner

# the number of mpi to use in analysis run
# note! this requires amber build with MPI
MPI_NUMBER = 20

# this variable specifyes interface residues indexes of wt_mp3 comlex
# these are used for all complexes decomposition
roi = '1-24,26-38,40-43,45,48,52,56,133-137,139-140,145-153,174-180,182-191,203-210,212,214-238'

class PrepareMMPBSA(remote_runner.Task):
    def __init__(self, complex_name, files_location, traj_length=1499, short=False, salt=150, mpi_number=MPI_NUMBER):

        from pyxmolpp2 import PdbFile
        from pathlib import Path

        self.complex_location = files_location + complex_name
        self.salt = salt
        self.mpi_number = mpi_number
        self.distance_cutoff = 6
        mmpbsa_dir = f"{self.complex_location}/8_mmpbsa_decomp"
        self.short = short
        if self.short:
            mmpbsa_dir += '_short'
        Path(mmpbsa_dir).mkdir(exist_ok=True)
        remote_runner.Task.__init__(self, wd=Path(mmpbsa_dir))
        self.mmpbsa_dir = self.wd

        pdbfile = PdbFile(f"{self.complex_location}/0_prepare/protein.pdb").frames()[0]
        molecules = [m.residues.size for m in pdbfile.molecules]
        print(f"There are {len(molecules)} molecules")
        molecules = molecules[:-1]
        self.ligand_n_reses = sum(molecules)
        self.complex_n_reses = pdbfile.residues.size
        print(f"Ligand is the first {self.ligand_n_reses} residues and {len(molecules)} molecule(s), "
              f"the whole complex is {self.complex_n_reses} residues")

        self.prep_topology()
        self.prep_trajectory(traj_length)
        self.residues_of_interest = roi
        # self.prep_residues_of_interest(traj_length)
        

    def run(self):
        from subprocess import call
        import shlex
        call(shlex.split(self.command), cwd=str(self.wd))
        self.save(self.mmpbsa_dir / 'state.dill')

    def prep_topology(self):
        from subprocess import call
        import pathlib, shlex

        # mbondi3 for igb=8, and mbondi2 for igb=2
        print("Preparing MMPBSA topology files...")

        antemmpbsa = f"ante-MMPBSA.py -p ../0_prepare/box.prmtop -c com_igb8.top -r rec_igb8.top -l lig_igb8.top " \
                     f"-s ':Na+,Cl-,WAT' -m ':1-{self.ligand_n_reses}' --radii=mbondi3"
        call(shlex.split(antemmpbsa), cwd=str(self.wd))

        antemmpbsa = f"ante-MMPBSA.py -p ../0_prepare/box.prmtop -c com_igb2.top -r rec_igb2.top -l lig_igb2.top " \
                     f"-s ':Na+,Cl-,WAT' -m ':1-{self.ligand_n_reses}' --radii=mbondi2"
        call(shlex.split(antemmpbsa), cwd=str(self.wd))

    def prep_trajectory(self, traj_length):
        from subprocess import call
        from pathlib import Path
        ncrst_files = ''
        if self.short:
            ncs = Path(f'{self.complex_location}/6_run_short').glob('*.ncrst')
            for i in range(0, traj_length + 1):
                ncrst_files += f'trajin {self.complex_location}/6_run_short/run_short{i:05d}.ncrst\n'
        else:
            ncs = Path(f'{self.complex_location}/6_run').glob('*.ncrst')
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
        # trajout {str(self.wd)}/mmpbsa_{skip_frames}.nc start {skip_frames + 1}
        with open(f"{str(self.wd)}/mmpbsa_trajectory.in", "w") as f:
            f.write(cpptrajin)
        print("Preparing trajectory file...")
        call(["cpptraj", "-i", f"{str(self.wd)}/mmpbsa_trajectory.in"])

    def prep_residues_of_interest(self, traj_length):
        from pyxmolpp2 import PdbFile, AtomPredicate, mName, Trajectory, AmberNetCDF
        from tqdm import tqdm
        import os
        print("Preparing residues of interest...")

        chain_counter = 0
        path_to_trj_ref = f"{self.complex_location}/0_prepare/box.pdb"
        path_to_xray_ref = f"{self.complex_location}/0_prepare/protein_named.pdb"
        xray_ref = PdbFile(path_to_xray_ref).frames()[0]
        trj_ref = PdbFile(path_to_trj_ref).frames()[0]
        chain_list = []
        for m_dest, m_src in zip(trj_ref.molecules, xray_ref.molecules):
            chain_list.append(m_src.name)
            m_dest.name = m_src.name
        chain_set = set(chain_list)
        heavy_atoms = {}
        for chain in chain_list:
            heavy_atoms[chain] = AtomPredicate(lambda a: not a.name.startswith("H")) & mName.is_in(chain)
        traj = Trajectory(trj_ref)
        traj.extend(AmberNetCDF(os.path.join(str(self.wd), f"mmpbsa.nc")))

        def extract_roi(traj_):
            residues_of_interest = []
            index = traj_[0].index
            print(index)
            for frame in tqdm(traj_, desc=f"{self.complex_location} residues of interest extraction"):
                if frame.index == index:
                    # create variables for selections from atoms of frame
                    frame_ats = frame.molecules.atoms
                    frame_by_chain = {}
                    for chain in chain_list:
                        frame_by_chain[chain] = frame_ats.filter(heavy_atoms[chain])
                for chain in chain_list[:-1]:
                    for atom_lig in frame_by_chain[chain]:
                        for atom_rbds in frame_by_chain[chain_list[-1]]:
                            if atom_lig.r.distance(atom_rbds.r) < self.distance_cutoff:
                                lId = str(atom_lig.residue.id)
                                rId = str(atom_rbds.residue.id)
                                if lId not in residues_of_interest:
                                    residues_of_interest.append(lId)
                                if rId not in residues_of_interest:
                                    residues_of_interest.append(rId)

            return residues_of_interest

        residues_of_interest = extract_roi(traj)
        residues_of_interest = sorted([int(r) for r in residues_of_interest])
        print(residues_of_interest)
        residues_of_interest.append(residues_of_interest[-1])
        ranges = []
        left = residues_of_interest[0]
        right = left
        for r1 in residues_of_interest[1:]:
            if right + 1 != r1 and right + 2 != r1:
                if left != right:
                    ranges.append(f"{left}-{right}")
                else:
                    ranges.append(str(left))
                left = r1
                right = left
            else:
                right = r1

        self.residues_of_interest = ', '.join(ranges)
        print(self.residues_of_interest)

    def change(self, igb_option, decomp):
        from pathlib import Path

        self.input_parameters = f"script{self.salt}_igb{igb_option}.in"
        self.wd = Path(f"{str(self.mmpbsa_dir)}/igb{igb_option}_salt{self.salt}")
        if decomp:
            self.wd = Path(str(self.wd) + '_decomp')
            self.wd.mkdir(exist_ok=True)
            with open(self.input_parameters) as f:
                input_text = f.readlines()
                input_text = ''.join(input_text)
                input_text = input_text.strip()
                decomp_lines = f"""
&decomp
idecomp=4, print_res="{self.residues_of_interest}"
dec_verbose=3,
/
"""
                input_text += f'\n{decomp_lines}'
                with open(self.wd / self.input_parameters, 'w') as f2:
                    f2.write(input_text)
        else:
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

    prefix = "../../2_MD_Amber/" # specifies the location of MD dirs

    # scpecify mutants for analysis
    mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

    for mutant in mutants:
        # NOTE: double check generated decomp files
        # change frames for trajectory prep
        task = PrepareMMPBSA(mutant, prefix)

        task.change(igb_option=2, decomp=True)
        task.save(task.wd / 'state.dill')

        task.change(igb_option=8, decomp=True)
        task.save(task.wd / 'state.dill')
