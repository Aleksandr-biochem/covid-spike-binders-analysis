from amber_runner.MD import *

class Prepare(Step):
    def run(self, md):
        from subprocess import call
        from pyxmolpp2 import PdbFile, rName
        import gemmi, propka.run
        from typing import Tuple, Iterator

        def extract_pkas(protein, conformation, parameters, target_ph):
            rename_residues = ['GLU', 'ASP', 'HIS', 'CYS']
            rename_map_acids = {
                'GLU': 'GLH',
                'ASP': 'ASH',
            }
            protonation = {}
            str_ = ""
            # printing pKa summary
            for residue_type in parameters.write_out_order:
                for group in protein.conformations[conformation].groups:
                    if group.residue_type == residue_type:
                        # print(group.coupled_titrating_group)
                        str__ = f"{group.residue_type:>9s} {group.atom.res_num:>4d} " \
                                f"{group.atom.chain_id:>9s} {group.pka_value:8.2f}\n"
                        str_ += str__
                        if group.residue_type in rename_residues:
                            if group.pka_value < target_ph:  # deprotonated
                                # protonation[(group.residue_type, group.atom.chain_id, group.atom.res_num)] = False
                                if group.residue_type == 'CYS':
                                    protonation[(group.residue_type, group.atom.chain_id, group.atom.res_num)] = 'CYM'
                            else:
                                # protonation[(group.residue_type, group.atom.chain_id, group.atom.res_num)] = True
                                if group.residue_type in rename_map_acids:  # protonated
                                    protonation[(group.residue_type, group.atom.chain_id, group.atom.res_num)] = \
                                        rename_map_acids[group.residue_type]
                                if group.residue_type == 'HIS':
                                    protonation[(group.residue_type, group.atom.chain_id, group.atom.res_num)] = 'HIP'
                        else:
                            print(str__.strip())
            return protonation

        def find_ss_bond_pairs(st: gemmi.Structure) -> Iterator[Tuple[int, int]]:
            cyx_residues = []
            for chain in st[0]:
                for residue in chain:
                    if residue.name == "CYX":
                        cyx_residues.append(residue)

            for i in range(len(cyx_residues)):
                ri = cyx_residues[i]
                sgi = ri["SG"][0]
                for j in range(i + 1, len(cyx_residues)):
                    rj = cyx_residues[j]
                    sgj = rj["SG"][0]
                    if sgi.pos.dist(sgj.pos) < 2.3:
                        yield ri.seqid.num, rj.seqid.num

        def get_residues(st: gemmi.Structure):
            assert len(st) == 1, "Structure MUST have one MODEL"
            for chain in st[0]:
                for residue in chain:
                    yield residue

        def gemmi_write_pdb(st: gemmi.Structure, *args, **kwargs):
            """ A workaround for https://github.com/project-gemmi/gemmi/issues/39 """
            types = []
            for residue in get_residues(st):
                types.append(residue.entity_type)
                residue.entity_type = gemmi.EntityType.Polymer

            st.write_pdb(*args, **kwargs)

            for residue, t in zip(get_residues(st), types):
                residue.entity_type = t

        def rename_residues_amber_to_standard(st):
            rename_map = {
                'GLH': 'GLU',
                'ASH': 'ASP',
                'HID': 'HIS',
                'HIP': 'HIS',
                'HIE': 'HIS',
                'CYX': 'CYS',
                'CYM': 'CYS',
                'G5': 'G',  # RNA 5-end
                'A5': 'A',
                'C5': 'C',
                'U5': 'U',
                'G3': 'G',  # RNA 3-end
                'A3': 'A',
                'C3': 'C',
                'U3': 'U',
                'DG5': 'DG',  # DNA 5-end
                'DA5': 'DA',
                'DC5': 'DC',
                'DT5': 'DT',
                'DG3': 'DG',  # DNA 3-end
                'DA3': 'DA',
                'DC3': 'DC',
                'DT3': 'DT',
            }

            for res in get_residues(st):
                res.name = rename_map.get(res.name, res.name)

        def build_tleaprc(pdb_fn, rc_fn, prefix='box', ss_bond_commands='', n_na=0, n_cl=0):
            print("Building parameters...")
            tleaprc = f"""source leaprc.protein.ff14SB
source leaprc.water.tip4pew
wbox = loadpdb {pdb_fn}
{ss_bond_commands}
solvateoct wbox TIP4PEWBOX 15
addions wbox Na+ {n_na}
addions wbox Cl- {n_cl}
saveamberparm wbox {prefix}.prmtop {prefix}.inpcrd
savepdb wbox {prefix}.pdb
quit"""

            with open(rc_fn, "w") as f:
                f.write(tleaprc)

        # 1: determine disulfides

        md.init_struct_pdb4amber = "pdb4amber.pdb"
        md.init_struct_propka = "propka.pdb"
        md.init_struct_proper_states = "proper_states.pdb"

        call(["pdb4amber", "-d", "--most-populous",
              "-i", md.init_struct, '-o', md.init_struct_pdb4amber, ],
             cwd=str(self.step_dir))
        st: gemmi.Structure = gemmi.read_pdb(str(self.step_dir / md.init_struct_pdb4amber), split_chain_on_ter=True)
        ss_bonds = find_ss_bond_pairs(st)
        ss_bond_commands = "\n".join(f"bond wbox.{i}.SG wbox.{j}.SG" for i, j in ss_bonds)
        ss_bond_commands = ""
        rename_residues_amber_to_standard(st)
        gemmi_write_pdb(st, str(self.step_dir / md.init_struct_pdb4amber), numbered_ter=False)

        # 2: determine protonation states with Propka

        st: gemmi.Structure = gemmi.read_pdb(str(self.step_dir / md.init_struct_pdb4amber), split_chain_on_ter=True)
        mol = propka.run.single(str(self.step_dir / md.init_struct_pdb4amber), write_pka=False)
        protein, conformation, parameters = mol, 'AVR', mol.version.parameters
        protonation_states = extract_pkas(protein, conformation, parameters, 7.5)
        print(f"Protonation states: {protonation_states}")
        chains = set([k[1] for k in protonation_states.keys()])
        for c in st[0]:
            if c.name in chains:
                for r in c:
                    key = (str(r.name), str(c.name), r.seqid.num)
                    if key in protonation_states.keys():
                        r.name = protonation_states[key]
        gemmi_write_pdb(st, str(self.step_dir / md.init_struct_propka), numbered_ter=False)
        call(["pdb4amber", "-d", "-i", md.init_struct_propka, '-o', md.init_struct_proper_states, ],
             cwd=str(self.step_dir))

        # ss_bond_commands = ''

        # 3: solvate to estimate charge and number of water molecules

        build_tleaprc(md.init_struct_proper_states, str(self.step_dir / "draft_tleap.rc"), 'draft_box', ss_bond_commands)

        call(["tleap", "-s", "-f", "draft_tleap.rc"], cwd=str(self.step_dir))

        pdb = PdbFile(str(self.step_dir / "draft_box.pdb")).frames()[0]
        residues = pdb.residues
        n_h2o = residues.filter(rName == "WAT").size
        n_na = residues.filter(rName == "Na+").size
        n_cl = residues.filter(rName == "Cl-").size

        V0 = 3.00050910373387e-26
        Na = 6.022140857e23
        salt_concentration = 0.15
        if n_na > 0:
            n_cl = int(salt_concentration * Na * n_h2o * V0)
            n_na += n_cl
        else:
            n_na = int(salt_concentration * Na * n_h2o * V0)
            n_cl += n_na

        # 4: re-solvate the protein with proper salt concentration

        build_tleaprc(md.init_struct_proper_states, str(self.step_dir / "tleap.rc"), 'box', ss_bond_commands, n_na, n_cl)
        call(["tleap", "-s", "-f", "tleap.rc"], cwd=str(self.step_dir))

        # 5: write the resultant box pdb-file
        md.sander.prmtop = self.step_dir / "box.prmtop"
        md.sander.inpcrd = self.step_dir / "box.inpcrd"
        # self.parm7_rst7_to_pdb(str(md.sander.prmtop), str(md.sander.inpcrd),
        #                        str(self.step_dir / "box.pdb"))

        # 6a: calculate the number of protein residues and dump only that part to protein.*  files with cpptraj

        pdb = PdbFile(str(self.step_dir / md.init_struct_proper_states)).frames()[0]
        nreses = pdb.residues.size
        print(f"Protein residues: {nreses}")

        cpptrajin = f"""parm {str(md.sander.prmtop)}
trajin {str(md.sander.inpcrd)}
strip !:1-{nreses} parmout {str(self.step_dir / "protein.prmtop")}
trajout {str(self.step_dir / "protein.inpcrd")} restart
trajout {str(self.step_dir / "protein.pdb")} pdb
run
"""

        with open(str(self.step_dir / "strip_solvent.in"), "w") as f:
            f.write(cpptrajin)
        call(["cpptraj", "-i", str(self.step_dir / "strip_solvent.in")])

        # 6b: calculate the number of protein atoms; to be used if only protein will be dumped in production

        num_atoms = PdbFile(str(self.step_dir / "protein.pdb")).n_atoms()
        print(f"Protein atoms: {num_atoms}")

        print("Building done.")

        # 7: simulation parameters

        md.min_1.input.cntrl(
            imin=1,
            maxcyc=1000,
            ncyc=1000,
            ntb=1,
            ntr=1,
            cut=10.0,
            restraint_wt=200,
            restraintmask=":1-{:d}".format(nreses),
        ) # restrain protein to move solvent and relax clashes
        md.min_2.input.cntrl(
            imin=1,
            maxcyc=1000,
            ncyc=1000,
            ntb=1,
            ntr=0,
            cut=10.0
        ) # minimize all
        md.heat.input.cntrl(
            imin=0,
            irest=0,
            ntx=1,
            ntb=1,
            ntr=1,
            ntc=2,
            ntf=2,
            tempi=0.0,
            temp0=298.0,
            ntt=1,
            nstlim=10000,
            dt=0.001,
            ntpr=50,
            ntwx=50,
            ntwr=50,
            ioutfm=1,
            restraint_wt=10,
            restraintmask=":1-{:d}".format(nreses),
        ) 
        md.equil.input.cntrl(
            imin=0,
            irest=1,  # continue simulation
            cut=10.5,  # cutoff radius
            ntx=5,  # read coordinates and velocities from restart file
            ntb=2,  # constant pressure
            iwrap=1,  # wrap periodic images into unit box
            ntt=11,  # w/ termostat Bussi
            tempi=298.0,
            temp0=298.0,
            ntp=1,  # isotropic barostat
            ntr=0,  # w/o restraints
            ntc=2,  # shake only for hydrogens
            ntf=2,  # bond force is shake
            nstlim=500000,  # number of steps
            dt=0.002,  # dt
            ntpr=500,  # md out print frequency
            ntwx=500,  # md crd output frequency
            ntwr=500000,  # md restart frequency
            ioutfm=1
        )
        md.production.input.cntrl(
            imin=0,
            irest=1,  # continue simulation
            cut=10.5,  # cutoff radius
            ntx=5,  # read coordinates and velocities from restart file
            ntb=2,  # constant pressure
            iwrap=1,  # wrap periodic images into unit box
            ntt=11,  # w/ termostat Bussi
            tempi=298.0,
            temp0=298.0,
            ntp=1,  # isotropic barostat
            ntr=0,  # w/o restraints
            ntc=2,  # shake only for hydrogens
            ntf=2,  # bond force is shake
            nstlim=500000,  # number of steps
            dt=0.002,  # dt
            ntpr=500,  # md out print frequency
            ntwx=500,  # md crd output frequency
            ntwr=500000,  # md restart frequency
            ioutfm=1,
            ntwprt=num_atoms # the amount of atoms to save, save just protein
        )

    def parm7_rst7_to_pdb(self, prmtop: str, inpcrd: str, pdb_name: str):
        from subprocess import call
        cpptrajin = f"""parm {prmtop}
trajin {inpcrd}
trajout {pdb_name} pdb include_ep
run
"""
        with open(str(self.step_dir / "cpptraj.in"), "w") as f:
            f.write(cpptrajin)
        call(["cpptraj", "-i", self.step_dir / "cpptraj.in"])


class Summary(Step):

    def run(self, md):
        import subprocess

        summary = md.mkdir(self.step_dir).absolute()

        subprocess.run([
            "process_mdout.perl",
            (md.heat.step_dir / f"{md.heat.name}.out").absolute(),
            (md.equil.step_dir / f"{md.equil.name}.out").absolute()
        ], cwd=str(summary))

        with remote_runner.utility.ChangeDirectory(self.step_dir):
            vol = self.read_summary(summary / "summary.VOLUME")
            temp = self.read_summary(summary / "summary.TEMP")
            density = self.read_summary(summary / "summary.DENSITY")
            etot = self.read_summary(summary / "summary.ETOT")
            eptot = self.read_summary(summary / "summary.EPTOT")
            ektot = self.read_summary(summary / "summary.EKTOT")

            import matplotlib.pyplot as plt

            plt.plot(vol.time, vol.value)
            plt.ylabel("volume, A^3")
            plt.xlabel("time, ps")
            plt.savefig("VOLUME.png")
            plt.close()

            plt.plot(temp.time, temp.value)
            plt.ylabel("temperature, K")
            plt.xlabel("time, ps")
            plt.savefig("TEMP.png")
            plt.close()

            plt.plot(density.time, density.value)
            plt.ylabel("density, g/cm^3")
            plt.xlabel("time, ps")
            plt.savefig("DENSITY.png")
            plt.close()

            plt.plot(ektot.time, ektot.value, label="EKTOT", color="red")
            plt.plot(eptot.time, eptot.value, label="EPTOT", color="green")
            plt.plot(etot.time, etot.value, label="ETOT", color="black")
            plt.ylabel("energy, Kcal/mol")
            plt.xlabel("time, ps")
            plt.savefig("ENERGY.png")
            plt.close()

    @staticmethod
    def read_summary(summary_path):
        import pandas as pd
        return pd.read_csv(summary_path, sep=r"\s+", names=["time", "value"])


class variant(MdProtocol):
    def __init__(self, path_to_initial_dir, mutant_dir, name):
        from amber_runner.command import LambdaStringArgument
        import os

        wd = Path(mutant_dir)
        self.mkdir(wd)
        print("Create directory: " + str(wd))
        MdProtocol.__init__(self, name=name, wd=wd)

        self.sander = PmemdCommand()
        self.sander.executable = ["/opt/amber20/bin/pmemd.cuda"]
        self.sander.allow_small_box = True

        self.init_struct = os.path.join(path_to_initial_dir,  f"{name}_minimized.pdb")

        self.sander.refc = LambdaStringArgument("-ref", lambda: self.sander.inpcrd)
        self.sander.mdinfo = LambdaStringArgument("-inf", lambda: f"{self.sander.output_prefix}.mdinfo")

        self.prepare = Prepare("prepare")
        self.min_1 = SingleSanderCall("min_1")
        self.min_2 = SingleSanderCall("min_2")
        self.heat = SingleSanderCall("heat")
        self.equil = SingleSanderCall("equil")
        self.summary = Summary("summary")
        self.production = RepeatedSanderCall("run", number_of_steps=1500)


if __name__ == '__main__':

    # specify path to folder with initial structures prepared at the step 0
    path_to_initial_dir = "../../0_prepare_structures/"

    # mutants for simulation
    mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

    for name in mutants:
        mutant_dir = f"../{name}"
        md = variant(path_to_initial_dir, mutant_dir, name)
        md.save(Path(mutant_dir + "/state.dill"))
