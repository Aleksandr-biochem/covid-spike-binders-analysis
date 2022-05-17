# The repository of covid-spike-binders-analysis project
**Authors:**
- Aleksandr Kovalenko *(All-Russian Institite of Plant Protection, Saint Petersburg, Russia)*
- Xenia Sukhanova *(ITMO University, Saint Petersburg, Russia)*  

**Supervisors:** Olga Lebedenko, Nikolai Skrynnikov *(BioNMR laboratory, Saint Petersburg State University, Saint Petersburg, Russia)*

### Background:
This repository contains results and materials for analysis reproduction for the covid-spike-binders-analysis project, which investigates mini-protein ligands binding to several SARS-CoV-2 RBD variants.

Mini-proteins 1 and 3 (MP1/3) are artificial ligands initially designed to bind the wild type SARS-CoV-2 RBD. However, the emergence of new covid strains hinders their effective use. 

**Goal:** This project aims to provide a workflow for the binding analysis of MP1/3 to different variants of SARS-CoV-2 and also for MPs structural adaptation to new strains.

**Objectives:**
- Analyze MP1/3 binding to several RBD variants via MD simulaition, MM-GBSA and decomposition analysis and investigate whether binding in strain-dependent
- Perform mutation of MP1/3 to improve ligand binding to recent SARS-CoV-2 variants delta+ and omicron
- Develop a collection of python sscripts for structural analysis, simulations and results processing and visualization 

### Methods: 
The workflow consists of 6 steps, highlighted at the diagram:

/images/workflow.png

0. Structures of different MP/RBD complex variants are prepared from Cryo-EM solved structures of MP1/3 with wild type RBD available in PDB. This step utilizes `Biobb` package for structure checking and mutation and `Phenix` software for structure physical relaxation
1. Mutated structures can be subjected to mutational scan with `FlexddG` method within `Rosetta` sortware, which allows to identify mutants with improved binding
2. For all structures MD simulations are performed with `Amber20`
3. RMSD metrics are assesed for simulated trajectories
4. Using conformations obtained from MD, MM-GBSA analysis i `Amber20` is performed to estimate binding free energies in different complexes 
5. For a detailed analysis of residue interactions contribution in binding MM-GBSA decomposition analysis is performed

These steps utilize two modeling softwares (Rosetta and Amber20) and several methods within them and provide a panel of evidence on MP binding to RBD.

### System requirements:
**Key packages and programs:**
- the majority of scripts are written on `Python3`, there are also `bash` scripts
- `slurm` cluster management and job scheduling system
- `Rosetta` and `Amber20` (amber assembly with MPI is highly preferrable)
- [pyxmolpp2](https://github.com/sizmailov/pyxmolpp2) python library (note thath this library works best under Linux)
- custom mini-library [md-utils](https://github.com/OOLebedenko/md-utils), which builds on pyxmolpp2 library and adds handy functions for structural analysis

**System:**
- The workflow was developed and tested within Lunix system
- Computations were run on clusters with CPUs (up to 20 CPU per MM-GBSA or FlexddG job) and GPUs (rtx_2080_ti cards, 1 card per MD simulation)

### Repository structure:  

There are two main folders with materials and results for **MP1** and **MP3** and an `md-utils` folder containing functions for analysis. Each MP folder contains 6 directories corresponding to 5 stages of the workflow. Each folder has 'scripts' subdir and its own .md instructions, by following which one can fully reproduce the analysis.

- **0_prepare_structures**. Scripts in this folder allow mutation of MP3/wild type RBD complex structure to desired variants. Mutated structures are subsequently used by FlexddG and Amber (steps 1 and 2)
- **1_FlexddG**. This folder provides instructions for launching saturational mutational scan and pair mutationnal scan. In addition scripts for output processing and visualisations are provided.
- **2_MD_Amber**. At this stage molecular dynamics jobs are prepared and launched for MP3/RBD complexes. There is also a script for trajectory conversion to .pdb for a convenient visualization in VMD
- **3_RMSD**. Here scripts for RMSD calculation and plotting are provided.
- **4_MM_GBSA**. This step provides instruction for MM-GBSA binding free energy estimation based on MD trajectories. Scripts for Amber output parsing and visualisation are provided too
- **5_MM_GBSA_decomposition**. In the same manner to the previous step, preparation and launching of MM-GBSA decomposition runs is described. A script for results parsing and differential contact map plotting is provided
