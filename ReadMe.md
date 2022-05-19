# Analysis and construction of SARS-CoV-2 neutralizing ligands with extensive Spike binding
*This repository contains results and materials for analysis reproduction of mini-protein ligands binding to several SARS-CoV-2 RBD variants.*
**Authors:**
- Aleksandr Kovalenko *(All-Russian Institite of Plant Protection, Saint Petersburg, Russia)*
- Xenia Sukhanova *(ITMO University, Saint Petersburg, Russia)*  

**Supervisors:** Olga Lebedenko, Nikolai Skrynnikov *(BioNMR laboratory, Saint Petersburg State University, Saint Petersburg, Russia)*

- [Background](#sec1) </br>
- [Methods](#sec2) </br>
- [System requirements](#sec3) </br>
- [Repository structure](#sec4) </br>
- [Results for MP1](#sec5) </br>
- [Results for MP3](#sec6) </br>
- [Conclusion](#sec7) </br>

<a name="sec1"></a>
### Background:
Mini-proteins 1 and 3 (MP1/3) are artificial ligands initially designed to bind the wild type SARS-CoV-2 RBD and impair its binding to human cells. However, the emergence of new covid strains hinders their effective use.  

**Goal:** This project aims to provide a workflow for the mini-pritein binding analysis of to different variants of SARS-CoV-2 RBDs and also for MPs structural adaptation to new strains.

**Objectives:**
- Analyze MP1/3 binding to several RBD variants via MD simulaition, MM-GBSA and decomposition analyses and investigate whether binding in strain-dependent
- Perform mutation of MP1/3 to improve ligand binding to recent SARS-CoV-2 variants delta+ and omicron
- Develop a collection of python scripts for structural analysis, simulations and results processing and visualization 

<a name="sec2"></a>
### Methods: 
The workflow consists of 6 steps:

![The workflow](/images/workflow.png)

0. Structures of different MP/RBD complex variants are prepared from Cryo-EM solved structures of MP1/3 with wild type RBD available from PDB. This step utilizes `Biobb` package for structure checking and mutation and `Phenix` software for complex physical relaxation prior to subsequent simulations
1. Mutated structures can be subjected to mutational scan with `FlexddG` method within `Rosetta` sortware, which allows to identify mutants with improved binding
2. For all structures MD simulations are performed with `Amber20`
3. RMSD metrics are assesed for simulated trajectories
4. Using conformations obtained from MD, MM-GBSA analysis in `Amber20` is performed to estimate binding free energies for different complexes 
5. For a detailed analysis of residue interactions contribution in binding, MM-GBSA decomposition analysis is performed

These steps utilize two modeling softwares (Rosetta and Amber20) and several methods within them to provide a panel of evidence on MP binding to RBD.

<a name="sec3"></a>
### System requirements:
**Key packages and programs:**
- the majority of scripts are written on `Python3` and there are also `bash` scripts
- `slurm` cluster management and job scheduling system
- `Rosetta` and `Amber20` (Amber assembly with MPI is highly preferrable)
- [pyxmolpp2](https://github.com/sizmailov/pyxmolpp2) python library (note thath this library works best under Linux)
- custom mini-library [md-utils](https://github.com/OOLebedenko/md-utils), which builds on pyxmolpp2 library and adds handy functions for structural analysis. Clone tmd-utils repository in this project to utilize additional functions

**System:**
- The workflow was developed and tested within Lunix system
- Computations were run on clusters with CPUs (up to 20 CPU per MM-GBSA or FlexddG job) and GPUs (rtx_2080_ti cards, 1 card per MD simulation)

<a name="sec4"></a>
### Repository structure:  

There are two main folders with materials and results for **MP1** and **MP3** and. Each MP folder contains 6 directories corresponding to the steps of the workflow. Each folder has 'scripts' subdir and its own .md instructions, by following which one can fully reproduce the analysis. The results are stored at `tables` and `results_plots` subdirs and are also discussed in .md files.

- **0_prepare_structures**. Scripts in this folder allow mutation of MP3/wild type RBD complex structure to desired variants. Mutated structures are subsequently used by FlexddG and Amber (steps 1 and 2)
- **1_FlexddG**. This folder provides instructions for launching saturational mutational scan and pair mutationnal scan. In addition scripts for output processing and visualisations are provided
- **2_MD_Amber**. At this step molecular dynamics jobs are prepared and launched for MP/RBD complexes. There is also a script for trajectory conversion to .pdb for a convenient visualization in VMD
- **3_RMSD**. Here scripts for RMSD calculation and plotting are provided
- **4_MM_GBSA**. This step provides instructions for MM-GBSA binding free energy estimation based on MD trajectories. Scripts for Amber20 output parsing and visualisation are provided too
- **5_MM_GBSA_decomposition**. In the same manner to the previous step, preparation and launching of MM-GBSA decomposition runs is described. A script for results parsing and differential contact map plotting is provided

<a name="sec5"></a>
### Results for MP1:

Within this work a panel of SARS-CoV-2 varints was analyzed: alpha, delta, delta+ and omicron.


<a name="sec6"></a>
### Results for MP3:

**MM-GBSA:**
MM-GBSA binding free energy estimations shows that MP3 binding changes differently across RBD variants (see ΔΔG bar chart below). Depending on the model chosen for MM-GBSA calculation (igb2 or 8) predictions for alpha and delta variants vary. However, a significant destabilization is observed for the complex with delta+ RDB (ΔΔG > 0) and stabilization is evident for omicron binding (ΔΔG < 0). 

![mmgbsa_mp3](/images/ddG_mmgbsa_mp3.png | width=100)

The differential contact map below, plotted from decomposition results, demostrates stabilization/destabilization of residue interactions in MP3/omicron complex in relation to the complex with wild type RDB. This chart highlights a strong interaction between Asp37 in MP3 and Arg493 in RBD, which leads to overall decrease in ΔG binding. Therefore, mutation of MP3 seems to be irrelevant in this case.

![decomposition_mp3_omicron](/images/contact_map_omicron+mp3_igb8.png | width=100)

The MP3/RBD delta+ complex was analyzed in FlexddG to design mutants with improved binding. The heatmap below shows FlexddG scores for each mutation. Mutations with scores below -1 are considered stabilising. Asp37Arg mutation appears to be one of the most stabilyzing and may lead to formation of new polar interactions. In addition to MP3(D37R), a binding of MP3(T10W;D37R) variant to delta+ was investigated by MD and MM-GBSA. The variant with two mutations was proposed based on pair mutational scan (see 1_FlexddG subdir for more details).

![mmgbsa_mp3](/images/heatmap_delta+.png | width=100)

MM-GBSA analysis for mutant MP3 varints demonstrates improved binding to delta+ RDB, although second mutation adds no significant inprovement compared to MP3(D37R) (see the ΔΔG bar chart above). Decomposition analysis for MP3(D37R)/delta+ complex shows that Arg37-Glu484 stabilizing interaction is formed, which proves the direct effect of introduced mutation. 

![mmgbsa_mp3](/images/contact_map_delta_p+mp3_d37r_igb8.png | width=100)

<a name="sec7"></a>
### Conclusions:
- The algorithm for binding analysis and mutational optimisation of MP ligands was assesed for a panel of SARS-CoV-2 RBD variants
- A collection of python scripts is provided for structure manipulation, simulations, binding analyses and results processing 
- The workflow can be further applied for newly emerging varints to facilotate design of improved inhibitors

**For MP1**:
- MP1 showed lowered binding to omicron and delta+ RBD. No improving mutations were found by FlexddG
- A structural redesign may be preferable in this case

**For MP3**
- MP3 exhibited lovered binding to delta+ variant, while binding to other variants is estimated same as to wild type or even more affine
- MP3(D37R) mutant with enhanced binding to delta+ RBD was proposed
- Omicron RBD is estimated to bind MP3 well without mutations
