# Analysis and construction of SARS-CoV-2 neutralizing ligands with extensive Spike binding
*This repository contains the analysis workflow and results regarding affinity of protein minibinders to receptor-binding domain (RBD) of SARS-CoV-2 variants*

**Authors:**
- Aleksandr Kovalenko *(All-Russian Institute of Plant Protection, Saint Petersburg, Russia)*
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
Protein minibinders (MP1 and MP3) have been designed *in silico* against RBD of wild-type spike protein to prevent  SARS-CoV-2 entry into cells. However, the emergence of new covid strains hinders their effective use.  

**Goal:** The main goal of this project is to develop the workflow for estimation of the binding affinity of MP1 and MP3 proteins to RBD of new SARS-CoV-2 variants (alpha, delta, delta+, omicron). The additional objective is to optimize the  MP1 and MP3 sequences for stronger interaction with RBD of the new SARS-CoV-2 strains.

**Objectives:**
- Estimate the binding free energy of MP1 and MP3 to RBD of the new SARS-CoV-2 strains by MM-GBSA protocol 
- Perform the pairwise residue energy decomposition to gain better understanding of recognition processes (based on MMPBSA.py from Amber20)
- Optimize MP1 and MP3 sequences for stronger interaction with RBD of the new SARS-CoV-2 strains (based on Flexx ddG from Rosetta)
- Develop a collection of python scripts for structural analysis, simulations and results processing and visualization 

<a name="sec2"></a>
### Methods: 
The workflow consists of 6 steps:

<p align="center">
  <img src="/images/workflow.png">
</p>

0. **Structure preparation.** We use the starting coordinates of MP1 and MP3 in complex with RBD of wild-type SARS-CoV2 from PDB that have medium resolution (cryo-EM structure with ~3 Å resolution). To produce a model of RBD of SARS-CoV-2 (alpha, delta, delta+, omicron), we performed *in silico* mutagenesis with `Biobb` package. Finally, the complexes with mutated RBD were relaxed by using `Phenix`
1. ***In silico* scanning saturation mutagenesis.** We use the `Flex ddG` protocol from `Rosetta` to model changes in binding free energies upon mutations in MP1 or MP3 in complex with RBD
2. **Molecular dynamics (MD) modelling.** All trajectories were recorded with `Amber20` package. The length of each MD trajectory is 1.5 μs
3. **Assessing the precision of MD models.** The RMSD profiles were calculated using in-house python script
4. **Estimation of the binding free energy.** For this purpose, MM-PBSA energy components of the complexes were calculated for 1000 equidistant snapshots extracted from the last 1 μs using `MMPBSA.py` from `Amber20`
5. **The pairwise residue energy decomposition**. For a detailed analysis of residue interactions contribution in binding, MM-GBSA decomposition analysis is performed (also with `MMPBSA.py`)

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
MM-GBSA binding free energy estimation demonstrates that MP3 binding changes differently across RBD variants (see ΔΔG bar chart below). Depending on the model chosen for MM-GBSA calculation (igb2 or 8) predictions for alpha and delta variants vary. However, a significant destabilization is observed for the complex with delta+ RDB (ΔΔG > 0) and stabilization is evident for omicron binding (ΔΔG < 0). 

<p align="center">
  <img src="/images/ddG_mmgbsa_mp3_igb8.png" width="600">
</p>

The differential contact map below, plotted from decomposition results, demostrates stabilization/destabilization of residue interactions in MP3/omicron complex in relation to the complex with wild type RDB. This chart highlights a strong interaction between Asp37 in MP3 and Arg493 in RBD, which leads to overall decrease in ΔG binding. Therefore, mutation of MP3 seems to be irrelevant in this case.

<p align="center">
  <img src="/images/contact_map_omicron+mp3_igb8.png" width="700">
</p>

<p align="center">
  <img src="/images/Omicron.png" width="700">
</p>

The MP3/RBD delta+ complex was analyzed in FlexddG to design mutants with improved binding. The heatmap below shows FlexddG scores for each mutation. Mutations with scores below -1 are considered stabilising. Asp37Arg mutation appears to be one of the most stabilyzing and may lead to formation of new polar interactions. In addition to MP3(D37R), a binding of MP3(T10W;D37R) variant to delta+ was investigated by MD and MM-GBSA. The variant with two mutations was proposed based on pair mutational scan (see 1_FlexddG subdir for more details).

<p align="center">
  <img src="/images/heatmap_delta+.png" width="700">
</p>

MM-GBSA analysis for mutant MP3 varints demonstrates improved binding to delta+ RDB, although second mutation adds no significant inprovement compared to MP3(D37R) (see the ΔΔG bar chart above). Decomposition analysis for MP3(D37R)/delta+ complex shows that Arg37-Glu484 stabilizing interaction is formed, which proves the direct effect of introduced mutation. 

<p align="center">
  <img src="/images/contact_map_delta_p+mp3_d37r_igb8.png" width="700">
</p>

<p align="center">
  <img src="/images/New_bond.png" width="700">
</p>

<a name="sec7"></a>
### Conclusions:
- The workflow for the analysis of binding affinity of MP1 and MP3 ligands to RBD of the new SARS-CoV-2 strains (alpha, delta, delta+, omicron) is developed
- Steps for the further optimization of MP1 and MP3 sequences are suggested 
- The collection of python scripts is provided for structure manipulation, simulations, binding analyses and results processing

**For MP1**:
- The binding of MP1 to RBD of SARS-CoV-2 omicron and delta+ is significantly weakened compared to RBD of SARS-CoV-2 wild-type 
- No improving mutations are found by FlexddG. The structural redesign may be preferable in this case

**For MP3**
- MP3 exhibits lovered binding to delta+ variant, while binding to other variants is estimated same as to wild type or even more affine
- MP3(D37R) mutant with enhanced binding to delta+ RBD is proposed
- Omicron RBD is estimated to bind MP3 well without mutations
