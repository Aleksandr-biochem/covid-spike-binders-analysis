# Analysis and construction of SARS-CoV-2 neutralizing ligands with extensive Spike binding
*This repository contains the workflow and results on the analysis of mini-protein ligangs affinity to receptor-binding domains (RBD) of SARS-CoV-2 variants*

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
- [Conclusionы](#sec7) </br>
- [Literature](#sec8) </br>

<a name="sec1"></a>
### Background:
Protein minibinders (MP1 and MP3) have been designed *in silico* against RBD of wild-type spike protein to prevent SARS-CoV-2 entry into cells. However, the emergence of new covid strains hinders MPs effective use.  

**Goal:** The main goal of this project is to develop a workflow for estimation of the MP1 and MP3 proteins binding affinity to RBD of the new SARS-CoV-2 variants (alpha, delta, delta+, omicron) and propose a way to optimize the MP1 and MP3 sequences for stronger interaction with new covid strains RBDs.

**Objectives:**
- Estimate the binding free energy (ΔG binding) of MP1 and MP3 to RBDs of the new SARS-CoV-2 strains by MM-GBSA protocol (based on `MMPBSA.py` from `Amber20`)
- Perform the pairwise residue energy decomposition to increase understanding of recognition processes (also with `MMPBSA.py`)
- Optimize MP1 and MP3 sequences for stronger interaction with RBDs of the new SARS-CoV-2 strains (based on `Flexx ddG` from `Rosetta`)
- Develop a collection of python scripts for structural analysis, simulations, results processing and visualization 

<a name="sec2"></a>
### Methods: 
The workflow consists of 6 steps:

<p align="center">
  <img src="/images/workflow.png">
</p>

These steps mainly employ two modeling softwares (`Rosetta` and `Amber20`) and several methods within them to provide a panel of detailed evidence on MP binding to RBD.

0. **Structure preparation.** We use the starting coordinates of MP1 and MP3 in complex with RBD of SARS-CoV2 wild-type from PDB that have medium resolution (cryo-EM structure with ~3 Å resolution). To produce a models of other SARS-CoV-2 RDBs (alpha, delta, delta+, omicron), we performed *in silico* mutagenesis and structure ckecking with `Biobb` package ([Bayarri G. et al 2022](https://doi.org/10.1093/bioinformatics/btac316)). Finally, the complexes with mutated RBDs were relaxed with `Phenix` software package suitable for work with electron cryo-microscopy structures ([Liebschner D. et al 2019](https://doi.org/10.1107/S2059798319011471))
1. ***In silico* scanning saturation mutagenesis of MP1 and MP3.** We use the `Flex ddG` protocol from `Rosetta`, which is an effective approach to model changes in binding free energies upon mutations in complex ([Barlow K. A. et al 2018](https://doi.org/10.1021/acs.jpcb.7b11367))
2. **Molecular dynamics (MD) modelling.** All trajectories were recorded with `Amber20` package, a widely used tools for molecular simulations ([Case D. A et al 2020](https://ambermd.org/doc12/Amber20.pdf)). The length of each MD trajectory is 1.5 μs, giving a reliable representation of complexes conformational states
4. **Assessing the precision of MD models.** The RMSD profiles were calculated using in-house python script.
5. **Estimation of the binding free energy.** For this purpose, MM-GBSA energy components of the complexes were calculated for 1000 equidistant snapshots extracted from the last 1 μs of MD jrajectory. We have utilized `MMPBSA.py` from `Amber20` as a modern and reliable tool for MM-GBSa analysis ([Miller III B. R. et al 2012](https://doi.org/10.1021/ct300418h))
6. **The pairwise residue energy decomposition.** For a detailed analysis of residue interactions contribution in binding, MM-GBSA decomposition analysis is performed with `MMPBSA.py` from `Amber20`

<a name="sec3"></a>
### System requirements:
**Key packages and programs:**
- the majority of scripts are written on `Python3` and there are also `bash` scripts
- `slurm` (20.11.8) cluster management and job scheduling system
- `Rosetta2`
- `Amber20` (a build with MPI is employed and is highly preferrable)
- [amber-runner](https://github.com/sizmailov/amber-runner) (0.0.8)
- [pyxmolpp2](https://github.com/sizmailov/pyxmolpp2) (1.6.0, note that this library works best under Linux)
- [biobb-structure-checking](https://github.com/bioexcel/biobb_structure_checking) (3.9.4)
- [Phenix](http://www.phenix-online.org) (1.19.2)
- in-house mini-library [md-utils](https://github.com/OOLebedenko/md-utils), which builds on `pyxmolpp2` library and adds handy functions for structural analysis. Clone `md-utils` repository into this project to utilize its additional functions
- other python libraries used for plotting and analysis are listed in requirements.txt

**System:**
- The workflow was developed and tested within Lunix system
- Computations were run on clusters with CPUs (up to 20 CPU per MM-GBSA or Flex ddG job) and GPUs (rtx_2080_ti cards, 1 card per MD simulation)

<a name="sec4"></a>
### Repository structure:  

There are two main folders with materials and results for **MP1** and **MP3**. Each MP folder contains 6 directories corresponding to the steps of the workflow. Each folder has 'scripts' subdir and its own .md instructions, by following which one can fully reproduce the analysis. The results are stored at `tables` and `results_plots` subdirs and are also discussed in .md files.

- **0_prepare_structures**. Scripts in this folder allow mutation of MP/wild-type RBD complex structure to desired variants. Mutated structures are subsequently used by Flex ddG and Amber20 (steps 1 and 2)
- **1_FlexddG**. This folder provides instructions for launching saturational mutational scan and pair mutationnal scan. In addition scripts for output processing and visualisations are provided
- **2_MD_Amber**. At this step molecular dynamics jobs are prepared and launched for MP/RBD complexes. There is also a script for trajectory conversion to .pdb for a convenient visualization in VMD
- **3_RMSD**. Here scripts for RMSD calculation and plotting are provided
- **4_MM_GBSA**. This step provides instructions for MM-GBSA binding free energy estimation based on MD trajectories. Scripts for output parsing and visualisation are provided too
- **5_MM_GBSA_decomposition**. In the same manner to the previous step, preparation and launching of MM-GBSA decomposition runs is described. A script for results parsing and differential contact map plotting is provided

<a name="sec5"></a>
### Results for MP1:

**MM-GBSA:**

Before modification of miniproteins its necessity was checked. This was obtained by MM-GBSA free binding energy estimations for complexes of MP1 with mutant RBD in comparison to wt complex. According to these result, MP1 shows worse affinity to all RBD variants then to wt regardless solvant model (see ΔΔG bar charts below). Among all variants significant decrease in free binding energy is noticed for delta+ and omicron RBDs, making the modification of MP1 highly important.

<p align="center">
  <img src="/images/MP1_dG_GBSA_igb8.png" width="600">
</p> 

The decomposition results below demonstrate that the revealed increased destabilization of omciron and delta+ complexes is mostly caused by local interactions. This detabilization is enhanced due to Lys 417 to Asn 417 substitution in RBD, which, probably, lead to the disappearance of important electrostatic interactions (salt bridges and hydrogen bonds) between 417 residue of RBD and Asp30 of MP1. As the destabiliation pattern has local behaviour, modification of MP1 is believed to solve the affinity issue.

<p align="center">
  <img src="/images/omicron_vs_wt_igb8.png" width="800">
</p> 

<p align="center">
  <img src="/images/delta_plus_vs_wt_igb8.png" width="800">
</p> 

And as it can be noticed via structural visualization, indeed, existing in wt complex Asp30-Lys317 MP1-RBD interactions is disabled in the case of omicron complex variant.

<p align="center">
  <img src="/images/wt_mp1.png" width="700">
</p> 

<p align="center">
  <img src="/images/omicron_mp1.png" width="700">
</p> 

To improve the binding of MP1 the MP1/RBD delta+ complex was analyzed in FlexddG to design MP1 mutants. The heatmap below shows FlexddG scores for each mutation. Mutations with scores below -1 are considered stabilising. Practically all mutations appear to be destabilizing. The only one leading to the increase of free binding energy is Glu3Trp mutation. Though it may lead to formation of new interactions, it also causes the introduction of sufficient aromatic group. Thus could disrupt the MP1 and complex structure. According to these findings, structural remodelling rather modification of MP1 is preferable to improve its binding to RBD.

<p align="center">
  <img src="/images/delta_plus_flexddG.png" width="850">
</p>

<a name="sec6"></a>
### Results for MP3:

**MM-GBSA:**
MM-GBSA binding free energy estimation demonstrates that MP3 binding changes differently across RBD variants (see ΔΔG bar chart below). Complexes with alpha and delta RBDs demonstrate an increased affinity an minor changes in ΔG respectively. A significant destabilization is observed for the complex with delta+ RDB (ΔΔG > 0) and stabilization is evident for omicron binding (ΔΔG < 0) in comparison to complex with wild-type RBD. 

<p align="center">
  <img src="/images/ddG_mmgbsa_mp3.png" width="600">
</p>

The differential contact map below, plotted from decomposition results, demostrates stabilization/destabilization of residue interactions in MP3/omicron complex in relation to the complex with wild type RDB. This chart highlights a strong interaction between Asp37 in MP3 and Arg493 in RBD, which leads to overall decrease in ΔG binding. Therefore, mutation of MP3 appears to be irrelevant in this case.

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

MM-GBSA analysis for mutant MP3 variants demonstrates improved binding to delta+ RDB, although second mutation adds no significant inprovement compared to MP3(D37R) (see the ΔΔG bar chart above). Decomposition analysis for MP3(D37R)/delta+ complex shows that Arg37-Glu484 stabilizing interaction is formed, which proves the direct effect of introduced mutation. 

<p align="center">
  <img src="/images/contact_map_delta_p+mp3_d37r_igb8.png" width="700">
</p>

<p align="center">
  <img src="/images/New_bond.png" width="700">
</p>

<a name="sec7"></a>
### Conclusions:
- The workflow for the analysis of the MP1 and MP3 ligands binding affinity to RBDs of the new SARS-CoV-2 strains (alpha, delta, delta+, omicron) is developed 
- Steps for the further optimization of MP1 and MP3 sequences are suggested 
- The collection of python scripts is provided for structure manipulation, simulations, binding analyses and results processing 

**For MP1**:
- The binding of MP1 to RBD of SARS-CoV-2 omicron and delta+ is significantly weakened compared to RBD of SARS-CoV-2 wild-type 
- No improving mutations are found by FlexddG for these RDB variants. The structural redesign may be preferable in this case

**For MP3**
- MP3 exhibits lovered binding to delta+ variant, while binding to other variants is estimated same as to wild type or even more affine
- MP3(D37R) mutant with enhanced binding to delta+ RBD is proposed
- Omicron RBD is estimated to bind MP3 well without mutations

<a name="sec8"></a>
### Literature
Bayarri, G., Andrio, P., Orozco, M., & Gelpí, J. L. (2022). BioExcel Building Blocks REST API (BioBB REST API), programmatic access to interoperable biomolecular simulation tools. Bioinformatics.

Liebschner, D., Afonine, P. V., Baker, M. L., Bunkóczi, G., Chen, V. B., Croll, T. I., ... & Adams, P. D. (2019). Macromolecular structure determination using X-rays, neutrons and electrons: recent developments in Phenix. Acta Crystallographica Section D: Structural Biology, 75(10), 861-877. 

Barlow, K. A., Ó Conchúir, S., Thompson, S., Suresh, P., Lucas, J. E., Heinonen, M., & Kortemme, T. (2018). Flex ddG: Rosetta ensemble-based estimation of changes in protein–protein binding affinity upon mutation. The Journal of Physical Chemistry B, 122(21), 5389-5399.

D.A. Case, K. Belfon, I.Y. Ben-Shalom, S.R. Brozell, D.S. Cerutti, T.E. Cheatham, III, V.W.D. Cruzeiro, T.A. Darden, R.E. Duke, G. Giambasu, M.K. Gilson, H. Gohlke, A.W. Goetz,R Harris, S. Izadi, S.A. Iz- mailov, K. Kasavajhala, A. Kovalenko, R. Krasny, T. Kurtzman, T.S. Lee, S. LeGrand, P. Li, C. Lin, J. Liu, T. Luchko, R. Luo, V. Man, K.M. Merz, Y. Miao, O. Mikhailovskii, G. Monard, H. Nguyen, A. Onufriev, F. Pan, S. Pantano, R. Qi, D.R. Roe, A. Roitberg, C. Sagui, S. Schott-Verdugo, J. Shen, C.L. Simmerling, N.R. Skrynnikov, J. Smith, J. Swails, R.C. Walker, J. Wang, L. Wilson, R.M. Wolf, X. Wu, Y. Xiong, Y. Xue, D.M. York and P.A. Kollman (2020), AMBER 2020, University of California, San Francisco.

Miller III, B. R., McGee Jr, T. D., Swails, J. M., Homeyer, N., Gohlke, H., & Roitberg, A. E. (2012). MMPBSA. py: an efficient program for end-state free energy calculations. Journal of chemical theory and computation, 8(9), 3314-3321.
