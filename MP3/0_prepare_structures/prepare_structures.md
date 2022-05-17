 # 0. Structure preparation
 
 Requirements:
 - Biobb_structure_checking package
 - Phenix (1.19.2-4158 or other)
   
 During this step structures of mini-protein3(MP3) with different RBD variants are prepared for subsequent work. The starting point is a 7JZM PDB structure of MP3 with wild type RBD solved by cryo-EM. It is located at "structures" folder.

In order to prepare structures navigate to 'scripts' dir and execute mutate_initial_structure.sh in terminal:
```sh
cd 0_prepare_structures/scripts
bash mutate_initial_structure.sh
```
- The example is set to complex with delta+ RBD
- The script contains specifications for all mutants preparation in comments
- Several structures can be prepared at once 
- Precomputed structures, minimized with Phenix, are available at 'structures/minimized_structures'. If you want to use them - move these up to '0_prepare_structures'