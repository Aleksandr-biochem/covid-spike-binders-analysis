 # 0. Structure preparation
 
 Requirements:
 - Biobb_structure_checking package
 - Phenix (1.19.2-4158 or other)
   
 During this step structures of LCB1 with different RBD variants are prepared for subsequent work.  The starting point is a 7JZU PDB structure of LCB with wild type RBD solved by cryo-EM. It is located at "structures" folder.

In order to prepare structures execute mutate.sh in terminal:
```sh
bash mutate.sh
```
Note that script contains specifications for all mutants preparation. Several structures can be prepared at once. By default structures will be saved in current working directory.

All expected resulting structures, minimized with Phenix, are stored at 'structures/minimized'.