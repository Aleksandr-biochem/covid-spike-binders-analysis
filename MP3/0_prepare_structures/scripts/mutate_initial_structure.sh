
# specify path to python environment
PATH_TO_VENV="/home/akovalenko/cov2/venv_cov2/"
# path to BioBB chech_structure.py
PATH_TO_BIOBB_STRUCTURE_CHECKING=$PATH_TO_VENV/lib/python3.8/site-packages/biobb_structure_checking/check_structure.py
# path to phenix
PATH_TO_PHENIX="/home/olebedenko/tools/phenix-1.19.2-4158/phenix_env.sh"

## make directory with temporary files
mkdir -p ../tmp

## remove ligands from original 7jzm.pdb of wild type RBD with MP3
## wt stands for wild type
python $PATH_TO_BIOBB_STRUCTURE_CHECKING -i ../structures/7jzm.pdb -o ../tmp/wt+mp3_intermediate_1.pdb --force_save ligands --remove All

## remove water
python $PATH_TO_BIOBB_STRUCTURE_CHECKING -i ../tmp/wt+mp3_intermediate_1.pdb -o ../tmp/wt+mp3_intermediate_2.pdb --force_save water --remove Yes

## fix missing heavy atoms 
python $PATH_TO_BIOBB_STRUCTURE_CHECKING -i ../tmp/wt+mp3_intermediate_2.pdb -o ../tmp/wt+mp3.pdb --force_save fixside --fix All

## mutate chains from original pdb to desired RBD and MP mutant
## an example is set to alpha variant RBD
python $PATH_TO_BIOBB_STRUCTURE_CHECKING -i ../tmp/wt+mp3.pdb -o ../tmp/delta_plus+mp3.pdb mutateside --mut B:Lys417Asn,B:Leu452Arg,B:Thr478Lys

#### below is a list of -mut specificastions for each covid and MP3 variant
# alpha                     B:Asn501Tyr
# delta 					B:Leu452Arg,B:Thr478Lys
# delta_plus 				B:Lys417Asn,B:Leu452Arg,B:Thr478Lys
# omicron 				    B:Lys417Asn,B:Thr478Lys,B:Glu484Ala,B:Asn501Tyr,B:Gly339Asp,B:Ser371Leu,B:Ser373Pro,B:Ser375Phe,B:Asn440Lys,B:Gly446Ser,B:Ser477Asn,B:Gln493Arg,B:Gly496Ser,B:Gln498Arg,B:Tyr505His
# delta_plus+MP3_D37R 	    A:Asp37Arg,B:Lys417Asn,B:Leu452Arg,B:Thr478Lys
# delta_plus+MP3_T10W_D37R  A:Asp37Arg,A:Thr10Trp,B:Lys417Asn,B:Leu452Arg,B:Thr478Lys
#### 

## minimize each mutated structure with phenix
source $PATH_TO_PHENIX
phenix.geometry_minimization ../tmp/wt+mp3.pdb write_geo_file=False directory=../
phenix.geometry_minimization ../tmp/delta_plus+mp3.pdb write_geo_file=False directory=../

# the output pdb files are a starting point for FlexddG or MD simulation 