PATH_TO_VENV="/home/xenia/cov2/venv/"
PATH_TO_BIOBB_STRUCTURE_CHECKING=$PATH_TO_VENV/lib/python3.8/site-packages/biobb_structure_checking/check_structure.py
PATH_TO_PHENIX="/home/olebedenko/tools/phenix-1.19.2-4158/phenix_env.sh"

mkdir -p tmp

#remove ligands
python $PATH_TO_BIOBB_STRUCTURE_CHECKING -i ../structures/7jzu.pdb -o ../tmp/wt+lcb1_intermediate_1.pdb --force_save ligands --remove All

#remove water
python $PATH_TO_BIOBB_STRUCTURE_CHECKING -i ../tmp/wt+lcb1_intermediate_1.pdb -o ../tmp/wt+lcb1_intermediate_2.pdb --force_save water --remove Yes

#fix missing heavy atoms 
python $PATH_TO_BIOBB_STRUCTURE_CHECKING -i ../tmp/wt+lcb1_intermediate_2.pdb -o ../tmp/wt+lcb1.pdb --force_save fixside --fix All

#mutate alpha, delta, delta_plus
##alpha
python $PATH_TO_BIOBB_STRUCTURE_CHECKING -i ../tmp/wt+lcb1.pdb -o ../tmp/alpha+lcb1.pdb mutateside --mut B:Asn501Tyr

# alpha                     B:Asn501Tyr
# delta 					B:Leu452Arg,B:Thr478Lys
# delta_plus 				B:Lys417Asn,B:Leu452Arg,B:Thr478Lys
# omicron 				    B:Lys417Asn,B:Thr478Lys,B:Glu484Ala,B:Asn501Tyr,B:Gly339Asp,B:Ser371Leu,B:Ser373Pro,B:Ser375Phe,B:Asn440Lys,B:Gly446Ser,B:Ser477Asn,B:Gln493Arg,B:Gly496Ser,B:Gln498Arg,B:Tyr505His
# delta_plus+LCB3_D37R 	    A:Asp37Arg,B:Lys417Asn,B:Leu452Arg,B:Thr478Lys
# delta_plus+LCB3_T10W_D37R A:Asp37Arg,A:Thr10Trp,B:Lys417Asn,B:Leu452Arg,B:Thr478Lys

# minimize structure with phenix
source $PATH_TO_PHENIX
phenix.geometry_minimization ../tmp/wt+lcb1.pdb write_geo_file=False directory=./
phenix.geometry_minimization ../tmp/alpha+lcb1.pdb write_geo_file=False directory=./


