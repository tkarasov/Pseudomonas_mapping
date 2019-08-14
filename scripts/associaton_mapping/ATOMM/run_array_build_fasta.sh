#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#  Import environment
#$ -V

temp_direc=$1
echo "Task id is $SGE_TASK_ID"
echo "temp_direc in run_array build is:"$temp_direc

bash /ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/associaton_mapping/ATOMM/build_fasta.sh  $temp_direc $SGE_TASK_ID
