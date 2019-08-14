#!/bin/bash
#$ -S /bin/bash
#
#  Reserve 1 CPU for this job
#$ -l h_vmem=32G
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
#$ -N prep_for_gemma.$1
#  Run job from current working directory
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -e error_prep_gemma.out
#$ -o output_prep_gemma.out
#$ -j y
#$ -cwd
#  Import environment
#$ -V

export PATH=/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/bin:$PATH
source ~/.bashrc
source activate /ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping

output_direc=$1
echo "Output directory in prep_for_gemma_atomm is:"$output_direc
code_direc=/ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/associaton_mapping/ATOMM/


python=/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/python

#the run_fill_fasta will pull together the many pandas dataframes and fills them and outputs to fasta file
$python $code_direc/prep_for_gemma_atomm.py $output_direc
