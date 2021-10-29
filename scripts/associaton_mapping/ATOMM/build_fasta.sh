#!/bin/bash
#$ -S /bin/bash
#
#  Reserve 1 CPU for this job
#$ -l h_vmem=8G
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
#$ -N build_fasta.$2
#  Run job from current working directory
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -e error_build_fasta.out
#$ -o output_build_fasta.out
#$ -j y
#$ -cwd
#  Import environment
#$ -V

export PATH=/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/bin:$PATH
source ~/.bashrc
source activate /ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping

temp_direc=$1
pd=$2
echo "The file is:"$pd
echo "temp_direc in build_fasta is:"$temp_direc
path_to_pd=$temp_direc"/temp_pd"$pd".cpk"


python=/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/bin/python

#the run_fill_fasta will pull together the many pandas dataframes and fills them and outputs to fasta file
$python /ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/associaton_mapping/ATOMM/run_fill_fasta.py $path_to_pd
