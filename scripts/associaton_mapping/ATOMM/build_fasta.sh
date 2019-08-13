#!/bin/sh
#
#  Reserve 16 CPUs for this job
#$ -pe parallel 16
#  Request 32G*4 of RAM
#$ -l h_vmem=8G
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
#$ -N build_fasta.$1
#  Run job from current working directory
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -e error_build_fasta.out
#$ -o output_build_fasta.out
#$ -j y
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#

# usage qsub -v curr_direc=/ebio/abt6_projects9/metagenomic_controlled/data/processed_reads/dc3000_infections/ /ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline/centrifuge/centrifuge_total_pipeline.sh


#centrifuge pipeline instructions
#change to the directory of reads
export PATH=/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/bin:$PATH
source /ebio/abt6/tkarasov/.bash_profile
source activate /ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping

python=/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/bin/python

$python /ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/associaton_mapping/ATOMM/run_fill_fasta.py
