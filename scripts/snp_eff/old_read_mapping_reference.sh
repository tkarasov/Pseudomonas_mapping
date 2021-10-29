#!/bin/sh
#  Reserve 8 CPUs for this job
#$ -pe parallel 8
#  Request 32G of RAM
#$ -l h_vmem=32G
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
#$ -N map_vifi.$read1
# contains both the real output and the error messages.
#$ -e /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/out_file/error_bowtie.out
#$ -o /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/out_file/output_bowtie.out
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#
## $ -V

# read1=$read1
# read2=$read2
# plate=$plate
# pos=$pos
# path=$path

# cd /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping
# mkdir tmp
# #ln -rs /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/references/GCF_900184295.1_Chr_1_genomic.fna tmp/
# ln -rs /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/references/JRXH01.fasta tmp
# bowtie2-build tmp/JRXH01.fasta tmp/reference
# mkdir tmp/alignments
# bowtie2 \
#  -x tmp/reference \
#  -1 $path/$read1 \
#  -2 $path/$read2 \
#   > tmp/alignments/$plate.$pos.sam

#   samtools view -Sb tmp/alignments/$plate.$pos.sam | samtools sort - > tmp/alignments/$plate.$pos.bam
#   samtools index tmp/alignments/$plate.$pos.bam   # creates f1_B.sorted.bam.bai


read1=$read1
read2=$read2
plate=$plate
pos=$pos
path=$path

cd /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/
mkdir tmp
#ln -rs /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/references/GCF_900184295.1_Chr_1_genomic.fna tmp/
ln -rs /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon.fasta tmp
#bowtie2-build tmp/p25.C2.contigs.second_polished.pilon.fasta tmp/reference
bwa index tmp/p25.C2.contigs.second_polished.pilon.fasta
mkdir tmp/alignments
bowtie2 \
 -x tmp/reference \
 -1 $path/$read1 \
 -2 $path/$read2 \
  > tmp/alignments/$plate.$pos.sam

  samtools view -Sb tmp/alignments/$plate.$pos.sam | samtools sort - > tmp/alignments/$plate.$pos.bam
  samtools index tmp/alignments/$plate.$pos.bam   # creates f1_B.sorted.bam.bai
