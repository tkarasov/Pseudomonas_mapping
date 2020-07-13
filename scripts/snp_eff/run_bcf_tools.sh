#!/bin/sh
#$ -pe parallel 8
#  Request 32G of RAM
#$ -l h_vmem=32G
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
#$ -N map_vifi.$read1
# contains both the real output and the error messages.
#$ -e /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/out_file/error_bcf_tools.out
#$ -o /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/out_file/error_bcf_tools.out
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#
## $ -V

cd /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping
bcftools=/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/bcftools
samtools=/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/samtools  

#ln -rs tmp/p25.C2.contigs.second_polished.pilon.fasta tmp/reference.fa
# $samtools faidx tmp/reference.fa

# Run samtools mpileup
mkdir tmp/variants

$bcftools mpileup -B --threads 16  --fasta-ref tmp/sequences.fa \
/ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/tmp/alignments/*.bam  \
> tmp/variants/raw_calls2.bcf

#mpileup -uf tmp/reference.fa  --min-MQ 60 


#cat tmp/variants/raw_calls.bcf | bcftools call -Ou -v -m - \
 #| bcftools filter -Ov -e 'QUAL<40 || DP<10 || GT!="1/1"' > tmp/variants/filtered_calls.vcf

#| bcftools norm -Ou -f "$REF" -d all - \
# Run bcftools call

$bcftools call --ploidy 1 --threads 16 -v -m tmp/variants/raw_calls2.bcf > tmp/variants/calls.vcf

$bcftools filter --exclude 'AVG(MQ)<30 || DP<5 || GT="het" || AN<80' \
tmp/variants/calls.vcf  > tmp/variants/filtered_calls2.vcf

# Now make another version of the file that doesn't have B728a
vcftools --remove-indv /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/tmp/alignments/B728a \
--stdout  --non-ref-ac-any 1 --recode --recode-INFO-all --vcf tmp/variants/filtered_calls2.vcf > \
tmp/variants/filtered_calls_no_B728.vcf


#cat tmp/variants/filtered_calls.vcf | sed -e 's/utg000001c:1.0-5963307.0_pilon/utg000001l_p25.C2/g' > tmp/variants/filtered_calls_renamed.vcf
# cat tmp/variants/filtered_calls.vcf \
# | sed -e 's/ENA|JRXH0100000/Contig_/g' \
# | sed -e 's/ENA|JRXH010000/Contig_/g' | sed -e 's/ENA|JRXH01000/Contig_/g' \
# |sed -e 's/ENA|JRXH0100/Contig_/g' | sed  's/|JRX.*,length/,length/ 
# ' > temp


# for i in {1..188}; do
# 	foo=`printf "%06d" $i` 
# 	echo $foo
# 	hm="|JRXH01${foo}.1"
# 	cat temp | sed "s/$hm//g" > temp.2
# 	echo $hm
# 	mv temp.2  temp
# done
#rm temp.2

# mv temp tmp/variants/filtered_calls_renamed.vcf
