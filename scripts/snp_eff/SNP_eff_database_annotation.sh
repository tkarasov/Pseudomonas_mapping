# The goal of this script is to build a database for Pseudomonas OTU5 genomes
#http://snpeff.sourceforge.net/SnpEff_manual.html#databases

#Configure a new genome in SnpEff's config file snpEff.config.
#Add genome entry to snpEff's configuration
#!/bin/bash
#$ -N run_mapping
#
# This one was written for mapping to a difference viridiflava genome


########################################
# Index reference sequence
########################################
conda activate mapping
cd /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/
mkdir tmp
#ln -rs /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/references/GCF_900184295.1_Chr_1_genomic.fna tmp/
ln -rs /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/data/Pseudomonas_viridiflava_p25_c2/sequences.fa.gz tmp
#bowtie2-build tmp/p25.C2.contigs.second_polished.pilon.fasta tmp/reference
bwa index tmp/sequences.fa.gz
mkdir tmp/alignments


cd /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff
#   
########################################
# Map list of Pseudomonas to p25.c2
########################################

#strain list of OTU5 to which we should map (taken from recombination analyses in previous study)
strains=("plate25.C11" "plate11.F10" "plate2.H5" "plate26.B10" "plate6.E4" "plate11.E1" "plate8.A7" \
	"plate7.E11" "plate11.C12" "plate3.F12" "plate1.D11" "plate8.B3" "plate9.B7" "plate7.D8" "plate12.B6" "plate6.C2"\
	"plate13.B6" "plate4.A6" "plate3.A3" "plate4.B3" "plate20.G1" "plate9.A7" "plate24.B1" "plate5.E12" \
	"plate12.D9" "plate5.D12" "plate3.A6" "plate20.D1" "plate3.D10" "plate20.C1" "plate20.E1" "plate5.D10"\
	 "plate22.A3" "plate5.A4" "plate7.D9" "plate24.H1" "plate2.C2" "plate3.C11" "plate8.C3" "plate3.B6" \
	 "plate13.H5" "plate1.A3" "plate12.G6" "plate12.G7" "plate4.G8" "plate3.F1" "plate2.H4" "plate26.G1" \
	 "plate23.B2" "plate23.A2" "plate5.B8" "plate13.H3" "plate26.B8" "plate6.F8" "plate6.B2" "plate22.D2" \
	 "plate26.G2" "plate22.E2" "plate26.D2" "plate5.F8" "plate3.C8" "plate22.G3" "plate1.H4" "plate11.H6" \
	 "plate11.H9" "plate12.B10" "plate24.H7" "plate8.G7" "plate26.D3" "plate11.G10" "plate9.B2" "plate9.D7" \
	 "plate27.F4" "plate3.F6" "plate23.C2" "plate26.F10" "plate11.A2" "plate1.A5" "plate23.A3" "plate8.H7" \
	 "plate3.B9" "plate13.F3" "plate26.H8" "plate24.G2" "plate2.D1" "plate9.C3" "plate8.E5" "plate26.B7" \
	 "plate4.H9" "plate22.A8" "plate13.C11" "plate6.H1" "plate9.G12" "plate13.D10" "plate22.C1" "plate25.B12" \
	  "plate7.F9" "plate5.H2" "plate11.H8" "plate21.C10" "plate13.H11" "plate24.G8" "plate23.B3" "plate6.A11" \
	  "plate23.A8" "plate7.E7" "plate23.D3")

WORD_LIST1=('pl1' 'pl2' 'pl6')
WORD_LIST2=('pl3' 'pl4' 'pl5')
WORD_LIST3=('pl07' 'pl08' 'pl09' 'pl11' 'pl12' 'pl13')
WORD_LIST4=('pl20' 'pl21' 'pl22' 'pl23' 'pl24' 'pl25' 'pl26' 'pl27')

for strain in "${strains[@]}" ; do
	plate=$(echo $strain | cut -f1 -d "." | sed -e 's/plate/pl/g')
	pos=$(echo $strain | cut -f2 -d ".")
	#echo $plate
	#echo $pos
if [[ " ${WORD_LIST1[@]} " =~ " ${plate} " ]];
	then curr_direc=/ebio/abt6_projects9/Pseudomonas_diversity/data/raw_reads/run_0040_0034_trimmed
	path=$curr_direc/illumina_ST-J00101_flowcellA_SampleId"$pos"__tkarasov_"$plate"_RunId0040_LaneId1/
	read1=`ls $path | grep R1` 
	read2=`ls $path | grep R2`
	qsub -v read1=$read1,read2=$read2,plate=$plate,pos=$pos,path=$path /ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/snp_eff/read_mapping_reference.sh 
fi
if [[ " ${WORD_LIST2[@]} " =~ " ${plate} " ]];
	then curr_direc=/ebio/abt6_projects9/Pseudomonas_diversity/data/raw_reads/run_0040_0034_trimmed
	path=$curr_direc/illumina_ST-J00101_flowcellA_SampleId"$pos"__tkarasov_"$plate"_*_03_2016_trimmed_RunId0034_LaneId7/ 
	read1=`ls $path| grep R1` 
	read2=`ls $path| grep R2` 
	echo "YEXXXXXXXXXXXXXXX"
	echo $read1
	qsub -v read1=$read1,read2=$read2,plate=$plate,pos=$pos,path=$path /ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/snp_eff/read_mapping_reference.sh
fi
if [[ " ${WORD_LIST3[@]} " =~ " ${plate} " ]];
	then curr_direc=/ebio/abt6_projects9/Pseudomonas_diversity/data/raw_reads/run_0041
	path=$curr_direc/illumina_ST-J00101_flowcellA_SampleId$"$pos"__tkarasov_"$plate"_RunId0041_LaneId2/
	read1=`ls $path | grep R1`
	read2=`ls $path | grep R2`
	qsub -v read1=$read1,read2=$read2,plate=$plate,pos=$pos,path=$path /ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/snp_eff/read_mapping_reference.sh
fi
if [[ " ${WORD_LIST4[@]} " =~ " ${plate} " ]] 
	then curr_direc=/ebio/abt6_projects9/Pseudomonas_diversity/data/raw_reads/run_0053
	path=$curr_direc/illumina_ST-J00101_flowcellA_SampleId"$pos"__tkarasov_"$plate"_RunId0053_LaneId*/
	read1=`ls $path | grep R1`
	read2=`ls $path | grep R2`
	qsub -v read1=$read1,read2=$read2,plate=$plate,pos=$pos,path=$path  /ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/snp_eff/read_mapping_reference.sh
fi
done

########################################
# We need one outgroup genome for pop-gen statistics. We are going to use B728a.
########################################
cd /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/references
wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR005/ERR005143/ID49_020708_20H04AAXX_R1.s_7_sequence.fastq 
read1=/ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/references/ID49_020708_20H04AAXX_R1.s_7_sequence.fastq

cd /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping
bwa mem tmp/sequences.fa.gz \
	$read1 \
 	 > tmp/alignments/B728a.sam

samtools view -Sb tmp/alignments/B728a.sam | samtools sort - > tmp/alignments/B728a.bam
samtools index tmp/alignments/B728a.bam

########################################
# Pause this next step until previous is done.
# Now build vcf: https://wurmlab.github.io/genomicscourse/2016-SIB/practicals/population_genetics/map_call,
########################################
qsub /ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/snp_eff/run_bcf_tools.sh

#bcftools view -g ^miss
########################################
# Check the validity of my edited vcf
# https://help.galaxyproject.org/t/snpeff-annotation-transcript-information-discordant-to-the-information-available-on-the-ensemble-website/2670/2
########################################
ln -rs /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/SnpSift.jar
ln -rs /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff.jar
#get chromsome names from database
java -Xmx4G -jar snpEff.jar chromosome-inf Pseudomonas_viridiflava_gca_001305955

java -Xmx4G -jar SnpSift.jar vcfCheck tmp/variants/filtered_calls_renamed.vcf

########################################
# Download genome of interest
########################################
java -Xmx4G -jar snpEff.jar  download 
java -Xmx4G -jar /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff.jar dump Pseudomonas_viridiflava_p25_c2

########################################
# Now run SNPeff
########################################
#to check the integrity of the vcf:
java -jar /ebio/abt6_projects8/Pseudomonas_mapping/Programs/jvarkit/dist/vcf2table.jar filtered_calls.vcf

cd /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/tmp/variants
java -jar /ebio/abt6_projects8/Pseudomonas_mapping/Programs/jvarkit/dist/vcf2table.jar my90_p25_c2.vcf


java -Xmx4G -jar /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff.jar  \
-c /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff_p25_c2.config \
-v Pseudomonas_viridiflava_p25_c2 /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/tmp/variants/filtered_calls2.vcf > \
 my90_B728a_p25_c2.vcf


java -Xmx4G -jar /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff.jar  \
-c /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff_p25_c2.config \
-v Pseudomonas_viridiflava_p25_c2 /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/tmp/variants/filtered_calls_no_B728.vcf > \
 my90_p25_c2.vcf

