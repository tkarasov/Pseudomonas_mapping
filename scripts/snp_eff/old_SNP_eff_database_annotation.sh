# The goal of this script is to build a database for Pseudomonas OTU5 genomes
#http://snpeff.sourceforge.net/SnpEff_manual.html#databases

#Configure a new genome in SnpEff's config file snpEff.config.
#Add genome entry to snpEff's configuration
#!/bin/bash
#$ -N run_mapping
#
# This one was written for mapping to a difference viridiflava genome


cd /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff
#   
########################################
# Map list of Pseudomonas to strain CFBP 1590 
# https://www.ncbi.nlm.nih.gov/genome/?term=viridiflava+uasws0038
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
# Pause this next step until previous is done.
# Now build vcf: https://wurmlab.github.io/genomicscourse/2016-SIB/practicals/population_genetics/map_call,
########################################


ln -rs tmp/JRXH01.fasta tmp/reference.fa
samtools faidx tmp/reference.fa

# Run samtools mpileup
mkdir tmp/variants
bcftools mpileup -B --threads 8  --fasta-ref tmp/reference.fa \
/ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/tmp/alignments/*.bam  \
> tmp/variants/raw_calls.bcf

#mpileup -uf tmp/reference.fa  --min-MQ 60 


#cat tmp/variants/raw_calls.bcf | bcftools call -Ou -v -m - \
 #| bcftools filter -Ov -e 'QUAL<40 || DP<10 || GT!="1/1"' > tmp/variants/filtered_calls.vcf

#| bcftools norm -Ou -f "$REF" -d all - \
# Run bcftools call
bcftools call --ploidy 1 -v -m tmp/variants/raw_calls.bcf > tmp/variants/calls.vcf

bcftools filter --exclude 'QUAL < 30 || DP<10 || GT="0/1"' \
tmp/variants/calls.vcf  > tmp/variants/filtered_calls.vcf

cat tmp/variants/filtered_calls.vcf \
| sed -e 's/ENA|JRXH0100000/Contig_/g' \
| sed -e 's/ENA|JRXH010000/Contig_/g' | sed -e 's/ENA|JRXH01000/Contig_/g' \
|sed -e 's/ENA|JRXH0100/Contig_/g' | sed  's/|JRX.*,length/,length/ 
' > temp


for i in {1..188}; do
	foo=`printf "%06d" $i` 
	echo $foo
	hm="|JRXH01${foo}.1"
	cat temp | sed "s/$hm//g" > temp.2
	echo $hm
	mv temp.2  temp
done
#rm temp.2

mv temp tmp/variants/filtered_calls_renamed.vcf


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
java -Xmx4G -jar snpEff.jar dump Pseudomonas_viridiflava_gca_001305955

########################################
# Now run SNPeff
########################################
cd /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/tmp/variants
java -Xmx4G -jar /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff.jar  \
-c /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff.config \
-v Pseudomonas_viridiflava_gca_001305955  -csvStats\
 /ebio/abt6_projects8/Pseudomonas_mapping/data/SNP_eff/read_mapping/tmp/variants/filtered_calls_renamed.vcf > \
 my90_gca001305955


