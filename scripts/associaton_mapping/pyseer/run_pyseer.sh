#!/bin/bash

# This script attempts to use pyseer to map the Pseudomonas results using pyseer
#https://pyseer.readthedocs.io/en/master/usage.html#annotating-k-mers

conda activate mapping
cd /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/pyseer
genome_loc=/ebio/abt6_projects9/Pseudomonas_diversity/data/deprecated_pseudomonas_genome_assembly/assembly_annotation
#pull all relevant genomes

for rec in `cat ./input/strain_phenot.txt | cut -f1`;
do
    plate=$(echo $rec | cut -f1 -d '.')
    pos=$(echo $rec | cut -f2 -d '.')
    echo $plate
   echo $pos
#    cp $genome_loc/$plate/$pos/*.contigs_renamed.fasta ./input/input_genomes/
#    cp $genome_loc/$plate/$pos/annotation/*.gff ./input/input_gff/
#	  my_gff=$genome_loc/$plate/$pos/annotation/*.gff
#	  my_genome=$genome_loc/$plate/$pos/*.contigs_renamed.fasta
	  printf "%s\t%s\t%s\n" ./input/input_genomes/$plate.$pos.pilon.contigs_renamed.fasta  ./input/input_gff/$plate.$pos.annotated.gff "draft" >> ./input/my_references.txt
 done

#make references.txt file


###########################################
# Count kmers in all my genomes with fsm-lite
# https://github.com/johnlees/seer/wiki/Tutorial#count
###########################################

#paste <(ls assemblies | cut -d "." -f 1) <(ls assemblies/*.fa) > fsm_files.txt
cd ./input/input_genomes
fsm-lite -v -l fsm_files.txt -t tmp_idx -s 2  -S 90 > fsm_kmers.txt
cp fsm_kmers.txt temp_fsm_kmers.txt
sed -i -e 's/.pilon.contigs_renamed.fasta//g' fsm_kmers.txt
split -d -n l/16 fsm_kmers.txt fsm_out
rm fsm_kmers.txt
gzip fsm_out*
gzip fsm_kmers.txt

###########################################
# Count unitigs
# https://github.com/johnlees/seer/wiki/Tutorial#count
###########################################
cd ./input/input_genomes
awk -F'|' '{print $0"\t"$0}' fsm_files.txt | awk 'BEGIN{FS=OFS="\t"} {gsub(/\.pilon.contigs_renamed.fasta/, "", $1)} 1' > unitig_files.txt
unitig-counter -strains unitig_files.txt -nb-cores 6 -out unitigs_output
gzip unitig_files.txt
###########################################
# Calculate kinship matrix
###########################################

mash sketch -s 10000 -o ../samples *.fasta
mash dist ../samples.msh ../samples.msh > poo
sed -i -e 's/.pilon.contigs_renamed.fasta//g' poo
awk 'BEGIN {OFS=FS="\t"} {gsub(/\./,"_",$1);gsub(/\./,"_",$2)}1' poo > poo2

cat poo2 | square_mash > ../mash.tsv
sed -i -e 's/_/./g' ../mash.tsv

###########################################
# Run pyseer as a fixed effect model using the mash distance matrix
# https://pyseer.readthedocs.io/en/master/usage.html#processing-k-mer-output
###########################################

cd /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/pyseer

#this is running pyseer with unitigs
pyseer --phenotypes ./input/strain_phenot.txt --kmers \
./input/input_genomes/unitigs_output/unitigs.txt.gz --distances ./input/mash.tsv \
--min-af 0.01 --max-af 0.99 --cpu 2  > ./output/thread2.assoc #./output/unitigs_pyseer.assoc

###########################################
# Filter on pvalue
###########################################
head -1 ./output/unitigs_pyseer.assoc  > ./output/unitigs_6_pyseer.assoc
awk '($3 + 0) < 1E-7' ./output/unitigs_pyseer.assoc  | awk '($2 + 0) < 8.5E-1' | awk '($2 + 0) > 1.5E-1'>> ./output/unitigs_6_pyseer.assoc


###########################################
# Mapping kmers to one reference
###########################################
phandango_mapper ./output/unitigs_6_pyseer.assoc /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon.fasta ./output/p25.c2_reference_1.plot
#./input/input_genomes/plate25.C2.pilon.contigs_renamed.fasta ./output/p25.c2_reference_1.plot

###########################################
# Annotating k-mers
###########################################
annotate_hits_pyseer ./output/unitigs_pyseer.assoc ./input/my_references.txt ./output/kmer_annotation.txt

python /ebio/abt6_projects8/Pseudomonas_mapping/Programs/pyseer/scripts/summarise_annotations.py kmer_annotation.txt
