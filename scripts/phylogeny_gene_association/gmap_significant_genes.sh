#!/bin/sh

dir=/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25_c2_liftover/

refgenome=/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon.fasta
#p25.C2.contigs.second_polished.pilon.fasta
gmap_build -D $dir/ -d p25_c2 $refgenome

fasta_record=/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/otu5_significant_genes.fasta

avrE=/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/presence_absence/p25.C2_avrE.fasta

output=/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/presence_absence/otu5_genes

cd /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/presence_absence/
gmap -D $dir -d p25_c2 -f samse -t 8 $fasta_record | samtools view -Shb| samtools sort -o otu5_significant_genes.bam
samtools index otu5_significant_genes.bam

gmap -D $dir -d p25_c2 -f samse -t 8 $avrE | samtools view -Shb| samtools sort -o p25_c2_avrE.bam
samtools index p25_c2_avrE.bam

#the next step is to
