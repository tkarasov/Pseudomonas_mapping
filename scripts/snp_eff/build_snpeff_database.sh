#!/bin/sh
#build a custom snpEff database
#https://www.biostars.org/p/359661/
cd /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/data
mkdir Pseudomonas_viridiflava_p25_c2
cd Pseudomonas_viridiflava_p25_c2
cp /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon.gff .
#cp /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon.fasta .  

#take the genome directly from the gff
zcat genes.gff.gz  | sed -n -e '/>utg000001l_p25.C2/,$p' > p25_gff_fasta.fa

#make multi-line fasta into single line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  < p25_gff_fasta.fa  > temp

#Now make into 
mv temp sequneces.fas

mv p25.C2.contigs.second_polished.pilon.gff genes.gff 
gzip genes.gff
mv p25.C2.contigs.second_polished.pilon.fasta sequences.fa
gzip sequences.fa 


cd /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/
java -jar snpEff.jar build -gff3 -v Pseudomonas_viridiflava_p25_c2 -c ./snpEff_p25_c2.config