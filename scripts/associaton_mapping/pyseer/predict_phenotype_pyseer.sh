#!/bin/bash

# The goal of this script is to see if we can predict the phenotype of a strain from its genotype
# https://pyseer.readthedocs.io/en/master/predict.html

conda activate mapping
cd /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/pyseer
genome_loc=/ebio/abt6_projects9/Pseudomonas_diversity/data/deprecated_pseudomonas_genome_assembly/assembly_annotation


###########################################
# Run pyseer as enet
# https://pyseer.readthedocs.io/en/master/usage.html#processing-k-mer-output
###########################################

cd /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/pyseer/

#this is running pyseer with unitigs
pyseer --phenotypes ./input/strain_phenot.txt --kmers ./input/input_genomes/unitigs_output/unitigs.txt.gz  --wg enet \
 --save-vars ./output/ma_snps --save-model ./output/unitig.lasso --cpu 4 --alpha 1  > ./output/selected.txt



###########################################
# Split dataset in two for regression
###########################################
head -47 ./input/strain_phenot.txt > ./input/train.pheno
cat <(head -1 ./input/train.pheno) <(tail -46 ./input/strain_phenot.txt) > ./input/test.pheno
cut -f 1 ./input/test.pheno | sed '1d' > ./input/test.samples


###########################################
# Run pyseer as enet
###########################################
pyseer --kmers ./input/input_genomes/unitigs_output/unitigs.txt.gz --phenotypes ./input/train.pheno --wg enet \
--load-vars ./output/ma_snps --alpha 1 --save-model ./output/test_lasso --cpu 4

###########################################
# Make predictions
###########################################
enet_predict --kmers ./input/input_genomes/unitigs_output/unitigs.txt.gz  --true-values ./input/test.pheno \
./output/test_lasso.pkl ./input/test.samples > ./output/test_predictions.txt
