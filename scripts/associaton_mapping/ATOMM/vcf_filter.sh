#!/bin/bash

#text stolen from here: https://github.com/johanzi/gwas_gemma#section-id-50

#Which accessions to keep?
vcftools --keep list_accessions_to_keep.txt --gzvcf file.vcf.gz --recode recode-INFO-all --out subset_80

#remove singletons
vcftools --singletons --vcf subset_80_biallelic_only_alt_no_indels.recode.vcf
vcftools --vcf subset_80_biallelic_only_alt_no_indels.recode.vcf \
            --exclude-positions out.singletons --recode --recode-INFO-all \
            --out subset_80_biallelic_only_alt_no_indels_no_singletons

#Compress and tabix the file
bgzip  subset_80_biallelic_only_alt_no_indels_no_singletons.recode.vcf && \
            tabix subset_80_biallelic_only_alt_no_indels_no_singletons.recode.vcf.gz

# Get list of accessions in vcf file
 bcftools query -l  subset_80.recode.vcf.gz > order_accession.txt
