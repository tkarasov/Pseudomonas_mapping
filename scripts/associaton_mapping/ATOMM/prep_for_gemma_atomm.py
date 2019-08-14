#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 16:13:49 2019

@author: tkarasov
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import glob
import sys
import pickle
import pandas as pd

def build_ped(phenotypes, strain_list, gene_pa):
	my_ped = []
	my_fam = []
	my_pheno = []
	keep = [line for line in phenotypes if line[0] in strain_list]
	genotypes = []
	for line in keep:
		my_fam.append([line[0], line[0], 0, 0, 0, line[1]])
		pa_keep = str([rec.seq for rec in gene_pa if rec.name == line[0]][0])
		pd_seq = []
		for snp in pa_keep:
			pd_seq.extend(str(int(snp) + 1)*2)
		all_dets = [line[0], line[0], 0, 0, 0, line[1]]
		all_dets.extend(pd_seq)
		my_ped.append(all_dets)
		my_pheno.append([line[0], line[0], line[1]])
	return(my_ped, my_fam, my_pheno)

def create_fasta(gc):
    gene_coordinates = []
    temp_fasta = {strain:SeqRecord.SeqRecord(seq="", id = strain) for strain in gc.columns}
    i=0
    for rec in gc.index:
       #gene_length = len(gc[STRAIN][rec])
       for strain in gc.columns:
           temp_fasta[strain].seq = Seq.Seq("".join([str(temp_fasta[strain].seq),gc[strain][rec]]))
           #gene_coordinates[rec] = [i, i+len(gc[strain][rec])]
    return temp_fasta

output_directory = sys.argv[1]#"/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/"
os.chdir(output_directory + "/temp_gc_pd")

#must merge all filled pd data frames together
filled_list = glob.glob("filled*")
all_pd = []
for rec in filled_list:
    all_pd.append(pickle.load(open(rec, 'rb')))

#now merge all dataframes
gc_filled = pd.concat([line for line in all_pd])

whole_fasta = create_fasta(gc_filled)

SeqIO.write(list(whole_fasta.values()), output_directory+"/concatenated_gene_SNP.fasta", "fasta")

#Now convert to vcf
os.chdir(output_directory)

os.system("/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/snp-sites -v -b concatenated_gene_SNP.fasta -o my_1524.vcf" )

#Now convert to plink
os.system("/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/plink --vcf my_1524.vcf --double-id --maf 0.01 --recode --out my_1524.ped")
os.system("/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/vcftools --vcf my_1524.vcf --out my_1524 --plink-tped")


#with open("snp_coordinates.txt", "w") as f:
 #   f.write(",".join(gene_coordinates))



phenotypes_all = [line.strip().split() for line in open("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/image_analysis/mean_pixels.txt").readlines()]

#phenotypes = [rec for rec in phenotypes_all if rec[0] in [thing.name for thing in gene_pa]]

#strain_list = [line[0] for line in phenotypes]

#my_ped, my_fam, my_pheno = build_ped(phenotypes, strain_list, gene_pa)


#handle = open("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/strain.fam", "w")
#for line in my_fam:
#	handle.write( '\t'.join([str(rec) for rec in line]))
#	handle.write("\n")


#from the vcf must build 

#convert bed file to tped
#cmd = ['/usr/bin/plink1', '--bfile', "gene_snp", "--recode12", "--output-missing-genotype", "0", "--transpose", "--out", "gene_snp"]
#plink --bfile [bed_prefix] (or --file [ped_prefix]) --recode12 --output-missing-genotype 0 --transpose --out [tped_prefix]
#result = subprocess.run(cmd) #, stdout=subprocess.PIPE, input=input)

##plink1 --file strain --allow-no-sex --make-bed --maf 0.01 -out strain
#cmd = ['plink1', '--file', 'strain', '--allow-no-sex', '--make-bed', '--maf', '0.01', '-out', 'strain']
#result = subprocess.run(cmd)

#Make ped that is compatible for gemma
#plink --file [file_prefix] --make-bed --out [bedfile_prefix]

#GEMMA reads four files
#genotype ped
#phenotype ped
#kinship