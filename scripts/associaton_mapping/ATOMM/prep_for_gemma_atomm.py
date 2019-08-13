#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 16:13:49 2019

@author: tkarasov
"""
import os

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
    
    
os.chdir("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_file/")

phenotypes_all = [line.strip().split() for line in open("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/image_analysis/mean_pixels.txt").readlines()]

phenotypes = [rec for rec in phenotypes_all if rec[0] in [thing.name for thing in gene_pa]]

strain_list = [line[0] for line in phenotypes]

my_ped, my_fam, my_pheno = build_ped(phenotypes, strain_list, gene_pa)


handle = open("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/strain.fam", "w")
for line in my_fam:
	handle.write( '\t'.join([str(rec) for rec in line]))
	handle.write("\n")


#from the vcf must build 

#convert bed file to tped
cmd = ['/usr/bin/plink1', '--bfile', "gene_snp", "--recode12", "--output-missing-genotype", "0", "--transpose", "--out", "gene_snp"]
#plink --bfile [bed_prefix] (or --file [ped_prefix]) --recode12 --output-missing-genotype 0 --transpose --out [tped_prefix]
result = subprocess.run(cmd) #, stdout=subprocess.PIPE, input=input)

##plink1 --file strain --allow-no-sex --make-bed --maf 0.01 -out strain
cmd = ['plink1', '--file', 'strain', '--allow-no-sex', '--make-bed', '--maf', '0.01', '-out', 'strain']
result = subprocess.run(cmd)

#Make ped that is compatible for gemma
plink --file [file_prefix] --make-bed --out [bedfile_prefix]

#GEMMA reads four files
#genotype ped
#phenotype ped
#kinship