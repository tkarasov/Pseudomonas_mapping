#!/usr/bin/python3 -*- coding: utf-8 -*-
"""
This script takes the output from panX gene clustering and rewrites fasta files and a SNP matrix that includes SNP information and P/A information.
The original program is here: https://github.com/Miaoyanwang/kinship_ATOMM
must be in conda environment mapping to have access to Bio and such
"""

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
from multiprocess import Process, Manager
import multiprocess as mp
import sys
import os
import pickle
#from subprocess import Popen, PIPE
from joblib import Parallel, delayed
import pandas as pd
from collections import ChainMap
import numpy as np
import copy

fna="GC00000001_100.fna"
def process_cluster(fna, gc_dict):
    i=0
    gene = fna.strip(".fna")
    for sequence in list(SeqIO.parse(gene_cluster_path+ "/"+ fna, "fasta")):
        i=i+1
        #print(i)
        strain = sequence.id.split("|")[0]
        gc_dict[(strain,gene)] = sequence.seq
        return(gc_dict)

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

        
def fill_missing(gc_pd):
    i = 0
    gc_filled = copy.deepcopy(gc_pd)
    gc_filled.ix[gene]
    for gene in gc_filled.index:
        fill_missing_sub(gene)
        print(gene)
    #processed_list = Parallel(n_jobs = 16)(delayed(fill_missing_sub)(gene) for gene in list(gc_filled.index))

def fill_missing_sub(gene):
#   for gene in gc_pd.index:
    gene_length = max([len(gc_pd[rec][gene]) for rec in gc_pd.columns if str(gc_pd[rec][gene])!='nan'])
    is_na = "-"*gene_length
    replace = [strain for strain in gc_pd.columns if str(gc_pd[strain][gene])=="nan"]
    
    for rec in replace:
        gc_filled[rec][gene] = is_na
    
    print(gene)
    #return gc_filled
        
def reorder():
    '''this function is meant to reorder the genes so that they ar
    e somewhat syntenic'''

def concatenate_ind(split):
    new_split = {}
    for strain_column in split.columns:
        full_genome = ''.join(split.columns[strain_column])
        new_split[strain_column] = full_genome
    return new_split

def concatenate_fasta_parallel(gc_filled):
    strains = gc_filled.columns
    split_strains = split_list(strains, 8)
    pool = mp.Pool(process = 8)
    results_concat = [pool.apply(concatenate_ind, args=split) for split in split_strains]
    
def gene_positions(gc_filled):
    gene_coordinates = []
    i=0
    for rec in gc.index:
       gene_length = len(gc[STRAIN][rec])
       for strain in gc.columns:
           gene_coordinates[rec] = [i, i+len(gc[strain][rec])]
    return temp_fasta, gene_coordinates   
    
def create_fasta(gc):
    gene_coordinates = []
    temp_fasta = {(strain, SeqRecord(Seq(""), id = strain)) for strain in gc.columns}
    i=0
    for rec in gc.index:
       gene_length = len(gc[STRAIN][rec])
       for strain in gc.columns:
           temp_fasta[strain].seq.append(gc[strain][rec])
           gene_coordinates[rec] = [i, i+len(gc[strain][rec])]
    return temp_fasta, gene_coordinates
        

# First must build the fasta file from the panX output. We need to build a data file that has the concatenated sequences from the pangenome and also records the positions of the g
panx_output = '/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/'

gense_cluster_path = panx_output + "geneCluster"

os.chdir("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/geneCluster")

sequenced_genome = os.listdir(panx_output+"input_GenBank")
strain_list = [rec.strip('.gbk') for rec in sequenced_genome] 
all_files = os.listdir(gene_cluster_path)
gene_list = [rec for rec in all_files if "fna" in rec]

gc_dict = {}
gc_pd=pd.DataFrame(columns=strain_list, index=[rec.strip('.fna') for rec in gene_list])

#now build a dictionary with the sequences for each gene in each strain

#num_cores = 8
pool = mp.Pool(processes = (16)) #mp.cpu_count() - 1))
gc_dict_args = [(fna, gc_dict) for fna in gene_list]
results = [pool.apply(process_cluster, args=rec) for rec in gc_dict_args]

i=0
for line in results:
    i=i+1
    print(i)
    for rec in list(line.keys()):
        gc_pd[rec[0]][rec[1]] = str(line[rec])


#results is a long list of dictionaries. Use chain map to make into one dictionary
#the following didn't work
#gc_full_dict = dict(ChainMap(*results))

#  


#gc_pd = pd.from_dict(gc_dict)
gc_pd.to_pickle("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/geneCluster.cpk")
gc_filled = fill_missing(gc_pd)
gc_filled.to_pickle("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/geneCluster_filled.cpk")
#gc_reordered = reorder(gc_filled)

whole_fasta = create_fasta(gc_reordered)

SeqIO.write(whole_fasta, "/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_file/fake_gene_SNP.fasta", "fasta")

#Now convert to vcf
os.chdir("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/")

os.system("/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/snp-sites -v -b fake_gene_SNP.fasta -o /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/my_1524.vcf " )

#Now convert to plink
os.system("/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/plink --vcf my_1524.vcf --double-id --maf 0.01 --recode --out my_1524.ped")
os.system("/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/vcftools
 --vcf my_1524.vcf --out my_1524 --plink-tped")


with open("snp_coordinates.txt", "w") as f:
    f.write(",".join(gene_coordinates))