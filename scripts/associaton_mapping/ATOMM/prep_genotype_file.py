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
panx_output = sys.argv[1]#'/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/'
output_directory = sys.argv[2]
core=sys.argv[3]

gene_cluster_path = panx_output + "/geneCluster"

os.chdir(gene_cluster_path)

sequenced_genome = os.listdir(panx_output+"/input_GenBank")

strain_list = [rec.strip('.gbk') for rec in sequenced_genome] 

all_files = os.listdir(gene_cluster_path)

gene_list = [rec for rec in all_files if "_na_aln" in rec and "reduced" not in rec]

#if we are only considering the core genome then we only need core_gene_list
if core=="Core":
    #first get core names
    core_names=[line.strip().split() for line in os.system('ls ' panx_output)]

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
        gc_pd[rec[0].split("-")[0]][rec[1]] = str(line[rec])


#results is a long list of dictionaries. Use chain map to make into one dictionary
#the following didn't work
#gc_full_dict = dict(ChainMap(*results))

#gc_pd = pd.from_dict(gc_dict)
gc_pd.to_pickle("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/geneCluster.cpk")

#now split gc into data frames of 500 each
os.chdir(output_direc)
os.system("mkdir "+output_directory+"/temp2_gc_pd")
num_chunks = np.round(len(gc_pd.index)/500.0)
sub_pd = np.array_split(gc_pd, num_chunks)

i=0
for pd in sub_pd:
   i=i+1
   pd.to_pickle(output_directory+"/temp2_gc_pd"+"/temp_pd"+str(i)+".cpk")



