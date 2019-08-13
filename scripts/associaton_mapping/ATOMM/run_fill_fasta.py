#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:38:27 2019

@author: tkarasov
"""
import sys
sys.path.append("/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/lib/python3.7/site-packages/")
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
from multiprocess import Process, Manager
import multiprocess as mp
import os
import pickle
#from subprocess import Popen, PIPE
from joblib import Parallel, delayed
import pandas as pd
from collections import ChainMap
import numpy as np
import copy



def fill_missing(gc_filled):
    #gc_filled.loc[gene]
    i = 0
    for gene in gc_filled.index:
        i = i+1
        print(i)
        fill_missing_sub(gene)
        #print(gene)
    #processed_list = Parallel(n_jobs = 16)(delayed(fill_missing_sub)(gene) for gene in list(gc_filled.index))

def fill_missing_sub(gene):
#   for gene in gc_pd.index:
    gene_length = max([len(gc_pd[rec][gene]) for rec in gc_pd.columns if str(gc_pd[rec][gene])!='nan'])
    is_na = "-"*gene_length
    replace = [strain for strain in gc_pd.columns if str(gc_pd[strain][gene])=="nan"]
    
    for rec in replace:
        gc_filled[rec][gene] = is_na

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
    temp_fasta = {strain:SeqRecord.SeqRecord(seq="", id = strain) for strain in gc.columns}
    i=0
    for rec in gc.index:
       #gene_length = len(gc[STRAIN][rec])
       for strain in gc.columns:
           temp_fasta[strain].seq = Seq.Seq("".join([str(temp_fasta[strain].seq),gc[strain][rec]]))
           #gene_coordinates[rec] = [i, i+len(gc[strain][rec])]
    return temp_fasta
#, gene_coordinates
        

os.chdir("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/")
gc_pd  = pickle.load(open("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/geneCluster.cpk", "rb")) 
gc_filled = copy.deepcopy(gc_pd)
fill_missing(gc_filled)
gc_filled.to_pickle("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/geneCluster_filled1.cpk")
#gc_reordered = reorder(gc_filled)

whole_fasta = create_fasta(gc_filled)

SeqIO.write(list(whole_fasta.values()), "/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files/fake_gene_SNP2.fasta", "fasta")
