#!/usr/bin/python3

'''this script is meant to take a list of gene numbers (based off of order in P/A aligned file and to return the information and sequence of that gene. The input is from panX output'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import numpy
import copy
import sys
import os
import matplotlib
from matplotlib import pyplot as plt
import pickle
import subprocess
import _pickle as cPickle


def load_sorted_clusters(path):  # Taken from Wei
    '''
    load gene clusters and sort 1st by abundance and then by clusterID
    '''
    geneClusterPath = '%s%s' % (path, 'protein_faa/diamond_matches/')
    geneCluster_dt = load_pickle(geneClusterPath + 'allclusters_postprocessed.cpk')
    from operator import itemgetter
    # sort by decreasing abundance (-v[0], minus to achieve decreasing)
    # followed by increasing strain count
    gc_items = geneCluster_dt.items()

    return sorted(geneCluster_dt.items(), key=lambda kv: (-itemgetter(0)(kv[1]), itemgetter(2)(kv[1])), reverse=False)
    # return sorted(geneCluster_dt.iteritems(),
    #            key=lambda (k,v): (-itemgetter(0)(v),itemgetter(2)(v)), reverse=False)


def load_pickle(filename):
    f = open(filename, "rb")
    p = cPickle.load(f)
    f.close()
    return(p)


# First must generate the sorted gene list and order the gene content
# my list of genes for which I want sequence and information
my_list = [4939, 5671]  # sys.argv[1]  # [5328, 5317, 5649]

path = "/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/"

fna = "/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/Ps_1524_all_fna.fna"
sorted_genelist = load_sorted_clusters(path)

record_dict = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))

# Dictionary of gene name (e.g. "p24.G9|DPNPJAGB_03821") and it's putative annotation and associated contig
gene_ids = load_pickle(path + "geneID_to_description.cpk")

all_recs = []
for rec in my_list: 
#[int(line.strip().split()[0]) for line in open(my_list).readlines()]:
    # choose the most relevant genoem
    GC = sorted_genelist[rec - 1][0]
    keep = sorted_genelist[rec - 1][1][1][0]
    look_up = gene_ids[keep]
    sequence = record_dict[keep]
    sequence.description = GC + ", " + str(rec) + ", " + look_up['annotation']
    all_recs.append(sequence)

SeqIO.write(all_recs, "/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/andy__significant_genes.fasta", "fasta")

# The output to return is a list of the genes, their putative function and a fasta with the sequence
