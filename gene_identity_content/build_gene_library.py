#!/usr/bin/python3

'''this script converts binary fasta to fake nt fasta'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import numpy
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
path = "/Users/tkarasov/work_main/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/"
sorted_genelist = load_sorted_clusters(path)
gene_ids = load_pickle(path + "geneID_to_description.cpk")
my_list = [5328, 5317, 5649]
tot = {}
for rec in my_list:
    keep = sorted_genelist[rec - 1][1][1][0]
    look_up = gene_ids[keep]
    tot[rec] = [keep, look_up]
