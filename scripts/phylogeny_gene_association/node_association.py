#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 10:30:30 2019

@author: tkarasov
"""
import os,sys,copy;import numpy as np
sys.path.append("/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/lib/python3.7/site-packages/")
from ete3 import Tree
import glob
import pickle
import pandas as pd
sys.path.append('./')
from treetime import *
import pickle as cPickle
#from treetime.treetime import treeanc as ta
import treetime.treeanc as ta
from treetime.gtr import GTR
#from treetime import io
from treetime import seq_utils
from Bio import Phylo, AlignIO

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


def infer_gene_gain_loss(rates = [1.0, 1.0], path_to_pangenome_dir = '/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/'):
    '''code nabbed and edited from panX'''
    
    # initialize GTR model with default parameters
    mu = np.sum(rates)
    gene_pi = np.array(rates)/mu
    gain_loss_model = GTR.custom(pi = gene_pi, mu=mu, W=np.ones((2,2)), alphabet = np.array(['0','1']))
    # add "unknown" state to profile
    gain_loss_model.profile_map['-'] = np.ones(2)
    #root_dir = os.path.dirname(os.path.realpath(__file__))

    # define file names for pseudo alignment of presence/absence patterns as in 001001010110
    #path_to_pangenome_dir='/ebio/ag-neher/share/users/wding/panX-refseq/data/Pseudomonadales'#sys.argv[1]
    nwk=path_to_pangenome_dir+"/vis/strain_tree.nwk"
    fasta=path_to_pangenome_dir+"/geneCluster/genePresence.aln"

    # instantiate treetime with custom GTR
    t = ta.TreeAnc(nwk, gtr =gain_loss_model, verbose=2)
    # fix leaves names since Bio.Phylo interprets numeric leaf names as confidence
    for leaf in t.tree.get_terminals():
        if leaf.name is None:
            leaf.name = str(leaf.confidence)
    t.aln = fasta
    t.tree.root.branch_length=0.0001
    t.reconstruct_anc(method='ml')
    for n in t.tree.find_clades():
        n.genepresence = n.sequence

    return t

def num_events(t, output_directory):
    mutation_dict={}
    for n in t.tree.find_clades():
        mut=n.mutations
        for mutation in mut:
            try:
                mutation_dict[mutation[1]]+=1
            except KeyError:
                mutation_dict[mutation[1]]=1
    handle=open(output_directory + "/pres_abs_events_gene.cpk", "wb")
    pickle.dump(mutation_dict, handle)
    handle.close()
            
def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]        

def calc_branch_length_in_tree(t, output_directory):
    '''takes output from running ancestral reconstruction in treetime, iterates through tree and if a gene present in node is present in parent, adds sum of branch length. Basically measure the amount of time (units in branch length) that gene has spent in the tree'''
    #do a breadth first search. Starting at the MRCA of tree, going to each child, store the branch length of parent to that child. Then go through every gene,
    # initialize
    parent=t.tree.root.sequence
    tracking_gene={}
    len_dict={}
    for i in range(0, len(t.tree.root.sequence)):
        tracking_gene[i]=t.tree.root.sequence[i]  
        len_dict[i]=0
    
    for n in t.tree.find_clades():
        try:
            p=get_parent(t.tree, n)
            pn=get_parent(t.tree, n).sequence
        except IndexError:
            continue
        for i in range(0, len(t.tree.root.sequence)):
            tracking_gene[i]=p.sequence[i] 
        branch=n.branch_length
        events=n.mutations
        #change in tracking_gene only those locations that have a putative mutation
        for mutation in events:
            print(tracking_gene[mutation[1]], mutation[2])
            tracking_gene[mutation[1]]=mutation[2]
        for gene in tracking_gene.keys():
            if tracking_gene[gene]=='1':
                len_dict[gene]=len_dict[gene]+branch

        #dump out branch lengths
        handle=open(output_directory + "/branch_gene.cpk", "wb")
        pickle.dump(len_dict, handle)
        handle.close()
        
def pull_gc(gene_number):
    #this function is to take the gene number output by panX and to return the original gene cluster

def pull_peg(value, genome, GC=True, gene_num=False):
    if GC==True:
        #if GC is true then must access GC file
    elif gene_num==True:
        #if peg then must 
    
panx_output = '/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/'# sys.argv[1] 

output_directory = '/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/presence_absence'#sys.argv[2] #


os.chdir(output_directory)
#import gene gain and losses from panX
path_to_pangenome_dir=panx_output
#t = path_to_pangenome_dir + '/vis/strain_tree.nwk'

t = infer_gene_gain_loss(path_to_pangenome_dir, rates = [1.0, 1.0])
calc_branch_length_in_tree(t, output_directory)#(path_to_pangenome_dir)
num_events(t, output_directory)

#get ancestor node 
intree = path_to_pangenome_dir+"/vis/strain_tree.nwk"

mytree = Tree(intree, format=1)

otu5_node = t.tree.common_ancestor("p4.C6", "p20.F8")#("p20.G1", "p26.F9")

all_mutations = t.get_mutations(otu5_node) #596

gains_otu5 = [mutation for mutation in all_mutations if mutation[0]=="0"] #586

gains_other = [mutation for mutation in all_mutations if mutation[0]=="1"] #10 genes?

fna = "/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/Ps_1524_all_fna.fna"

sorted_genelist = load_sorted_clusters(path_to_pangenome_dir)

record_dict = SeqIO.to_dict(SeqIO.parse(fna, "fasta")) #this takes a long time to load

# Dictionary of gene name (e.g. "p24.G9|DPNPJAGB_03821") and it's putative annotation and associated contig
gene_ids = load_pickle(path_to_pangenome_dir + "geneID_to_description.cpk")

all_recs = []
for rec in gains_otu5: 
#[int(line.strip().split()[0]) for line in open(my_list).readlines()]:
    # choose the most relevant genoem
    GC = sorted_genelist[rec - 1][0]
    keep = sorted_genelist[rec - 1][1][1][0]
    look_up = gene_ids[keep]
    sequence = record_dict[keep]
    sequence.description = GC + ", " + str(rec) + ", " + look_up['annotation']
    all_recs.append(sequence)

SeqIO.write(all_recs, "/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/otu5_significant_genes.fasta", "fasta")

