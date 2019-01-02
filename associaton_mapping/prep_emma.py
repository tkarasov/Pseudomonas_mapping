#/usr/bin/python3
from Bio import SeqIO
import numpy
import sys
import os
import matplotlib 
from matplotlib import pyplot as plt
import pickle
import subprocess


def load_sorted_clusters(path): #Taken from Wei
    '''
    load gene clusters and sort 1st by abundance and then by clusterID
    '''
    geneClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    geneCluster_dt=load_pickle(geneClusterPath+'allclusters_postprocessed.cpk')
    from operator import itemgetter
    # sort by decreasing abundance (-v[0], minus to achieve decreasing)
    # followed by increasing strain count
    gc_items = geneCluster_dt.items()

    return sorted(geneCluster_dt.items(), key=lambda kv: (-itemgetter(0)(kv[1]),itemgetter(2)(kv[1])), reverse=False)
    #return sorted(geneCluster_dt.iteritems(),
    #            key=lambda (k,v): (-itemgetter(0)(v),itemgetter(2)(v)), reverse=False)


def load_pickle(filename):
    f = open(filename,"rb")
    p = pickle.load(f)
    f.close()
    return(p)

#http://zzz.bwh.harvard.edu/plink/binary.shtml
def build_map(sorted_genelist):
	my_map = []
	pos = 0
	i=1
	for rec in sorted_genelist:
		snp_id = rec[0]
		snp_pos = i
		pos = ["chr1", snp_id, 0, i]
		i = i+1
		my_map.append(pos)
	#for i in range(len(mat)):
	#	vals = [1, i, i]
	#	vals.extend(mat[i])
	#	mat_file.append(vals)
	return my_map

'''

def build_bim(my_bed, sorted_genelist):
	my_bim = []
	i = 0
	for line in my_bed:
		new =[1, sorted_genelist[i][0], '0', line[1], '1', '0']
		my_bim.append(new)
		i = i + 1
	return my_bim
'''

#https://training.h3abionet.org/postgraduate_workshop_2014/wp-content/uploads/2014/04/H3ABionet_2014_GWAS_2_Plink_Data_Format.pdf
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


#emmax-kin isn't working so well, but we could build a kinship matrix from the mash output
def build_kin(mash, strain_list):
	temp = [rec.strip(".fasta") for rec in mash.columns]
	mash.columns = temp
	mash.index = temp
	mash_subset = 1 - mash[strain_list].loc[strain_list]
	 mash_subset.to_csv("mash_kin.kinf", sep = " ", header = False, index = False)           



path = "/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/"

gene_pa = list(SeqIO.parse("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/geneCluster/genePresence.aln", 'fasta'))

mash = pandas.read_csv("my_addition.dist", sep = "\t", index_col = 0)

#gene_info = pickle.load(open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/geneID_to_description.cpk", "rb"))

sorted_genelist = load_sorted_clusters(path)
phenotypes_all = [line.strip().split() for line in open("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/image_analysis/mean_pixels.txt").readlines()]
phenotypes = [rec for rec in phenotypes_all if rec[0] in [thing.name for thing in gene_pa]]
strain_list = [line[0] for line in phenotypes]


#build bed_file precursors
my_map = build_map(sorted_genelist)
#my_bim = build_bim(my_bed, sorted_genelist)
my_ped,my_fam, my_pheno = build_ped(phenotypes, strain_list, gene_pa)

#build kinship matrix


handle = open("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/mapping/emma/strain.map", "w")
for line in my_map:
	handle.write( '\t'.join([str(rec) for rec in line]))
	handle.write("\n")

handle.close()

handle = open("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/mapping/emma/strain.ped", "w")
for line in my_ped:
	handle.write( '\t'.join([str(rec) for rec in line]))
	handle.write("\n")

handle.close()

handle = open("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/mapping/emma/strain.fam", "w")
for line in my_fam:
	handle.write( '\t'.join([str(rec) for rec in line]))
	handle.write("\n")

handle.close()

handle = open("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/mapping/emma/strain.pheno", "w")
for line in my_pheno:
	handle.write( '\t'.join([str(rec) for rec in line]))
	handle.write("\n")

handle.close()
#build bim 
'''Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
Variant identifier
Position in morgans or centimorgans (safe to use dummy value of '0')
Base-pair coordinate (normally 1-based, but 0 ok; limited to 231-2)
Allele 1 (corresponding to clear bits in .bed; usually minor)
Allele 2 (corresponding to set bits in .bed; usually major)
'''



#build ped_file


#build fam file Family ID ('FID')
'''Within-family ID ('IID'; cannot be '0')
Within-family ID of father ('0' if father isn't in dataset)
Within-family ID of mother ('0' if mother isn't in dataset)
Sex code ('1' = male, '2' = female, '0' = unknown)
Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
'''

#make bed file

#convert bed file to tped
cmd = ['/usr/bin/plink1', '--bfile', "strain", "--recode12", "--output-missing-genotype", "0", "--transpose", "--out", "strain"]
#plink --bfile [bed_prefix] (or --file [ped_prefix]) --recode12 --output-missing-genotype 0 --transpose --out [tped_prefix]
result = subprocess.run(cmd) #, stdout=subprocess.PIPE, input=input)

##plink1 --file strain --allow-no-sex --make-bed --maf 0.01 -out strain
cmd = ['plink1', '--file', 'strain', '--allow-no-sex', '--make-bed', '--maf', '0.01', '-out', 'strain']
result = subprocess.run()



#create kinship matrix: emmax-kin did not work well. Instead I used
#cmd = ['emmax-kin',  "strain.tped"]
#result = subprocess.run(cmd)




#run emmax
cmd = ['emmax', '-v', '-d', '10', '-t', 'strain.tped' '-p' 'strain.ped'  '-o' 'strain_out']
result = subprocess.run()



sed -i -e 's/-nan/0/g' strain.BN.kinf
emmax -v -d 10 -t strain -o strain.out -k mash_kin.kinf -p strain.pheno                                                                                                             
                         



GC00000053_48 -19853.69 4.262302e-25