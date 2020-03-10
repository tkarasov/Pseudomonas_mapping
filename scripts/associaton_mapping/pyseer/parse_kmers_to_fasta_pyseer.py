#!/usr/bin/py3

#  The goal of this script is to take the .assoc file after running a pyseer model,
# pull the associated kmers/unitigs and write a fasta file in which the fasta header has the variant frequencies and other info. The script also overlaps with bam
import os


class Kmer:
    def __init__(self, seq, af, pval, beta, i):
        self.seq = seq
        self.af = af
        self.pval = pval
        self.beta = beta
        self.name = "var" + str(i)


def rewrite_kmer(line, i):
    kmer = Kmer(line[0], line[1], line[3], line[4], i)
    return kmer


association_file = [line.strip().split() for line in open("./output/unitigs_6_pyseer.assoc").readlines()[1:]]

i = 0
kmer_dict = {}
myfile = open("./output/kmers.fasta", "a")
for line in association_file:
    i = i + 1
    kmer = rewrite_kmer(line, i)
    var1 = (">var" + str(i) + "; allele_freq=" + kmer.af + "; pval=" + kmer.pval + "; beta=" + kmer.beta)
    var2 = (kmer.seq)
    kmer_dict[kmer.name] = kmer.pval
    myfile.write(var1 + '\n' + var2 + "\n")

myfile.close()


##############
# map to genome with bwa
##############
os.system("bwa mem /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon.fasta ./output/kmers.fasta | \
	sed -e 's/utg000001c:1.0-5963307.0_pilon/utg000001l_p25.C2/g' | samtools sort  > ./output/kmers_bwa.bam")

##############
# Intersect with genome gff
##############
os.system("bedtools intersect -bed -b ./output/kmers_bwa.bam \
	-a /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon_no_fasta.gff \
	> ./output/intersect_gene.bed")

os.system("bedtools intersect -bed -a ./output/kmers_bwa.bam \
	-b /ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon_no_fasta.gff \
	> ./output/intersect_var.bed")

genes = [line.strip().split('\t') for line in open("./output/intersect_gene.bed").readlines()]

vars = [line.strip().split('\t') for line in open("./output/intersect_var.bed").readlines()]

for i in range(0, len(genes)):
    genes[i][0] = vars[i][3].strip(";")
    genes[i].append(vars[i][1])
    genes[i].append(kmer_dict[vars[i][3].strip(";")])

myfile = open("./output/kmers_genes_pvalues.txt", "w")
for line in genes:
    myfile.write(line[0] + "\t" + line[8] + "\t" + line[9] + line[10] + "\n")

myfile.close()
