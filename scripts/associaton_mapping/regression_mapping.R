library(coxme)
library(snpStats)
library(rrBLUP)
library(rjson)
setwd("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/mapping/emma/")
#pheno = read.table("strain.pheno")
kinship = read.table("mash_kin.kinf")
rownames( kinship ) = phenotype$member
#phenotype = snp_data$fam[,c("member","affected")]
#genotype = snp_data$genotypes[,c("dimnames[[2]]")]

#snp_data = read.plink(bed = "strain.bed")
map =read.table("strain.map")
genotpe = read.table("strain.ped")
phenotype = genotpe[,c(1,6)]
colnames(phenotype) = c("member", "affected")
snp_info = t(genotpe[ , c(7:dim(genotpe)[2]) ][ , c(TRUE,FALSE) ])

#reassign : https://cran.r-project.org/web/packages/rrBLUP/rrBLUP.pdf
snp_info[snp_info == 1] <- -1
snp_info[snp_info == 2] <- 1
snp_info = as.data.frame(snp_info)
colnames( snp_info ) = phenotype$member
genot_final = cbind( map[ , c(2,1,4) ], snp_info )
pdf("rrBLUP_output_kinship.pdf")
scores = GWAS(phenotype, genot_final, min.MAF = 0.05, plot = TRUE)
dev.off()

#qqplot
pdf("rrBLUP_qqplot.pdf")
exp.pvalues<-(rank((10^(-scores$affected)), ties.method="first")+.5)/(length(scores$affected)+1)
my.pvalues <- scores$affected
plot(-log10(exp.pvalues), (my.pvalues), asp=1, pch = 20)
abline(0,1,color = "RED")
dev.off()


#look at genes
cluster_info <- fromJSON( file = "/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/vis/geneCluster.json", simplify = TRUE, collapse =' ')


