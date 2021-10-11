library(dplyr)
#this script reads in pyseer and tnseq results and looks at the overlap


gwa = read.table("/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/pyseer/output/kmers_genes_pvalues.txt", sep = "\t"
                 )

colnames(gwa) = c("var", "gene", "pos_start", "pval_gwa")

#now must calculate minimum per gene
genes <- as.character(unique(gwa$gene))
gwa$min_pval <-NA
for(gene in genes){
  is_gene <- gwa[which(gwa$gene==gene),]$pval_gwa
  min_pva <- min(is_gene, na.rm = T)
  gwa[which(gwa$gene==gene),]$min_pval <- min_pva
  
}



tnseq = read.table("/ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process/edgeR-wt-vs-mut.csv", sep=",",
                   header=T, row.names = 1)
tnseq$gene = rownames(tnseq)

merged = full_join(tnseq, gwa)
plot(log10(merged$PValue), log10(merged$pval_gwa))

plot(merged$pos_start, -log10(merged$FDR))

sig_ass <- merged[merged$FDR<0.0000000001,] #& merged$pval_gwa<0.0000000000001,]

sim_unif <-runif(0,1, n= length(gwa$pval_gwa))

pdf("/ebio/abt6_projects8/Pseudomonas_mapping/data/fig_misc/qqplot_pyseer.pdf", useDingbats = FALSE, font = "ArialMT")
qqplot(y = -log10(gwa$pval_gwa), x= -log10(sim_unif))
abline(0,1, col = "RED")
dev.off()