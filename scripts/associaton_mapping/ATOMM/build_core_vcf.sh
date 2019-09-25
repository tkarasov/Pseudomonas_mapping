
code_direc=/ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/associaton_mapping/ATOMM/

panX_directory="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524"

output_direc="/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files"

cd $output_direc
/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/snp-sites -v -b /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/all_concat_1524.fasta -o all_concat_1524.vcf

/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/snp-sites -v /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/all_concat_1524.fasta -o all_concat_1524_snps_only.vcf

/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/snp-sites -v /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/all_concat_1524.fasta -c -o all_concat_1524_snps_only_no_ast.vcf

cp /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/all_concat.fasta $output_direc

/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/mapping/bin/snp-sites  -c -v $output_direc/all_concat.fasta -o all_concat_otu5_no_ast.vcf


vcftools --hap-r2 --vcf all_concat_otu5_no_ast.vcf --out all_concat_snps_rsquared
