
code_direc=/ebio/abt6_projects8/Pseudomonas_mapping/code_Pseudomonas_mapping_git/scripts/associaton_mapping/ATOMM/

panX_directory="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524"

output_direc="/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/SNP_files"

python $code_direc/prep_genotype.py $panX_directory

qsub -t 1-$num_chunks $code_direc/run_array_build_fasta.sh $output_direc/temp_gc_pd

qsub -hold_jid "build_fasta" $code_direc/run_prep_for_gemma_atomm.sh $output_direc
