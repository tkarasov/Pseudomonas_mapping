cd /ebio/abt6_projects8/Pseudomonas_mapping/Programs/flo/flo_p25.c2
#original assembly
/ebio/abt6_projects9/Pseudomonas_diversity/data/submitted_genomes/p25.C2.fasta.gz

#original annotation
/ebio/abt6_projects9/Pseudomonas_diversity/data/deprecated_pseudomonas_genome_assembly/assembly_annotation/plate25/C2/annotation/

#new assembly
/ebio/abt6_projects9/PseudomonasEvolution/data/ONT_genomes/assembly/qcat/pomoxis/p25.C2/4Xracon_meduka/medaka/pilon/p25.C2.contigs.second_polished.pilon.fasta


# https://github.com/wurmlab/flo
# liftover of genome annotation using the git

/ebio/abt6_projects8/Pseudomonas_mapping/Programs/flo/gff_remove_feats.rb /ebio/abt6_projects9/Pseudomonas_diversity/data/deprecated_pseudomonas_genome_assembly/assembly_annotation/plate25/C2/annotation/plate25.C2.annotation.gff > plate25.C2.annotation_longest_transcripts.gff

rake -f /ebio/abt6_projects8/Pseudomonas_mapping/Plsrograms/flo/Rakefile

#rake seems to fail after building the chain file. but perhaps I can use the chain file and then do the annotation lifting from there
mkdir run/plate25.C2.annotation
liftOver -gff /ebio/abt6_projects9/Pseudomonas_diversity/data/deprecated_pseudomonas_genome_assembly/assembly_annotation/plate25/C2/annotation/plate25.C2.annotation.gff run/liftover.chn run/plate25.C2.annotation/lifted.gff3 run/plate25.C2.annotation/unlifted.gff3
Reading liftover chains
Mapping coordinates
WARNING: -gff is not recommended.
Use 'ldHgGene -out=<file.gp>' and then 'liftOver -genePred <file.gp>'
Expecting at least 8 words line 10653 of /ebio/abt6_projects9/Pseudomonas_diversity/data/deprecated_pseudomonas_genome_assembly/assembly_annotation/plate25/C2/annotation/plate25.C2.annotation.gff
rake aborted!
Command failed with status (255): [liftOver -gff /ebio/abt6_projects9/Pseudom...]
/usr/lib/ruby/vendor_ruby/rake/file_utils.rb:67:in `block in create_shell_runner'
/usr/lib/ruby/vendor_ruby/rake/file_utils.rb:57:in `sh'
/ebio/abt6_projects8/Pseudomonas_mapping/Programs/flo/Rakefile:45:in `block (2 levels) in <top (required)>'
/ebio/abt6_projects8/Pseudomonas_mapping/Programs/flo/Rakefile:40:in `each'
/ebio/abt6_projects8/Pseudomonas_mapping/Programs/flo/Rakefile:40:in `block in <top (required)>'
/usr/lib/ruby/vendor_ruby/rake/task.rb:271:in `block in execute'
/usr/lib/ruby/vendor_ruby/rake/task.rb:271:in `each'
/usr/lib/ruby/vendor_ruby/rake/task.rb:271:in `execute'
/usr/lib/ruby/vendor_ruby/rake/task.rb:213:in `block in invoke_with_call_chain'
/usr/lib/ruby/2.5.0/monitor.rb:226:in `mon_synchronize'
/usr/lib/ruby/vendor_ruby/rake/task.rb:193:in `invoke_with_call_chain'
/usr/lib/ruby/vendor_ruby/rake/task.rb:182:in `invoke'
/usr/lib/ruby/vendor_ruby/rake/application.rb:160:in `invoke_task'
/usr/lib/ruby/vendor_ruby/rake/application.rb:116:in `block (2 levels) in top_level'
/usr/lib/ruby/vendor_ruby/rake/application.rb:116:in `each'
/usr/lib/ruby/vendor_ruby/rake/application.rb:116:in `block
