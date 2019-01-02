

#script modified from Hawk (https://github.com/atifrahman/HAWK/blob/master/ecoli_analysis/countKmers) for counting kmers in genomes
CORES=30 #number of cores to use for blast searches
KMERSIZE=31 # RD:61

#modified from NIKS script

dir=/home/atif/1000_genomes/Ecoli		#directory for read files 
hawkDir=/home/atif/hawk-0.8.3-beta-linux			#directory where hawk is installed
jellyfishDir=/home/atif/jellyfish-Hawk/bin		#directory where jellyfish is installed
sortDir=/home/atif/coreutils/deps/bin		#directory where parallel sort is installed

cd ${dir}

for file in `cat links.txt*`
do
	OUTPREFIX=$file
	
	mkdir ${file}
	
	cd ${file}

	fastq-dump /home/atif/ncbi/sra/${file}.sra

	mkdir ${OUTPREFIX}_kmers

	${jellyfishDir}/jellyfish count -C -o ${OUTPREFIX}_kmers/tmp -m ${KMERSIZE} -t ${CORES} -s 20G *.fastq 	#change if gzipped

	COUNT=$(ls ${OUTPREFIX}_kmers/tmp* |wc -l)

	if [ $COUNT -eq 1 ]
	then
 		mv ${OUTPREFIX}_kmers/tmp_0 ${OUTPREFIX}_kmers_jellyfish
	else
		${jellyfishDir}/jellyfish merge -o ${OUTPREFIX}_kmers_jellyfish ${OUTPREFIX}_kmers/tmp*
	fi
	rm -rf ${OUTPREFIX}_kmers
	
	COUNT=$(ls ${OUTPREFIX}_kmers_jellyfish |wc -l)

	if [ $COUNT -eq 1 ]
	then

		${jellyfishDir}/jellyfish histo -f -o ${OUTPREFIX}.kmers.hist.csv -t ${CORES} ${OUTPREFIX}_kmers_jellyfish
		awk '{print $2"\t"$1}' ${OUTPREFIX}.kmers.hist.csv > ${OUTPREFIX}_tmp
		mv ${OUTPREFIX}_tmp ${OUTPREFIX}.kmers.hist.csv

		awk -f ${hawkDir}/countTotalKmer.awk ${OUTPREFIX}.kmers.hist.csv >> ${dir}/total_kmer_counts.txt

		CUTOFF=1 
		echo $CUTOFF > ${OUTPREFIX}_cutoff.csv


		${jellyfishDir}/jellyfish dump -c -L `expr $CUTOFF + 1` ${OUTPREFIX}_kmers_jellyfish > ${OUTPREFIX}_kmers.txt 
		sort --parallel=${CORES} -n -k 1 ${OUTPREFIX}_kmers.txt > ${OUTPREFIX}_kmers_sorted.txt
	
		rm ${OUTPREFIX}_kmers_jellyfish	
		rm ${OUTPREFIX}_kmers.txt		
			
		echo "${dir}/${OUTPREFIX}/${OUTPREFIX}_kmers_sorted.txt" >> ${dir}/sorted_files.txt
		
	fi

	rm *.fastq

	cd ..

done