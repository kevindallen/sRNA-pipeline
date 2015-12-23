untrimmed=(
tanriverdi_seq/IonXpressRNA_013
tanriverdi_seq/IonXpressRNA_014
tanriverdi_seq/IonXpressRNA_015
tanriverdi_seq/IonXpressRNA_016
);

#AGATCGGAAGAGCACACGTCT
#TGGAATTCTCGGGTGCCAAG

write_stats_file="tanriverdi_seq/run_stats.csv"
ext=".fastq"
counter=1

printf "\n****************************************************************\n"
printf "************ Small RNA Sequencing Analysis Pipeline ************\n"
printf "****************************************************************\n\n"

for fq in ${untrimmed[*]}
do
	printf "  ($counter of ${#untrimmed[*]}) Trimming adapters (25%%) . . . . .          \r"
	#cutadapt -a AGATCGGAAGAGCACACGTCT -o $fq.trim1.fq --untrimmed-output $fq.untrim.fq --minimum-length 24 --too-short-output $fq.short.fq $fq$ext &> $fq.cutadapt.out
	#cutadapt -u 4 -o $fq.trim2.fq $fq.trim1.fq >/dev/null
	cutadapt -u 4 -o $fq.trim2.fq $fq$ext >/dev/null
	cutadapt -u -4 -o $fq.trim3.fq $fq.trim2.fq >/dev/null
	printf "  ($counter of ${#untrimmed[*]}) Aligning $item via bowtie2 (50%%) . . . . . \r"
	#bowtie -S $index $fq.trim3.fq $fq.sam
	bowtie2 -x /usr/local/bowtie2-2.2.6/human_hairpin -U $fq.trim3.fq -S $fq.sam &> $fq.bowtie2.out
	#echo "${bowtie2_out}"
	printf "  ($counter of ${#untrimmed[*]}) Running samtools (75%%) . . . . .           \r"
	samtools view -bS -o $fq.bam $fq.sam
	samtools sort $fq.bam $fq.sorted
	samtools index $fq.sorted.bam
	samtools idxstats $fq.sorted.bam > $fq.counts
	printf "  ($counter of ${#untrimmed[*]}) Running mirUtils (100%%) . . . . .          \r"
	mirUtils mbaseMirStats --organism=hsa --out-prefix=$fq $fq.bam &> $fq.mirUtils.out
	
	reads_processed="$(cat $fq.cutadapt.out | grep -oP 'Total reads processed:\s*\K[,\d]*')"
	reads_with_adapters="$(cat $fq.cutadapt.out | grep -oP 'Reads with adapters:\s*\K[,\d]*')"
	reads_too_short="$(cat $fq.cutadapt.out | grep -oP 'Reads that were too short:\s*\K[,\d]*')"
	reads_written="$(cat $fq.cutadapt.out | grep -oP 'Reads written \(passing filters\):\s*\K[,\d]*')"
	
	bowtie_reads="$(cat $fq.bowtie2.out | grep -oP '\d*(?=\sreads; of these:)')"
	aligned_0_times="$(cat $fq.bowtie2.out | grep -oP '\d*(?=\s\([\d\.%]*\)\saligned 0 times)')"
	aligned_1_time="$(cat $fq.bowtie2.out | grep -oP '\d*(?=\s\([\d\.%]*\)\saligned exactly 1 time)')"
	aligned_g1_times="$(cat $fq.bowtie2.out | grep -oP '\d*(?=\s\([\d\.%]*\)\saligned >1 times)')"
	bowtie_aln_rate="$(cat $fq.bowtie2.out | grep -oP '[\d\.]*(?=%\soverall alignment rate)')"
	
	mirUtils_num_alignments="$(cat $fq.mirUtils.out | grep -oP '\d*(?=\salignments)')"
	hairpin_mirs="$(cat $fq.mirUtils.out | grep -oP '\d*(?= hairpin mirs to )')"
	group_mirs="$(cat $fq.mirUtils.out | grep -oP '\d*(?= group mirs to )')"
	family_mirs="$(cat $fq.mirUtils.out | grep -oP '\d*(?= family mirs to )')"
	cluster_mirs="$(cat $fq.mirUtils.out | grep -oP '\d*(?= cluster mirs to )')"
	
	python write_stats.py $write_stats_file $fq \
	$reads_processed $reads_with_adapters $reads_too_short $reads_written \
	$bowtie_reads $aligned_0_times $aligned_1_time $aligned_g1_times $bowtie_aln_rate \
	$mirUtils_num_alignments $hairpin_mirs $group_mirs $family_mirs $cluster_mirs
	
	printf "  FastQ file $counter of ${#untrimmed[*]} processed successfully!              \n"
	counter=$(($counter+1))
done
printf "\n"