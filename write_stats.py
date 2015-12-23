import sys
import csv
import os

write_stats_file = str(sys.argv[1])
fastq_file = str(sys.argv[2])
reads_processed = int(sys.argv[3].replace(',',''))
reads_with_adapters = int(sys.argv[4].replace(',',''))
reads_too_short = int(sys.argv[5].replace(',',''))
reads_written = int(sys.argv[6].replace(',',''))
bowtie_reads = int(sys.argv[7].replace(',',''))
aligned_0_times = int(sys.argv[8].replace(',',''))
aligned_1_time = int(sys.argv[9].replace(',',''))
aligned_g1_times = int(sys.argv[10].replace(',',''))
miRBase_aln_rate = float(sys.argv[11])
mirUtils_num_alignments = int(sys.argv[12].replace(',',''))
hairpin_mirs = int(sys.argv[13].replace(',',''))
group_mirs = int(sys.argv[14].replace(',',''))
family_mirs = int(sys.argv[15].replace(',',''))
cluster_mirs = int(sys.argv[16].replace(',',''))

percent_reads_too_short = float(reads_too_short * 100) / float(reads_processed)
reads_without_adapters = reads_processed - reads_with_adapters
percent_reads_without_adapters = float(reads_without_adapters * 100) / float(reads_processed)
miRBase_total_aln_rate = float(miRBase_aln_rate * bowtie_reads) / float(reads_processed)

with open(write_stats_file, 'a') as f:
    writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    if os.stat(write_stats_file).st_size == 0:
    	writer.writerow(["filename", "reads_processed",
    		"percent_reads_without_adapters", "percent_reads_too_short", "reads_written",
    		"bowtie_reads", "miRBase_aln_rate", "miRBase_total_aln_rate", "hairpin_mirs", "group_mirs"])
    writer.writerow([fastq_file, reads_processed,
    	percent_reads_without_adapters, percent_reads_too_short, reads_written,
    	bowtie_reads, miRBase_aln_rate, miRBase_total_aln_rate, hairpin_mirs, group_mirs])