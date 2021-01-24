#!/usr/bin/bash

# Edited by Filipe Moreira, Jaque Goes and Ingra Claro
# It's a dull version of the previouslly working pipeline. 
# make sure to run this on the directory containing the fastq_pass/ diretory 
# Usage: $ bash get_consensus_sequences.sh ~/path/to/ref/reference.fasta ~/path/to/bedfile		# ~/path/to/ref/reference.fasta reference genome with absolute path
													# ~/path/to/ref/bedfile primers bed file with absolute path
eval "$(conda shell.bash hook)"

conda activate nano-pipe										# This line and the previous one activate a required conda environment

ref=$(head -n 1 $1 | sed 's/>//g'); seq=$(tail -n 1 $1); n=$(echo $seq | wc -c); refpos=($ref:1-$n);	# get reference string to pass on nanopolish variants

echo $refpos

cat fastq_pass/*fastq > all_barcodes.fastq						# concatenate fastq files

porechop -i all_barcodes.fastq --format fastq -b demux -t 10 --require_two_barcodes 					# demux -> # alternativa: porechop -i all_barcodes.fastq --format fastq -b demux -t 10 --require_two_barcodes 

nanopolish index -d fast5_pass/ all_barcodes.fastq 							# index the reads -> required for 'nanopolish variants' -> raw files in fast5

mkdir safe
sed 's/fast5_/..\/fast5_/g' all_barcodes.fastq.index.readdb > all_barcodes.fastq.index.readdb2;	# These 4 lines edit the index map file to correct the path of fast5 files
mv all_barcodes.fastq.index.readdb safe
mv all_barcodes.fastq.index.readdb2 all_barcodes.fastq.index.readdb

cd demux_qcat												# enter he directory containing demuxed reads / barcode separated

for i in *fastq												# for each fastq file
do

	j=$(echo $i | sed 's/.fastq//g')								# salva o nome do barcode sem a extensao fastq

	echo "########## $j ##########"

	minimap2 -x map-ont -a $1 $i | samtools sort -o $j.sorted.bam					# map the reads against reference genome

	samtools index $j.sorted.bam									# indexing bam file
	
	align_trim.py --start --normalise 200 $2 --report $j.alignreport.txt < $j.sorted.bam 2> $j.alignreport.er | samtools view -bS - | samtools sort -T $j -o $j.trimmed.sorted.bam											# trimming bam file

	align_trim.py --normalise 200 $2 --report $j.alignreport.txt < $j.sorted.bam 2> $j.alignreport.er | samtools view -bS - | samtools sort -T $j -o $j.primertrimmed.sorted.bam										# Also trimming primers

	samtools index $j.trimmed.sorted.bam								# Index the trimmed bam files	

	samtools index $j.primertrimmed.sorted.bam							# again... 

	nanopolish variants --min-flanking-sequence 10 --fix-homopolymers -x 1000000 --progress -t 10 --reads ../all_barcodes.fastq -o $j.vcf -b $j.trimmed.sorted.bam -g $1 -w $refpos --snps --ploidy 1													# variant calling

	nanopolish variants --fix-homopolymers -x 1000000 --progress -t 10 --reads ../all_barcodes.fastq -o $j.primertrimmed.vcf -b $j.primertrimmed.sorted.bam -g $1 -w $refpos --snps --ploidy 1													# again...

	vcfextract $j > $j.variants.tab									# create a variants tab
	
	margin_cons.py $1 $j.vcf $j.primertrimmed.sorted.bam a > $j.consensus.fasta			# (finally) obtain consensus sequences

done

mkdir ../CONSENSUS_SEQUENCES										# make a directory

cp *consensus* ../CONSENSUS_SEQUENCES/									# copy consensus sequences for the later

cd ../													# get out of demux/





