#!/usr/bin/bash

# get_consensus_sequences.sh #

# Usage: $ bash get_consensus_sequences.sh ~/path/to/library ~/path/to/ref/reference.fasta ~/path/to/bed/ref.bed	number-of-threads coverage-for-normalization
# Example: bash ~/Desktop/SCRIPTS/get_consensus_sequences.sh ~/Desktop/LIBRARIES/SARS-CoV-2_LNCC_Library1_20200323/2_LNCC_Library1_20200323/20200322_1246_MN25989_FAL86072_b61e7753/  ~/Desktop/PRIMER_SCHEMES/nCoV-2019.reference.fasta ~/Desktop/PRIMER_SCHEMES/nCoV-2019.scheme.bed 8 20

# Edited by Filipe Moreira, Jaque Goes and Ingra Claro
# It's a simple pipeline to obtain consensus sequences from a MiniION sequencing run

# Edit v1.1 
# This version requires five arguments: 
# 1 - Path for the main library directory (containing the fastq_pass/ directory created during MiniION sequencing);
# 2 - Path for the reference genome (.fasta file)
# 3 - Path for the .bed file for each reference, specifying primer positions
# 4 - Maximum number of threads to be used in the pipeline
# 5 - Coverage used for normalization with 'align_trim.py'

# Warning: Do not try to run this pipeline without coducting some changes on the set up highlighted on the README.txt file
# These changes are the following:
# 1 - Create the necessary conda env
# 2 - Install qcat
# 2 - Get the latest version of nanopolish
# 3 - Change the way python scripts import dependencies (README.txt)
# 4 - Runtime may be optimized if nanopolish index is able to use a sequencing summary file (README.txt)


# The next two lines activate a required conda environment
eval "$(conda shell.bash hook)"
conda activate nano-pipe

# Go to the main library directory 
cd $1

# Get reference string to pass on nanopolish variants
ref=$(head -n 1 $2 | sed 's/>//g'); seq=$(tail -n 1 $2); n=$(echo $seq | wc -c); refpos=($ref:1-$n);

# Concatenate fastq files
echo "Concatenating fastq files"
cat fastq_pass/*fastq > all_barcodes.fastq

# Edit the sequencing summary file speeding up reads indexing
#cp *summary* summary.txt
#R --no-save < ~/Desktop/SCRIPTS/reformat_summary.R 

# Index the reads -> required for 'nanopolish variants' -> raw files in fast5
echo "Indexing reads"
nanopolish index -d /home/artic/Desktop/LIBRARIES/nCov19_SP_Library7_20200331/nCov19_SP_Library7_20200331/fast5_pass/ all_barcodes.fastq --verbose
#nanopolish index -d fast5_pass/ all_barcodes.fastq --verbose --sequencing-summary formated_summary_file.txt


# Correct the path for fast5 files
sed -i -E 's/fast5_pass\//..\/fast5_pass\//' all_barcodes.fastq.index.readdb

# demux
echo "Demuxing"
porechop --verbosity 1 -t 2 -i all_barcodes.fastq -b demux --barcode_threshold 70 --discard_middle --discard_unassigned --barcode_diff 5 --require_two_barcodes	  	  

# Enter demuxed reads directory
cd demux

# Iterating over all fastq files (all separate barcodes)
echo "Iterating over fastq files to obtain consensus sequences"
for i in *fastq
do

	# Gets barcode name without .fastq extension
	j=$(echo $i | sed 's/.fastq//g')
	
	# Print barcode on the screen
	echo "########## $j ##########"

	# Map reads against reference genome
	minimap2 -x map-ont -a $2 $i | samtools sort -o $j.sorted.bam

	# Index bam file
	samtools index $j.sorted.bam
	
	# trim bam file
	align_trim.py --start --normalise $5 $3 --report $j.alignreport.txt < $j.sorted.bam 2> $j.alignreport.er | samtools view -bS - | samtools sort -T $j -o $j.trimmed.sorted.bam
	
	# Also trimming primers
	align_trim.py --normalise $5 $3 --report $j.alignreport.txt < $j.sorted.bam 2> $j.alignreport.er | samtools view -bS - | samtools sort -T $j -o $j.primertrimmed.sorted.bam

	# Index the trimmed bam files
	samtools index $j.trimmed.sorted.bam
	
	# again...
	samtools index $j.primertrimmed.sorted.bam 

	# Variant calling
	nanopolish variants --min-flanking-sequence 10 --fix-homopolymers -x 1000000 --progress -t $4 --reads ../all_barcodes.fastq -o $j.vcf -b $j.trimmed.sorted.bam -g $2 -w $refpos --snps --ploidy 1

	# again...
	nanopolish variants --fix-homopolymers -x 1000000 --progress -t $4 --reads ../all_barcodes.fastq -o $j.primertrimmed.vcf -b $j.primertrimmed.sorted.bam -g $2 -w $refpos --snps --ploidy 1

	# Create variants tab
	vcfextract $j > $j.variants.tab
	
	# (Finally) obtain consensus sequences
	margin_cons.py $2 $j.vcf $j.primertrimmed.sorted.bam a > $j.consensus.fasta

done

# Make a directory containing all consensus sequence for the library
mkdir ../CONSENSUS_SEQUENCES/
cat *consensus* > consensus_sequences.fasta
cp consensus_sequences.fasta ../CONSENSUS_SEQUENCES/
cd ../
