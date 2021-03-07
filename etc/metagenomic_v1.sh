#!/bin/bash

THREADS="12" #12
NUMCALLER="46" #NVIDIAGeForceRTX2060=60; NVIDIAGeForceRTX2080=46
GPUPERDEVICE="64" #NVIDIAGeForceRTX2060=32; NVIDIAGeForceRTX2080=64

RAWDIR="/media/rklab/Data/DENV/DENV_FTA_1/"
HACDEMUXDIR="/media/rklab/Data/DENV/WO_require_barcodes_both_ends/DENV_FTA_1_HAC_Demux/"

BARCODEKIT="EXP-NBD104" #EXP-NBD104
FLOWCELL="FLO-MIN106" #FLO-MIN106
SEQKIT="SQK-LSK109" #SQK-LSK109

LIBNAME="DENV_FTA_1"

SAMPLEID1="barcode01"
SAMPLEID2="barcode02"
SAMPLEID3="barcode03"
SAMPLEID4="barcode04"
SAMPLEID5="barcode05"
SAMPLEID6="barcode06"
SAMPLEID7="barcode07"
SAMPLEID8="barcode08"
SAMPLEID9="barcode09"
SAMPLEID10="barcode10"
SAMPLEID11="barcode11"
SAMPLEID12="barcode12"

REFSEQ="GCF_000001405.39_GRCh38.p13_genomic.fna.gz"

guppy_basecaller -r -i ${RAWDIR} -s ${HACDEMUXDIR} --flowcell ${FLOWCELL} --kit ${SEQKIT} --barcode_kits ${BARCODEKIT} --require_barcodes_both_ends --trim_barcodes --qscore_filtering --num_barcode_threads ${THREADS} -x auto --gpu_runners_per_device ${NUMCALLER} --num_callers ${NUMCALLER} --verbose_logs
source activate pycoqc
pycoQC -q -f ${HACDEMUXDIR}/sequencing_summary.txt -o ${LIBNAME}_qc.html --report_title ${LIBNAME}
[ ! -d "01_guppy" ] && mkdir -p 01_guppy
[ -d "${HACDEMUXDIR}/pass/barcode01" ] && cat ${HACDEMUXDIR}/pass/barcode01/*.fastq > 01_guppy/${SAMPLEID1}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode02" ] && cat ${HACDEMUXDIR}/pass/barcode02/*.fastq > 01_guppy/${SAMPLEID2}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode03" ] && cat ${HACDEMUXDIR}/pass/barcode03/*.fastq > 01_guppy/${SAMPLEID3}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode04" ] && cat ${HACDEMUXDIR}/pass/barcode04/*.fastq > 01_guppy/${SAMPLEID4}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode05" ] && cat ${HACDEMUXDIR}/pass/barcode05/*.fastq > 01_guppy/${SAMPLEID5}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode06" ] && cat ${HACDEMUXDIR}/pass/barcode06/*.fastq > 01_guppy/${SAMPLEID6}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode07" ] && cat ${HACDEMUXDIR}/pass/barcode07/*.fastq > 01_guppy/${SAMPLEID7}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode08" ] && cat ${HACDEMUXDIR}/pass/barcode08/*.fastq > 01_guppy/${SAMPLEID8}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode09" ] && cat ${HACDEMUXDIR}/pass/barcode09/*.fastq > 01_guppy/${SAMPLEID9}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode10" ] && cat ${HACDEMUXDIR}/pass/barcode10/*.fastq > 01_guppy/${SAMPLEID10}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode11" ] && cat ${HACDEMUXDIR}/pass/barcode11/*.fastq > 01_guppy/${SAMPLEID11}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode12" ] && cat ${HACDEMUXDIR}/pass/barcode12/*.fastq > 01_guppy/${SAMPLEID12}.fastq
[ ! -d "02_cutadapt" ] && mkdir -p 02_cutadapt
[ ! -d "03_nanofilt" ] && mkdir -p 03_nanofilt
[ ! -d "04_prinseq" ] && mkdir -p 04_prinseq
[ ! -d "05_samtools-sorted" ] && mkdir -p 05_samtools-sorted
[ ! -d "06_samtools-unmapped" ] && mkdir -p 06_samtools-unmapped
[ ! -d "07_racon" ] && mkdir -p 07_racon
source activate metagenomic
for i in $(find ./01_guppy/ -type f -name "*.fastq" | while read o; do basename $o | rev | cut -c 7- | rev; done | sort | uniq); do
	cutadapt -j ${THREADS} -g "GTTTCCCACTGGAGGATA" -e 0.2 --trimmed-only -o 02_cutadapt/${i}.cutadapt.fastq 01_guppy/${i}.fastq;
	NanoFilt -l 200 --headcrop 10 < 02_cutadapt/${i}.cutadapt.fastq > 03_nanofilt/${i}.nanofilt.fastq;
	prinseq++ -threads ${THREADS} -lc_dust=0.1 -rm_header -fastq 03_nanofilt/${i}.nanofilt.fastq -out_good 04_prinseq/${i}.good.fastq -out_bad 04_prinseq/${i}.bad.fastq;
	minimap2 -ax map-ont -t ${THREADS} ${REFSEQ} 04_prinseq/${i}.good.fastq | samtools sort -@ 12 -o 05_samtools-sorted/${i}.sorted.bam -;	
	samtools view -bS -f 4 05_samtools-sorted/${i}.sorted.bam > 06_samtools-unmapped/${i}.unmapped.sam -@ 12;
	samtools fastq -f 4 06_samtools-unmapped/${i}.unmapped.sam > 06_samtools-unmapped/${i}.unmapped.fastq -@ 12;
	minimap2 -ax ava-ont -t ${THREADS} 06_samtools-unmapped/${i}.unmapped.fastq 06_samtools-unmapped/${i}.unmapped.fastq > 07_racon/${i}.overlap.sam;
	racon -t ${THREADS} -f -u 06_samtools-unmapped/${i}.unmapped.fastq 07_racon/${i}.overlap.sam 06_samtools-unmapped/${i}.unmapped.fastq > 07_racon/${i}.racon.fasta;
done