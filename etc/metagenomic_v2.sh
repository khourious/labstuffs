#!/bin/bash

CONFIG="dna_r9.4.1_450bps_hac.cfg" #dna_r9.4.1_450bps_hac.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" #barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg
MINSCORE="8"
TRIMADAPTER="18"

THREADS="12"
NUMCALLER="46" #NVIDIAGeForceRTX2060=60; NVIDIAGeForceRTX2080=46
GPUPERDEVICE="64" #NVIDIAGeForceRTX2060=32; NVIDIAGeForceRTX2080=64

RAWDIR="/media/rklab/Data/DENV/DENV_FTA_1/"
HACDEMUXDIR="/media/rklab/Data/DENV/WO_require_barcodes_both_ends/DENV_FTA_1_"
RUNNAME="DENV_FTA_1"

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

guppy_basecaller -r -i ${RAWDIR} -s ${HACDEMUXDIR} -c ${CONFIG} --arrangements_files ${ARRANGEMENTS} --qscore_filtering --min_qscore ${MINSCORE} --require_barcodes_both_ends --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER} --num_barcode_threads ${THREADS} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --num_callers ${NUMCALLERS} --verbose_logs

source activate pycoqc
pycoQC -q -f ${HACDEMUXDIR}/sequencing_summary.txt -o ${RUNNAME}_qc.html --report_title ${RUNNAME}
[ ! -d "hac_demux" ] && mkdir -p hac_demux
[ -d "${HACDEMUXDIR}/pass/barcode01" ] && cat ${HACDEMUXDIR}/pass/barcode01/*.fastq > ${SAMPLEID1}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode02" ] && cat ${HACDEMUXDIR}/pass/barcode02/*.fastq > ${SAMPLEID2}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode03" ] && cat ${HACDEMUXDIR}/pass/barcode03/*.fastq > ${SAMPLEID3}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode04" ] && cat ${HACDEMUXDIR}/pass/barcode04/*.fastq > ${SAMPLEID4}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode05" ] && cat ${HACDEMUXDIR}/pass/barcode05/*.fastq > ${SAMPLEID5}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode06" ] && cat ${HACDEMUXDIR}/pass/barcode06/*.fastq > ${SAMPLEID6}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode07" ] && cat ${HACDEMUXDIR}/pass/barcode07/*.fastq > ${SAMPLEID7}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode08" ] && cat ${HACDEMUXDIR}/pass/barcode08/*.fastq > ${SAMPLEID8}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode09" ] && cat ${HACDEMUXDIR}/pass/barcode09/*.fastq > ${SAMPLEID9}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode10" ] && cat ${HACDEMUXDIR}/pass/barcode10/*.fastq > ${SAMPLEID10}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode11" ] && cat ${HACDEMUXDIR}/pass/barcode11/*.fastq > ${SAMPLEID11}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode12" ] && cat ${HACDEMUXDIR}/pass/barcode12/*.fastq > ${SAMPLEID12}.fastq
source activate metagenomic
for i in $(find ./ -type f -name "*.fastq" | while read o; do basename $o | rev | cut -c 7- | rev; done | sort | uniq); do
	minimap2 -ax map-ont -t ${THREADS} ${REFSEQ} ${i}.fastq | samtools sort -@ 12 -o ${i}.sorted.bam -;
	samtools view -bS -f 4 ${i}.sorted.bam > ${i}.unmapped.sam -@ 12;
	samtools fastq -f 4 ${i}.unmapped.sam > ${i}.unmapped.fastq -@ 12;
	minimap2 -ax ava-ont -t ${THREADS} ${i}.unmapped.fastq ${i}.unmapped.fastq > {i}.overlap.sam;
	racon -t ${THREADS} -f -u ${i}.unmapped.fastq ${i}.overlap.sam ${i}.unmapped.fastq > .racon.fasta;
done