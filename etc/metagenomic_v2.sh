#!/bin/bash

# author: Laise de Moraes <laisepaixao@live.com>
# institution: Universidade Federal da Bahia, Brazil
# URL: https://github.com/lpmor22
# date: 07 APR 2021

# conda create -y -n pycoqc -c conda-forge -c bioconda -c defaults python=3.6 pycoQC
# mamba create -y -n pycoqc -c conda-forge -c bioconda -c defaults python=3.6 pycoQC
# conda create -y -n minimap2 -c conda-forge -c bioconda -c defaults minimap2 samtools
# mamba create -y -n minimap2 -c conda-forge -c bioconda -c defaults minimap2 samtools
# conda create -y -n racon -c conda-forge -c bioconda -c defaults nanopolish racon
# mamba create -y -n racon -c conda-forge -c bioconda -c defaults nanopolish racon

CONFIG="" #dna_r9.4.1_450bps_hac.cfg
ARRANGEMENTS="" #barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg
MINSCORE="" #8
TRIMADAPTER="" #18

THREADS="" #12
NUMCALLERS="" #NVIDIAGeForceRTX2060=60; NVIDIAGeForceRTX2080=64
GPUPERDEVICE="" #NVIDIAGeForceRTX2060=32; NVIDIAGeForceRTX2080=46

RAWDIR=""
HACBASECALLDIR=""
DEMUXDIR=""
CONCATENATED=""
READLEVELDIR=""
CONTIGLEVELDIR=""
RUNNAME=""

REFSEQ=""
# Homo sapiens: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13

guppy_basecaller -r -i ${RAWDIR} -s ${HACBASECALLDIR} -c ${CONFIG} --arrangements_files ${ARRANGEMENTS} --qscore_filtering --min_qscore ${MINSCORE} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --num_callers ${NUMCALLERS} --verbose_logs

guppy_barcoder -r -i ${HACBASECALLDIR} -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER} -t ${THREADS} -x auto

source activate pycoqc
pycoQC -q -f ${HACBASECALLDIR}/sequencing_summary.txt -b ${DEMUXDIR}/barcoding_summary.txt -o ${RUNNAME}_qc.html --report_title ${RUNNAME}

[ ! -d "${CONCATENATED}" ] && mkdir -vp ${CONCATENATED}
[ ! -d "${READLEVELDIR}" ] && mkdir -vp ${READLEVELDIR}

for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
	[ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > ${CONCATENATED}/${i}.fastq;
done

for i in $(find ${CONCATENATED} -type f -name "*.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	source activate minimap2
	minimap2 -ax map-ont -t ${THREADS} ${REFSEQ} ${CONCATENATED}/${i}.fastq | samtools sort -@ 12 -o ${READLEVELDIR}/${i}.sorted.bam -;
	samtools view -bS -f 4 ${READLEVELDIR}/${i}.sorted.bam > ${READLEVELDIR}/${i}.unmapped.sam -@ 12;
	samtools fastq -f 4 ${READLEVELDIR}/${i}.unmapped.sam > ${READLEVELDIR}/${i}.unmapped.fastq -@ 12;
	minimap2 -ax ava-ont -t ${THREADS} ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.overlap.sam;
	source activate racon
	racon -t ${THREADS} -f -u ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.overlap.sam ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.corrected.fasta;
done
