#!/bin/bash

# conda create -y -n pycoqc -c conda-forge -c bioconda -c defaults python=3.6 pycoQC
# mamba create -y -n pycoqc -c conda-forge -c bioconda -c defaults python=3.6 pycoQC
# conda create -y -n minimap2 -c conda-forge -c bioconda -c defaults minimap2 samtools
# mamba create -y -n minimap2 -c conda-forge -c bioconda -c defaults minimap2 samtools
# conda create -y -n racon -c conda-forge -c bioconda -c defaults nanopolish racon
# mamba create -y -n racon -c conda-forge -c bioconda -c defaults nanopolish racon

CONFIG="dna_r9.4.1_450bps_hac.cfg" #dna_r9.4.1_450bps_hac.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" #barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg
MINSCORE="8"
TRIMADAPTER="18"

THREADS="12"
NUMCALLERS="60" #NVIDIAGeForceRTX2060=60; NVIDIAGeForceRTX2080=46
GPUPERDEVICE="32" #NVIDIAGeForceRTX2060=32; NVIDIAGeForceRTX2080=64

RAWDIR="/media/lpmor22/2019-12-19/DENV_Seq_Data/DENV_FTA_1/DENV_Run1_data"
HACBASECALLDIR="/media/lpmor22/2019-12-19/DENV_Seq_Data/DENV_FTA_1/HAC_BASECALL"
DEMUXDIR="/media/lpmor22/2019-12-19/DENV_Seq_Data/DENV_FTA_1/DEMUX"
CONCATENATED="/media/lpmor22/2019-12-19/DENV_Seq_Data/DENV_FTA_1/DEMUX_concatenated"
READLEVELDIR="/media/lpmor22/2019-12-19/DENV_Seq_Data/DENV_FTA_1/READ_level_analysis"
CONTIGLEVELDIR="/media/lpmor22/2019-12-19/DENV_Seq_Data/DENV_FTA_1/CONTIG_level_analysis"
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

REFSEQ="/media/lpmor22/2019-12-19/DENV_Seq_Data/DENV_FTA_1/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13

guppy_basecaller -r -i ${RAWDIR} -s ${HACBASECALLDIR} -c ${CONFIG} --arrangements_files ${ARRANGEMENTS} --qscore_filtering --min_qscore ${MINSCORE} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --num_callers ${NUMCALLERS} --verbose_logs

guppy_barcoder -r -i ${HACBASECALLDIR} -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER} -t ${THREADS} -x auto

source activate pycoqc
pycoQC -q -f ${HACBASECALLDIR}/sequencing_summary.txt -b ${DEMUXDIR}/barcoding_summary.txt -o ${RUNNAME}_qc.html --report_title ${RUNNAME}

[ ! -d "${CONCATENATED}" ] && mkdir -vp ${CONCATENATED}
[ ! -d "${READLEVELDIR}" ] && mkdir -vp ${READLEVELDIR}
[ ! -d "${CONTIGLEVELDIR}" ] && mkdir -vp ${CONTIGLEVELDIR}

[ -d "${DEMUXDIR}/barcode01" ] && cat ${DEMUXDIR}/barcode01/*.fastq > ${CONCATENATED}/${SAMPLEID1}.fastq
[ -d "${DEMUXDIR}/barcode02" ] && cat ${DEMUXDIR}/barcode02/*.fastq > ${CONCATENATED}/${SAMPLEID2}.fastq
[ -d "${DEMUXDIR}/barcode03" ] && cat ${DEMUXDIR}/barcode03/*.fastq > ${CONCATENATED}/${SAMPLEID3}.fastq
[ -d "${DEMUXDIR}/barcode04" ] && cat ${DEMUXDIR}/barcode04/*.fastq > ${CONCATENATED}/${SAMPLEID4}.fastq
[ -d "${DEMUXDIR}/barcode05" ] && cat ${DEMUXDIR}/barcode05/*.fastq > ${CONCATENATED}/${SAMPLEID5}.fastq
[ -d "${DEMUXDIR}/barcode06" ] && cat ${DEMUXDIR}/barcode06/*.fastq > ${CONCATENATED}/${SAMPLEID6}.fastq
[ -d "${DEMUXDIR}/barcode07" ] && cat ${DEMUXDIR}/barcode07/*.fastq > ${CONCATENATED}/${SAMPLEID7}.fastq
[ -d "${DEMUXDIR}/barcode08" ] && cat ${DEMUXDIR}/barcode08/*.fastq > ${CONCATENATED}/${SAMPLEID8}.fastq
[ -d "${DEMUXDIR}/barcode09" ] && cat ${DEMUXDIR}/barcode09/*.fastq > ${CONCATENATED}/${SAMPLEID9}.fastq
[ -d "${DEMUXDIR}/barcode10" ] && cat ${DEMUXDIR}/barcode10/*.fastq > ${CONCATENATED}/${SAMPLEID10}.fastq
[ -d "${DEMUXDIR}/barcode11" ] && cat ${DEMUXDIR}/barcode11/*.fastq > ${CONCATENATED}/${SAMPLEID11}.fastq
[ -d "${DEMUXDIR}/barcode12" ] && cat ${DEMUXDIR}/barcode12/*.fastq > ${CONCATENATED}/${SAMPLEID12}.fastq

for i in $(find ${CONCATENATED} -type f -name "*.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	source activate minimap2
	minimap2 -ax map-ont -t ${THREADS} ${REFSEQ} ${CONCATENATED}/${i}.fastq | samtools sort -@ 12 -o ${READLEVELDIR}/${i}.sorted.bam -;
	samtools view -bS -f 4 ${READLEVELDIR}/${i}.sorted.bam > ${READLEVELDIR}/${i}.unmapped.sam -@ 12;
	samtools fastq -f 4 ${READLEVELDIR}/${i}.unmapped.sam > ${READLEVELDIR}/${i}.unmapped.fastq -@ 12;
	minimap2 -ax ava-ont -t ${THREADS} ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.overlap.sam;
	source activate racon
	racon -t ${THREADS} -f -u ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.overlap.sam ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.corrected.fasta;
done