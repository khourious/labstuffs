#!/bin/bash

# author: Laise de Moraes <laisepaixao@live.com>
# institution: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# URL: https://lpmor22.github.io
# date: 29 MAY 2021

CONFIG="dna_r9.4.1_450bps_hac.cfg" #dna_r9.4.1_450bps_hac.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" #barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg

THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}')"
GPUPERDEVICE="64" # NVIDIAGeForceRTX2060=32; NVIDIAGeForceRTX2080=64
NUMCALLERS="46" # NVIDIAGeForceRTX2060=60; NVIDIAGeForceRTX2080=46

RAWDIR="/path/directory"
HACDEMUXDIR="/path/directory"
RUNNAME="run"

SAMPLEID0="" # unclassified
SAMPLEID1=""
SAMPLEID2=""
SAMPLEID3=""
SAMPLEID4=""
SAMPLEID5=""
SAMPLEID6=""
SAMPLEID7=""
SAMPLEID8=""
SAMPLEID9=""
SAMPLEID10=""
SAMPLEID11=""
SAMPLEID12=""
SAMPLEID13=""
SAMPLEID14=""
SAMPLEID15=""
SAMPLEID16=""
SAMPLEID17=""
SAMPLEID18=""
SAMPLEID19=""
SAMPLEID20=""
SAMPLEID21=""
SAMPLEID22=""
SAMPLEID23=""
SAMPLEID24=""

guppy_basecaller -r -i ${RAWDIR} -s ${HACDEMUXDIR} -c ${CONFIG} --arrangements_files ${ARRANGEMENTS} --qscore_filtering --require_barcodes_both_ends --trim_barcodes --num_barcode_threads ${THREADS} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --num_callers ${NUMCALLERS} --verbose_logs
source activate pycoqc
pycoQC -q -f ${HACDEMUXDIR}/sequencing_summary.txt -o ${RUNNAME}_qc.html --report_title ${RUNNAME}
[ -d "${HACDEMUXDIR}/pass/unclassified" ] && cat ${HACDEMUXDIR}/pass/unclassified/*.fastq > ${SAMPLEID0}.fastq
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
[ -d "${HACDEMUXDIR}/pass/barcode13" ] && cat ${HACDEMUXDIR}/pass/barcode13/*.fastq > ${SAMPLEID13}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode14" ] && cat ${HACDEMUXDIR}/pass/barcode14/*.fastq > ${SAMPLEID14}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode15" ] && cat ${HACDEMUXDIR}/pass/barcode15/*.fastq > ${SAMPLEID15}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode16" ] && cat ${HACDEMUXDIR}/pass/barcode16/*.fastq > ${SAMPLEID16}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode17" ] && cat ${HACDEMUXDIR}/pass/barcode17/*.fastq > ${SAMPLEID17}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode18" ] && cat ${HACDEMUXDIR}/pass/barcode18/*.fastq > ${SAMPLEID18}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode19" ] && cat ${HACDEMUXDIR}/pass/barcode19/*.fastq > ${SAMPLEID19}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode20" ] && cat ${HACDEMUXDIR}/pass/barcode20/*.fastq > ${SAMPLEID20}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode21" ] && cat ${HACDEMUXDIR}/pass/barcode21/*.fastq > ${SAMPLEID21}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode22" ] && cat ${HACDEMUXDIR}/pass/barcode22/*.fastq > ${SAMPLEID22}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode23" ] && cat ${HACDEMUXDIR}/pass/barcode23/*.fastq > ${SAMPLEID23}.fastq
[ -d "${HACDEMUXDIR}/pass/barcode24" ] && cat ${HACDEMUXDIR}/pass/barcode24/*.fastq > ${SAMPLEID24}.fastq
