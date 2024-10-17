#!/bin/bash

THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

CONFIG="dna_r9.4.1_450bps_hac.cfg" #dna_r9.4.1_450bps_sup.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" #barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg

GPUPERDEVICE="64" # NVIDIAGeForceRTX2060=32; NVIDIAGeForceRTX2080=64
NUMCALLERS="46" # NVIDIAGeForceRTX2060=60; NVIDIAGeForceRTX2080=46

RAWDIR="/media/rklab/Data/library15-2019-11-15"
HACDEMUXDIR="/media/rklab/Data/Lb_BA144_HAC_Demux"
RUNNAME="Lb_BA144"

SAMPLEID0="Lb_BA144" # unclassified

guppy_basecaller -r -i ${RAWDIR} -s ${HACDEMUXDIR} -c ${CONFIG} --arrangements_files ${ARRANGEMENTS} --qscore_filtering --require_barcodes_both_ends --trim_barcodes --num_barcode_threads ${THREADS} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --num_callers ${NUMCALLERS} --verbose_logs

source activate pycoqc
pycoQC -q -f ${HACDEMUXDIR}/sequencing_summary.txt -o ${RUNNAME}_qc.html --report_title ${RUNNAME}

[ -d "${HACDEMUXDIR}/pass/unclassified" ] && cat ${HACDEMUXDIR}/pass/unclassified/*.fastq > ${SAMPLEID0}.fastq
