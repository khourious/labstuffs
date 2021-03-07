#!/bin/bash

RAWDIR="/path/directory"
HACDEMUXDIR="/path/directory"
RUNNAME="run"
SAMPLELISTDIR="/path/directory/sample_list.txt"

THREADS="12"
GPUPERDEVICE="64" # NVIDIAGeForceRTX2060=32; NVIDIAGeForceRTX2080=64
NUMCALLERS="46" # NVIDIAGeForceRTX2060=60; NVIDIAGeForceRTX2080=46

guppy_basecaller -r -i ${RAWDIR} -s ${HACDEMUXDIR} -c dna_r9.4.1_450bps_hac.cfg --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" --require_barcodes_both_ends --qscore_filtering --num_barcode_threads ${THREADS} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --num_callers ${NUMCALLERS} --verbose_logs

source activate pycoqc
pycoQC -q -f ${HACDEMUXDIR}/sequencing_summary.txt -o ${RUNNAME}_qc.html --report_tittle ${RUNNAME}

source activate artic-ncov2019
for barcode in $(find ${HACDEMUXDIR}/pass -type d -name "barcode*" | while read output; do basename $output; done); do artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ${HACDEMUXDIR}/pass/$barcode --prefix ${RUNNAME}; done
cat ${SAMPLELISTDIR} | for file in $(find -type f -name "*.fastq*" | while read output; do basename $output; done | rev | cut -c 7- | rev | sort | uniq); do read line; mv -v "$file".fastq "$line".fastq; done
for sample in $(find -type f -name "*.fastq" | rev | cut -c 7- | rev); do artic minion --normalise 200 --threads ${THREADS} --scheme-directory $HOME/softwares/artic-ncov2019/primer_schemes --read-file $sample.fastq --fast5-directory ${RAWDIR} --sequencing-summary ${HACDEMUXDIR}/sequencing_summary.txt nCoV-2019/V3 $sample --medaka; done
cat *.consensus.fasta > ${RUNNAME}_consensus_genomes.fasta