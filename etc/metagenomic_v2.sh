#!/bin/bash

# author: Laise de Moraes <laisepaixao@live.com>
# institution: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# URL: https://lpmor22.github.io
# date: 30 MAY 2021

CONFIG="dna_r9.4.1_450bps_sup.cfg" #dna_r9.4.1_450bps_sup.cfg
ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" #barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg
MINSCORE="10" #8
TRIMADAPTER="18" #18

THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"
NUMCALLERS="60" #NVIDIAGeForceRTX1650=64; NVIDIAGeForceRTX2060=60; NVIDIAGeForceRTX2080=64
GPUPERDEVICE="32" #NVIDIAGeForceRTX1650=16; NVIDIAGeForceRTX2060=32; NVIDIAGeForceRTX2080=46

RAWDIR="/mnt/c/Dropbox/NGS_Kalabrick/NGS_LIBRARY02_20210521/NGS_LIBRARY02_20210521"
SUPBASECALLDIR="/mnt/c/Dropbox/NGS_Kalabrick/NGS_LIBRARY02_20210521/SUP_BASECALL"
DEMUXDIR="/mnt/c/Dropbox/NGS_Kalabrick/NGS_LIBRARY02_20210521/DEMUX"
DEMUXCATDIR="/mnt/c/Dropbox/NGS_Kalabrick/NGS_LIBRARY02_20210521/DEMUX_CAT"
READLEVELDIR="/mnt/c/Dropbox/NGS_Kalabrick/NGS_LIBRARY02_20210521/LEVEL_READS"
CONTIGLEVELDIR="/mnt/c/Dropbox/NGS_Kalabrick/NGS_LIBRARY02_20210521/LEVEL_CONTIGS"
ASSEMBLYDIR="/mnt/c/Dropbox/NGS_Kalabrick/NGS_LIBRARY02_20210521/ASSEMBLY"
RUNNAME="NGS_LIBRARY02_20210521"

# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz
HUMANREFSEQ="/mnt/c/Dropbox/NGS_Kalabrick/REFS/GRCh38.p13.genome.fa.gz"

CHIKVREFSEQ="/mnt/c/Dropbox/NGS_Kalabrick/REFS/CHIKV_ASIAN.fasta"
DENV1REFSEQ="/mnt/c/Dropbox/NGS_Kalabrick/REFS/DENV_1.fasta"
DENV2REFSEQ="/mnt/c/Dropbox/NGS_Kalabrick/REFS/DENV_2.fasta"
DENV3REFSEQ="/mnt/c/Dropbox/NGS_Kalabrick/REFS/DENV_3.fasta"
DENV4REFSEQ="/mnt/c/Dropbox/NGS_Kalabrick/REFS/DENV_4.fasta"
ZIKVREFSEQ="/mnt/c/Dropbox/NGS_Kalabrick/REFS/ZIKV_ASIAN.fasta"

guppy_basecaller -r -i ${RAWDIR} -s ${SUPBASECALLDIR} -c ${CONFIG} --arrangements_files ${ARRANGEMENTS} --min_qscore ${MINSCORE} -x auto --gpu_runners_per_device ${GPUPERDEVICE} --num_callers ${NUMCALLERS} --verbose_logs

guppy_barcoder -r -i ${SUPBASECALLDIR} -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER} -t ${THREADS} -x auto

source activate pycoqc

pycoQC -q -f ${SUPBASECALLDIR}/sequencing_summary.txt -b ${DEMUXDIR}/barcoding_summary.txt -o ${RUNNAME}_qc.html --report_title ${RUNNAME}

[ ! -d "${DEMUXCATDIR}" ] && mkdir -vp ${DEMUXCATDIR}
[ ! -d "${READLEVELDIR}" ] && mkdir -vp ${READLEVELDIR}
[ ! -d "${ASSEMBLYDIR}" ] && mkdir -vp ${ASSEMBLYDIR}

for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
	[ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > ${DEMUXCATDIR}/${i}.fastq
done

for i in $(find ${DEMUXCATDIR} -type f -name "*.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	source activate minimap2
	minimap2 -ax map-ont -t ${THREADS} ${HUMANREFSEQ} ${DEMUXCATDIR}/${i}.fastq | samtools sort -@ 12 -o ${READLEVELDIR}/${i}.sorted.bam -
	samtools index -@ ${THREADS} ${READLEVELDIR}/${i}.sorted.bam
	samtools view -bS -f 4 ${READLEVELDIR}/${i}.sorted.bam > ${READLEVELDIR}/${i}.unmapped.sam -@ 12
	samtools fastq -f 4 ${READLEVELDIR}/${i}.unmapped.sam > ${READLEVELDIR}/${i}.unmapped.fastq -@ 12
	minimap2 -ax ava-ont -t ${THREADS} ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.overlap.sam
	source activate racon
	racon -t ${THREADS} -f -u ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.overlap.sam ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.corrected.fasta
done

for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	kraken2 --db /mnt/x/kraken2db/k2_standard_8gb_20201202 --threads ${THREADS} --report ${READLEVELDIR}/${i}_report.txt --report-minimizer-data --output ${READLEVELDIR}/${i}_output.txt ${READLEVELDIR}/${i}.corrected.fasta
done

source activate minimap2	

for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	minimap2 -t ${THREADS} -ax map-ont ${CHIKVREFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.chikv.sorted.bam -
	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.chikv.sorted.bam > ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
	samtools mpileup -A -B -Q 0 --reference ${CHIKVREFSEQ} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.chikv -n N -i ${i}
done

for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	minimap2 -t ${THREADS} -ax map-ont ${DENV1REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv1.sorted.bam -
	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv1.sorted.bam > ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
	samtools mpileup -A -B -Q 0 --reference ${DENV1REFSEQ} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv1 -n N -i ${i}
done

for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	minimap2 -t ${THREADS} -ax map-ont ${DENV2REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv2.sorted.bam -
	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv2.sorted.bam > ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
	samtools mpileup -A -B -Q 0 --reference ${DENV2REFSEQ} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv2 -n N -i ${i}
done

for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	minimap2 -t ${THREADS} -ax map-ont ${DENV3REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv3.sorted.bam -
	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv3.sorted.bam > ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
	samtools mpileup -A -B -Q 0 --reference ${DENV3REFSEQ} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv3 -n N -i ${i}
done

for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	minimap2 -t ${THREADS} -ax map-ont ${DENV4REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv4.sorted.bam -
	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv4.sorted.bam > ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
	samtools mpileup -A -B -Q 0 --reference ${DENV4REFSEQ} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv4 -n N -i ${i}
done

for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	minimap2 -t ${THREADS} -ax map-ont ${ZIKVREFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.zikv.sorted.bam -
	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.zikv.sorted.bam > ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
	samtools mpileup -A -B -Q 0 --reference ${ZIKVREFSEQ} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.zikv -n N -i ${i}
done

source activate plot

for i in $(find ${ASSEMBLYDIR} -type f -name "*.sorted.mapped.bam" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	fastcov ${ASSEMBLYDIR}/barcode*.chikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_chikv.pdf
	fastcov -l ${ASSEMBLYDIR}/barcode*.chikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_chikv_log.pdf
	fastcov ${ASSEMBLYDIR}/barcode*.denv1.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv1.pdf
	fastcov -l ${ASSEMBLYDIR}/barcode*.denv1.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv1_log.pdf
	fastcov ${ASSEMBLYDIR}/barcode*.denv2.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv2.pdf
	fastcov -l ${ASSEMBLYDIR}/barcode*.denv2.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv2_log.pdf
	fastcov ${ASSEMBLYDIR}/barcode*.denv3.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv3.pdf
	fastcov -l ${ASSEMBLYDIR}/barcode*.denv3.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv3_log.pdf
	fastcov ${ASSEMBLYDIR}/barcode*.denv4.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv4.pdf
	fastcov -l ${ASSEMBLYDIR}/barcode*.denv4.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv4_log.pdf
	fastcov ${ASSEMBLYDIR}/barcode*.zikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_zikv.pdf
	fastcov -l ${ASSEMBLYDIR}/barcode*.zikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_zikv_log.pdf
done