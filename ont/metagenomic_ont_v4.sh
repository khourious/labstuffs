#!/bin/bash

# Set the path of RAW fast5 files
RAWDIR=""

# Set the path of kraken2 database
KRAKEN2DB=""

MYSHELL=$(echo $SHELL | awk -F/ '{print $NF}')

# install GUPPY
if [[ -z "$(which guppy_basecaller)" ]]; then
    version=5.0.11
    cd
    curl https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_"$version"_linux64.tar.gz -o ont-guppy.tar.gz
    tar -vzxf ont-guppy.tar.gz
    rm -rf ont-guppy.tar.gz
    echo 'export PATH=$HOME/ont-guppy/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
    export PATH=$HOME/ont-guppy/bin:/usr/local/share/rsi/idl/bin:$PATH
else
    guppy_basecaller --version
fi

# install CONDA and create metagenomic environments
if [[ -z "$(which conda)" ]]; then
    cd
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -bfp miniconda3
    rm Miniconda3-latest-Linux-x86_64.sh
    echo 'export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
    export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH
    conda install -y -c conda-forge mamba
    mamba update -y -n base conda
    mamba create -y -n minimap2 -c conda-forge -c anaconda -c bioconda -c defaults minimap2 samtools
    mamba create -y -n plot -c conda-forge -c anaconda -c bioconda -c defaults pysam numpy pandas seaborn
    mamba create -y -n pycoqc -c aleg -c conda-forge -c anaconda -c bioconda -c defaults python=3.6 pycoqc
    mamba create -y -n racon -c conda-forge -c bioconda -c defaults nanopolish racon
else
    if [[ -z "$(which mamba)" ]]; then
        conda install -y -c conda-forge mamba
        mamba update -y -n base conda
        mamba create -y -n minimap2 -c conda-forge -c anaconda -c bioconda -c defaults minimap2 samtools
        mamba create -y -n plot -c conda-forge -c anaconda -c bioconda -c defaults pysam numpy pandas seaborn
        mamba create -y -n pycoqc -c aleg -c conda-forge -c anaconda -c bioconda -c defaults python=3.6 pycoqc
        mamba create -y -n racon -c conda-forge -c bioconda -c defaults nanopolish racon
    else
        mamba update -y -n base conda
        mamba create -y -n minimap2 -c conda-forge -c anaconda -c bioconda -c defaults minimap2 samtools
        mamba create -y -n plot -c conda-forge -c anaconda -c bioconda -c defaults pysam numpy pandas seaborn
        mamba create -y -n pycoqc -c aleg -c conda-forge -c anaconda -c bioconda -c defaults python=3.6 pycoqc
        mamba create -y -n racon -c conda-forge -c bioconda -c defaults nanopolish racon
    fi
fi

bg() {

    start=$(date +%s.%N)

    THREADS=$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')

    CONFIG=dna_r9.4.1_450bps_sup.cfg #dna_r9.4.1_450bps_hac.cfg
    ARRANGEMENTS=barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg

    TRIMADAPTER=18

    RUNNAME=$(basename "$RAWDIR")
    REFSEQS="$RAWDIR"/REFSEQS
    BASECALLDIR="$RAWDIR"/BASECALL
    DEMUXDIR="$RAWDIR"/DEMUX
    DEMUXCATDIR="$RAWDIR"/DEMUX_CAT
    READLEVELDIR="$RAWDIR"/LEVEL_READS
    CONTIGLEVELDIR="$RAWDIR"/LEVEL_CONTIGS
    ASSEMBLYDIR="$RAWDIR"/ASSEMBLY

    [ ! -d "$REFSEQS" ] && mkdir -vp "$REFSEQS"
    [ ! -d "$DEMUXCATDIR" ] && mkdir -vp "$DEMUXCATDIR"
    [ ! -d "$READLEVELDIR" ] && mkdir -vp "$READLEVELDIR"
    [ ! -d "$ASSEMBLYDIR" ] && mkdir -vp "$ASSEMBLYDIR"

    wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -q \
    -O "$RAWDIR"/REFSEQS/GRCh38.p13.genome.fa.gz

    HUMANREFSEQ="$RAWDIR"/REFSEQS/GRCh38.p13.genome.fa.gz

    guppy_basecaller -r -x auto --verbose_logs -c "$CONFIG" -i "$RAWDIR" -s "$BASECALLDIR" \
    -q 1 --min_qscore 10 --chunk_size 1000 --num_callers "$THREADS" --gpu_runners_per_device 1

    end=$(date +%s.%N)

    runtime=$(python -c "print(${end} - ${start})")

    echo "" && echo "Done. The runtime was $runtime seconds." && echo ""

}

bg &>>metagenomic_ont_v4_log.txt &


    # guppy_barcoder -r -i ${BASECALLDIR} -s ${DEMUXDIR} --arrangements_files ${ARRANGEMENTS} --require_barcodes_both_ends --trim_barcodes --num_extra_bases_trim ${TRIMADAPTER} -t ${THREADS} -x auto

    # source activate pycoqc

    # pycoQC -q -f ${BASECALLDIR}/sequencing_summary.txt -b ${DEMUXDIR}/barcoding_summary.txt -o ${RAWDIR}/${RUNNAME}_qc.html --report_title ${RUNNAME}


# for i in $(find ${DEMUXDIR} -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
    # [ -d "${DEMUXDIR}/${i}" ] && cat ${DEMUXDIR}/${i}/*.fastq > ${DEMUXCATDIR}/${i}.fastq
# done

# for i in $(find ${DEMUXCATDIR} -type f -name "*.fastq" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
    # source activate minimap2
    # minimap2 -ax map-ont -t ${THREADS} ${HUMANREFSEQ} ${DEMUXCATDIR}/${i}.fastq | samtools sort -@ 12 -o ${READLEVELDIR}/${i}.sorted.bam -
    # samtools index -@ ${THREADS} ${READLEVELDIR}/${i}.sorted.bam
    # samtools view -bS -f 4 ${READLEVELDIR}/${i}.sorted.bam > ${READLEVELDIR}/${i}.unmapped.sam -@ 12
    # samtools fastq -f 4 ${READLEVELDIR}/${i}.unmapped.sam > ${READLEVELDIR}/${i}.unmapped.fastq -@ 12
    # minimap2 -ax ava-ont -t ${THREADS} ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.overlap.sam
    # source activate racon
    # racon -t ${THREADS} -f -u ${READLEVELDIR}/${i}.unmapped.fastq ${READLEVELDIR}/${i}.overlap.sam ${READLEVELDIR}/${i}.unmapped.fastq > ${READLEVELDIR}/${i}.corrected.fasta
# done

# for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
	# kraken2 --db ${KRAKENDB} --threads ${THREADS} --report ${READLEVELDIR}/${i}_report.txt --report-minimizer-data --output ${READLEVELDIR}/${i}_output.txt ${READLEVELDIR}/${i}.corrected.fasta
# done

# CHIKVREFSEQ=""
# DENV1REFSEQ=""
# DENV2REFSEQ=""
# DENV3REFSEQ=""
# DENV4REFSEQ=""
# ZIKVREFSEQ=""

#source activate minimap2

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${CHIKVREFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.chikv.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.chikv.sorted.bam > ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${CHIKVREFSEQ} ${ASSEMBLYDIR}/${i}.chikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.chikv -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV1REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv1.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv1.sorted.bam > ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV1REFSEQ} ${ASSEMBLYDIR}/${i}.denv1.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv1 -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV2REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv2.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv2.sorted.bam > ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV2REFSEQ} ${ASSEMBLYDIR}/${i}.denv2.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv2 -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV3REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv3.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv3.sorted.bam > ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV3REFSEQ} ${ASSEMBLYDIR}/${i}.denv3.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv3 -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${DENV4REFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.denv4.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.denv4.sorted.bam > ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${DENV4REFSEQ} ${ASSEMBLYDIR}/${i}.denv4.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.denv4 -n N -i ${i}
#done

#for i in $(find ${READLEVELDIR} -type f -name "*.fasta" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	minimap2 -t ${THREADS} -ax map-ont ${ZIKVREFSEQ} ${READLEVELDIR}/${i}.corrected.fasta | samtools sort -@ ${THREADS} -o ${ASSEMBLYDIR}/${i}.zikv.sorted.bam -
#	samtools view -@ ${THREADS} -h -F 4 -b ${ASSEMBLYDIR}/${i}.zikv.sorted.bam > ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
#	samtools index -@ ${THREADS} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam
#	samtools mpileup -A -B -Q 0 --reference ${ZIKVREFSEQ} ${ASSEMBLYDIR}/${i}.zikv.sorted.mapped.bam | ivar consensus -p ${ASSEMBLYDIR}/${i}.zikv -n N -i ${i}
#done

#source activate plot

#for i in $(find ${ASSEMBLYDIR} -type f -name "*.sorted.mapped.bam" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do
#	fastcov ${ASSEMBLYDIR}/barcode*.chikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_chikv.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.chikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_chikv_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.denv1.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv1.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.denv1.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv1_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.denv2.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv2.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.denv2.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv2_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.denv3.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv3.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.denv3.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv3_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.denv4.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv4.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.denv4.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_denv4_log.pdf
#	fastcov ${ASSEMBLYDIR}/barcode*.zikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_zikv.pdf
#	fastcov -l ${ASSEMBLYDIR}/barcode*.zikv.sorted.mapped.bam -o ${ASSEMBLYDIR}/assembly_zikv_log.pdf
#done
