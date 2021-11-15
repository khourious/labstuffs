#!/bin/bash

while getopts ":i:l:s:t:" opt; do
        case $opt in
                i ) INPUT="$OPTARG";;
                l ) LEVEL="$OPTARG";;
                s ) SAMPLESHEET="$OPTARG";;
                t ) THREADS="$OPTARG";;
        esac
done
shift $((OPTIND -1))

# LEVEL 0 = install conda + mamba + enviroments + subread
# LEVEL 1 = copy RAW files + fastQC/multiQC lanes
# LEVEL 2 = cat lanes + fastQC/multiQC runs
# LEVEL 3 = cat runs + fastQC/multiQC merged runs
# LEVEL 4 = trimmomatic
# LEVEL 5 = STAR
# LEVEL 6 = featureCounts

MYSHELL=$(echo $SHELL | awk -F/ '{print $NF}')

if [[ -z "$THREADS" ]]; then
    THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"
fi

if [[ "$LEVEL" = "0" ]]; then
    if [[ -z "$(which conda)" ]]; then
        cd
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -bfp miniconda3
        rm Miniconda3-latest-Linux-x86_64.sh
        echo 'export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
        export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH
        conda install -y -c conda-forge mamba
        mamba update -y -c conda-forge -c anaconda -c bioconda -c defaults -n base conda
        mamba create -y -n r -c r -c conda-forge -c anaconda -c bioconda -c defaults bioconductor-deseq2 bioconductor-enhancedvolcano bioconductor-genefilter bioconductor-vsn boost-cpp r-ashr r-base=4.1 r-biocmanager r-calibrate r-devtools r-ggplot2 r-gplots r-pheatmap r-rcolorbrewer r-remotes
        mamba create -y -n rnaseq -c conda-forge -c anaconda -c bioconda -c defaults fastqc multiqc trimmomatic star
    else
        if [[ -z "$(which mamba)" ]]; then
            conda install -y -c conda-forge mamba
            mamba update -y -c conda-forge -c anaconda -c bioconda -c defaults -n base conda
            if [[ -z "$(conda env list | grep "r|rnaseq")" ]]; then
                mamba create -y -n r -c r -c conda-forge -c anaconda -c bioconda -c defaults bioconductor-deseq2 bioconductor-enhancedvolcano bioconductor-genefilter bioconductor-vsn boost-cpp r-ashr r-base=4.1 r-biocmanager r-calibrate r-devtools r-ggplot2 r-gplots r-pheatmap r-rcolorbrewer r-remotes
                mamba create -y -n rnaseq -c conda-forge -c anaconda -c bioconda -c defaults fastqc multiqc trimmomatic star
            else
                mamba update -y -n r -c r -c conda-forge -c anaconda -c bioconda -c defaults --all
                mamba update -y -n rnaseq -c conda-forge -c anaconda -c bioconda -c defaults --all
            fi
        else
            mamba update -y -c conda-forge -c anaconda -c bioconda -c defaults -n base conda
            if [[ -z "$(conda env list | grep "r|rnaseq")" ]]; then
                mamba create -y -n r -c r -c conda-forge -c anaconda -c bioconda -c defaults bioconductor-deseq2 bioconductor-enhancedvolcano bioconductor-genefilter bioconductor-vsn boost-cpp r-ashr r-base=4.1 r-biocmanager r-calibrate r-devtools r-ggplot2 r-gplots r-pheatmap r-rcolorbrewer r-remotes
                mamba create -y -n rnaseq -c conda-forge -c anaconda -c bioconda -c defaults fastqc multiqc trimmomatic star
            else
                mamba update -y -n r -c r -c conda-forge -c anaconda -c bioconda -c defaults --all
                mamba update -y -n rnaseq -c conda-forge -c anaconda -c bioconda -c defaults --all
            fi
        fi
    fi
    if [[ -z "$(which featureCounts)" ]]; then
        version=2.0.3
        cd
        wget https://ufpr.dl.sourceforge.net/project/subread/subread-"$version"/subread-"$version"-source.tar.gz
        tar -xvf subread-"$version"-source.tar.gz
        rm subread-"$version"-source.tar.gz
        mv subread-"$version"-source subread-source
        cd subread-source/src
        make -f Makefile.Linux
        echo 'export PATH=$HOME/subread-source/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
        export PATH=$HOME/subread-source/bin:/usr/local/share/rsi/idl/bin:$PATH
    else
        featureCounts -v
    fi
fi

bg() {

    start=$(date +%s.%N)

    source activate rnaseq

    if [[ -n "$INPUT" ]]; then
        RAWDIR="$INPUT"
        RAWRUNDIR1="$RAWDIR"/FASTQ_Generation_2018-08-17_15_33_27Z-116577544
        RAWRUNDIR2="$RAWDIR"/FASTQ_Generation_2018-10-17_13_51_20Z-130391988
        RAWRUNDIR3="$RAWDIR"/FASTQ_Generation_2018-12-11_13_17_10Z-141945716
        RAWRUNDIR4="$RAWDIR"/FASTQ_Generation_2018-12-13_14_45_41Z-143225097
        RAWRUNDIR5="$RAWDIR"/FASTQ_Generation_2018-12-15_10_26_29Z-143716575
        ANALYSIS="$RAWDIR"/ANALYSIS
        if [[ ! -d "$ANALYSIS" ]]; then
            mkdir "$ANALYSIS" "$ANALYSIS"/ALIGN "$ANALYSIS"/INDEX "$ANALYSIS"/QC_RUN{1..5} "$ANALYSIS"/QC_RUNS_MERGED \
            "$ANALYSIS"/REFERENCE "$ANALYSIS"/RUN{1..5} "$ANALYSIS"/RUNS_MERGED "$ANALYSIS"/TRIMMED
        fi
    fi

    if [[ "$LEVEL" == 1 ]]; then
        find "$RAWRUNDIR1" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN1 {} +
        find "$RAWRUNDIR2" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN2 {} +
        find "$RAWRUNDIR3" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN3 {} +
        find "$RAWRUNDIR4" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN4 {} +
        find "$RAWRUNDIR5" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN5 {} +
        for i in $(find "$ANALYSIS" -maxdepth 1 -mindepth 1 -type d -name "RUN*[1-5]" | awk -F/ '{print $NF}'); do
            fastqc -t "$THREADS" "$ANALYSIS"/"$i"/*.fastq.gz -o "$ANALYSIS"/QC_"$i"
            multiqc -s -i "ArbovirusFiocruzBA-83677594 "$i" lanes" -ip --no-data-dir -n "$ANALYSIS"/"$i"_lanes_multiqc_report "$ANALYSIS"/QC_"$i"/*
        done
    fi

    if [[ "$LEVEL" == 2 ]]; then
        for i in $(find "$ANALYSIS"/RUN{1..5} -type f -name "*.fastq.gz" | awk -F/ '{print $NF}' | awk -F_ '{print $1"_"$2}' | sort -u); do
            cat "$ANALYSIS"/RUN1/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN1/"$i"_RUN1_R1.fastq.gz
            cat "$ANALYSIS"/RUN1/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN1/"$i"_RUN1_R2.fastq.gz
            cat "$ANALYSIS"/RUN2/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN2/"$i"_RUN2_R1.fastq.gz
            cat "$ANALYSIS"/RUN2/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN2/"$i"_RUN2_R2.fastq.gz
            cat "$ANALYSIS"/RUN3/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN3/"$i"_RUN3_R1.fastq.gz
            cat "$ANALYSIS"/RUN3/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN3/"$i"_RUN3_R2.fastq.gz
            cat "$ANALYSIS"/RUN4/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN4/"$i"_RUN4_R1.fastq.gz
            cat "$ANALYSIS"/RUN4/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN4/"$i"_RUN4_R2.fastq.gz
            cat "$ANALYSIS"/RUN5/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN5/"$i"_RUN5_R1.fastq.gz
            cat "$ANALYSIS"/RUN5/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN5/"$i"_RUN5_R2.fastq.gz
        done
        rm -rf "$ANALYSIS"/RUN{1..5}/*_001.fastq.gz
        for i in $(find "$ANALYSIS" -maxdepth 1 -mindepth 1 -type d -name "RUN*[1-5]" | awk -F/ '{print $NF}'); do
            fastqc -t "$THREADS" "$ANALYSIS"/"$i"/*RUN* -o "$ANALYSIS"/QC_"$i"
            multiqc -s -i "$i" -b "ArbovirusFiocruzBA-83677594 "$i" merged lanes" -ip --no-data-dir \
            -n "$ANALYSIS"/"$i"_merged_lanes_multiqc_report "$ANALYSIS"/QC_"$i"/*RUN*
        done
    fi

    if [[ "$LEVEL" == 3 ]]; then
        for i in $(find "$ANALYSIS"/RUN1 -type f -name "*.fastq.gz" | awk -F/ '{print $NF}' | awk -F_ '{print $1"_"$2}' | sort -u); do
            cat "$ANALYSIS"/RUN*/"$i"_RUN*_R1.fastq.gz > "$ANALYSIS"/RUNS_MERGED/"$i"_MERGED_R1.fastq.gz
            cat "$ANALYSIS"/RUN*/"$i"_RUN*_R2.fastq.gz > "$ANALYSIS"/RUNS_MERGED/"$i"_MERGED_R2.fastq.gz
        done
        rm -rf "$ANALYSIS"/RUN{1..5}
        fastqc -t "$THREADS" "$ANALYSIS"/RUNS_MERGED/*.fastq.gz -o "$ANALYSIS"/QC_RUNS_MERGED
        multiqc -s -i "ArbovirusFiocruzBA-83677594 RUNS_MERGED" -ip --no-data-dir \
        -n "$ANALYSIS"/RUNS_MERGED_multiqc_report "$ANALYSIS"/QC_RUNS_MERGED/*MERGED*
    fi

    if [[ "$LEVEL" == 4 ]]; then
        for i in $(find "$ANALYSIS"/RUNS_MERGED -type f -name "*.fastq.gz" | awk -F/ '{print $NF}' | awk -F_ '{print $1"_"$2}' | sort -u); do
            trimmomatic PE -threads "$THREADS" -phred33 \
            "$ANALYSIS"/RUNS_MERGED/"$i"_MERGED_R1.fastq.gz "$ANALYSIS"/RUNS_MERGED/"$i"_MERGED_R2.fastq.gz \
            "$ANALYSIS"/TRIMMED/"$i"_R1_PAIRED.fastq.gz "$ANALYSIS"/TRIMMED/"$i"_R1_UNPAIRED.fastq.gz \
            "$ANALYSIS"/TRIMMED/"$i"_R2_PAIRED.fastq.gz "$ANALYSIS"/TRIMMED/"$i"_R2_UNPAIRED.fastq.gz \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36
        done
    fi

    if [[ "$LEVEL" == 5 ]]; then
        wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -q \
        -O "$ANALYSIS"/REFERENCE/GRCh38.p13.genome.fa.gz
        wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz -q \
        -O "$ANALYSIS"/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz
        gunzip "$ANALYSIS"/REFERENCE/GRCh38.p13.genome.fa.gz
        gunzip "$ANALYSIS"/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz
        STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$ANALYSIS"/INDEX --genomeFastaFiles "$ANALYSIS"/REFERENCE/GRCh38.p13.genome.fa
        for i in $(find "$ANALYSIS"/TRIMMED -type f -name "*.fastq.gz" | awk -F/ '{print $NF}' | awk -F_ '{print $1"_"$2}' | sort -u); do
            STAR --runThreadN "$THREADS" --runMode alignReads --genomeDir "$ANALYSIS"/INDEX \
            --sjdbGTFfile "$ANALYSIS"/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf \
            --readFilesIn "$ANALYSIS"/TRIMMED/"$i"_R1_PAIRED.fastq.gz "$ANALYSIS"/TRIMMED/"$i"_R2_PAIRED.fastq.gz \
            --readFilesCommand zcat --outFileNamePrefix "$ANALYSIS"/ALIGN/"$i"_ --outSAMtype BAM Unsorted --outReadsUnmapped Fastx
        done
    fi

    if [[ "$LEVEL" == 6 ]]; then
        featureCounts -a "$ANALYSIS"/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf -o "$ANALYSIS"/counts_CHIKV_CTRL.txt \
        -s 2 -G "$ANALYSIS"/REFERENCE/GRCh38.p13.genome.fa -p --countReadPairs -C -T "$THREADS" \
        $(cat "$SAMPLESHEET" | awk -F, -v ALIGN="$ANALYSIS"/ALIGN/ '$2~/CHIKV|CTRL/ {print ALIGN $1"_Aligned.out.bam"}' | xargs)
        featureCounts -a "$ANALYSIS"/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf -o "$ANALYSIS"/counts_CHIKV_DENV.txt \
        -s 2 -G "$ANALYSIS"/REFERENCE/GRCh38.p13.genome.fa -p --countReadPairs -C -T "$THREADS" \
        $(cat "$SAMPLESHEET" | awk -F, -v ALIGN="$ANALYSIS"/ALIGN/ '$2~/CHIKV|DENV/ {print ALIGN $1"_Aligned.out.bam"}' | xargs)
        featureCounts -a "$ANALYSIS"/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf -o "$ANALYSIS"/counts_CHIKV_ZIKV.txt \
        -s 2 -G "$ANALYSIS"/REFERENCE/GRCh38.p13.genome.fa -p --countReadPairs -C -T "$THREADS" \
        $(cat "$SAMPLESHEET" | awk -F, -v ALIGN="$ANALYSIS"/ALIGN/ '$2~/CHIKV|ZIKV/ {print ALIGN $1"_Aligned.out.bam"}' | xargs)
        featureCounts -a "$ANALYSIS"/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf -o "$ANALYSIS"/counts_DENV_CTRL.txt \
        -s 2 -G "$ANALYSIS"/REFERENCE/GRCh38.p13.genome.fa -p --countReadPairs -C -T "$THREADS" \
        $(cat "$SAMPLESHEET" | awk -F, -v ALIGN="$ANALYSIS"/ALIGN/ '$2~/DENV|CTRL/ {print ALIGN $1"_Aligned.out.bam"}' | xargs)
        featureCounts -a "$ANALYSIS"/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf -o "$ANALYSIS"/counts_DENV_ZIKV.txt \
        -s 2 -G "$ANALYSIS"/REFERENCE/GRCh38.p13.genome.fa -p --countReadPairs -C -T "$THREADS" \
        $(cat "$SAMPLESHEET" | awk -F, -v ALIGN="$ANALYSIS"/ALIGN/ '$2~/DENV|ZIKV/ {print ALIGN $1"_Aligned.out.bam"}' | xargs)
        featureCounts -a "$ANALYSIS"/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf -o "$ANALYSIS"/counts_ZIKV_CTRL.txt \
        -s 2 -G "$ANALYSIS"/REFERENCE/GRCh38.p13.genome.fa -p --countReadPairs -C -T "$THREADS" \
        $(cat "$SAMPLESHEET" | awk -F, -v ALIGN="$ANALYSIS"/ALIGN/ '$2~/ZIKV|CTRL/ {print ALIGN $1"_Aligned.out.bam"}' | xargs)
    fi

    end=$(date +%s.%N)

    runtime=$(python -c "print(${end} - ${start})")

    echo "" && echo "Done. The runtime was $runtime seconds." && echo ""

}

bg $1 $2 $3 &>>rnaseq_arbovirus_zkproject_log.txt &
