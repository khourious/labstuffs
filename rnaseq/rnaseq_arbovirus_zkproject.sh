#!/bin/bash

# if [[ -z "$(which conda)" ]]; then
    # cd
    # wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    # bash Miniconda3-latest-Linux-x86_64.sh -bfp miniconda3
    # rm Miniconda3-latest-Linux-x86_64.sh
    # MYSHELL=$(echo $SHELL | awk -F/ '{print $NF}')
    # echo 'export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
    # export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH
    # conda install -y -c conda-forge mamba
    # mamba update -y -n base conda
    # mamba create -y -n rnaseq -c conda-forge -c anaconda -c bioconda -c defaults fastqc multiqc trimmomatic star subread
# else
    # if [[ -z "$(which mamba)" ]]; then
        # conda install -y -c conda-forge mamba
        # mamba update -y -n base conda
        # mamba create -y -n rnaseq -c conda-forge -c anaconda -c bioconda -c defaults fastqc multiqc trimmomatic star subread
    # else
        # mamba update -y -n base conda
        # mamba create -y -n rnaseq -c conda-forge -c anaconda -c bioconda -c defaults fastqc multiqc trimmomatic star subread
    # fi
# fi

bg() {

    start=$(date +%s.%N)

    source activate rnaseq

    THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

    RAWDIR=/mnt/x/rnaseq_arbovirus/ArbovirusFiocruzBA-83677594
    RAWRUNDIR1="$RAWDIR"/FASTQ_Generation_2018-08-17_15_33_27Z-116577544
    RAWRUNDIR2="$RAWDIR"/FASTQ_Generation_2018-10-17_13_51_20Z-130391988
    RAWRUNDIR3="$RAWDIR"/FASTQ_Generation_2018-12-11_13_17_10Z-141945716
    RAWRUNDIR4="$RAWDIR"/FASTQ_Generation_2018-12-13_14_45_41Z-143225097
    RAWRUNDIR5="$RAWDIR"/FASTQ_Generation_2018-12-15_10_26_29Z-143716575
    ANALYSIS="$RAWDIR"/ANALYSIS

    # mkdir "$ANALYSIS" "$ANALYSIS"/ALIGN "$ANALYSIS"/INDEX \
    # "$ANALYSIS"/QC_RUN{1..5} "$ANALYSIS"/QC_RUNS_MERGED \
    # "$ANALYSIS"/REFERENCE "$ANALYSIS"/RUN{1..5} \
    # "$ANALYSIS"/RUNS_MERGED "$ANALYSIS"/TRIMMED

    # find "$RAWRUNDIR1" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN1 {} +
    # find "$RAWRUNDIR2" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN2 {} +
    # find "$RAWRUNDIR3" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN3 {} +
    # find "$RAWRUNDIR4" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN4 {} +
    # find "$RAWRUNDIR5" -type f -name '*.fastq.gz' -exec cp -vat "$ANALYSIS"/RUN5 {} +

    # for i in $(find "$ANALYSIS" -maxdepth 1 -mindepth 1 -type d -name "RUN*[1-5]" | awk -F/ '{print $NF}'); do
        # fastqc -t "$THREADS" "$ANALYSIS"/"$i"/*.fastq.gz -o "$ANALYSIS"/QC_"$i"
        # multiqc -s -i "ArbovirusFiocruzBA-83677594 "$i" lanes" -ip --no-data-dir -n "$ANALYSIS"/"$i"_lanes_multiqc_report "$ANALYSIS"/QC_"$i"/*
    # done

    # for i in $(find "$ANALYSIS"/RUN{1..5} -type f -name "*.fastq.gz" | awk -F/ '{print $NF}' | awk -F_ '{print $1"_"$2}' | sort -u); do
        # cat "$ANALYSIS"/RUN1/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN1/"$i"_RUN1_R1.fastq.gz
        # cat "$ANALYSIS"/RUN1/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN1/"$i"_RUN1_R2.fastq.gz
        # cat "$ANALYSIS"/RUN2/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN2/"$i"_RUN2_R1.fastq.gz
        # cat "$ANALYSIS"/RUN2/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN2/"$i"_RUN2_R2.fastq.gz
        # cat "$ANALYSIS"/RUN3/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN3/"$i"_RUN3_R1.fastq.gz
        # cat "$ANALYSIS"/RUN3/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN3/"$i"_RUN3_R2.fastq.gz
        # cat "$ANALYSIS"/RUN4/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN4/"$i"_RUN4_R1.fastq.gz
        # cat "$ANALYSIS"/RUN4/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN4/"$i"_RUN4_R2.fastq.gz
        # cat "$ANALYSIS"/RUN5/"$i"_L00*_R1_001.fastq.gz > "$ANALYSIS"/RUN5/"$i"_RUN5_R1.fastq.gz
        # cat "$ANALYSIS"/RUN5/"$i"_L00*_R2_001.fastq.gz > "$ANALYSIS"/RUN5/"$i"_RUN5_R2.fastq.gz
    # done

    rm -rf "$ANALYSIS"/RUN{1..5}/*_001.fastq.gz

    for i in $(find "$ANALYSIS" -maxdepth 1 -mindepth 1 -type d -name "RUN*[1-5]" | awk -F/ '{print $NF}'); do
        fastqc -t "$THREADS" "$ANALYSIS"/"$i"/*RUN* -o "$ANALYSIS"/QC_"$i"
        multiqc -s -i "$i" -b "ArbovirusFiocruzBA-83677594 "$i" merged lanes" -ip --no-data-dir -n "$ANALYSIS"/"$i"_MERGED_LANES_multiqc_report "$ANALYSIS"/QC_"$i"/*RUN*
    done

    mv "$ANALYSIS"/RUN{1..5}/* "$ANALYSIS"/RUNS_MERGED

    rm -rf "$ANALYSIS"/RUN{1..5}

    # for i in $(find "$ANALYSIS"/RUNS_MERGED -type f -name "*.fastq.gz" | awk -F/ '{print $NF}' | awk -F_ '{print $1"_"$2}' | sort -u); do
        # cat "$ANALYSIS"/RUNS_MERGED/"$i"_RUN*_R1.fastq.gz > "$ANALYSIS"/RUNS_MERGED/"$i"_MERGED_R1.fastq.gz
        # cat "$ANALYSIS"/RUNS_MERGED/"$i"_RUN*_R2.fastq.gz > "$ANALYSIS"/RUNS_MERGED/"$i"_MERGED_R2.fastq.gz
    # done

    # rm -rf "$ANALYSIS"/RUNS_MERGED/*RUN*




    # for i in $(find "$RAWDIR"/ANALYSIS -maxdepth 1 -mindepth 1 -type d -name "RUNS_MERGED" | while read o; do basename $o; done); do
        # fastqc -t "$THREADS" "$RAWDIR"/ANALYSIS/"$i"/*.fastq.gz -o "$RAWDIR"/ANALYSIS/QC_"$i"
        # multiqc -s -i "ArbovirusFiocruzBA-83677594 "$i"" -ip --no-data-dir -n "$RAWDIR"/ANALYSIS/"$i"_multiqc_report "$RAWDIR"/ANALYSIS/QC_"$i"/*MERGED*
    # done

    # for i in $(find "$RAWDIR"/ANALYSIS/RUNS_MERGED -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1,2 | sort -u); do
        # trimmomatic PE -threads "$THREADS" -phred33 "$RAWDIR"/ANALYSIS/RUNS_MERGED/"$i"_MERGED_R1.fastq.gz "$RAWDIR"/ANALYSIS/RUNS_MERGED/"$i"_MERGED_R2.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R1_PAIRED.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R1_UNPAIRED.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R2_PAIRED.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"R2_UNPAIRED.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36
    # done

    # wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -q -O "$RAWDIR"/ANALYSIS/REFERENCE/GRCh38.p13.genome.fa.gz
    # wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz -q -O "$RAWDIR"/ANALYSIS/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz

    # gunzip "$RAWDIR"/ANALYSIS/REFERENCE/GRCh38.p13.genome.fa.gz
    # gunzip "$RAWDIR"/ANALYSIS/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz

    # STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$RAWDIR"/ANALYSIS/INDEX --genomeFastaFiles "$RAWDIR"/ANALYSIS/REFERENCE/GRCh38.p13.genome.fa

    # for i in $(find "$RAWDIR"/ANALYSIS/TRIMMED -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1,2 | sort -u); do
        # STAR --runThreadN "$THREADS" --runMode alignReads --genomeDir "$RAWDIR"/ANALYSIS/INDEX --sjdbGTFfile "$RAWDIR"/ANALYSIS/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf --readFilesIn "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R1_PAIRED.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R2_PAIRED.fastq.gz --readFilesCommand zcat --outFileNamePrefix "$RAWDIR"/ANALYSIS/ALIGN/"$i"_ --outSAMtype BAM Unsorted --outReadsUnmapped Fastx
    # done

    end=$(date +%s.%N)

    runtime=$(python -c "print(${end} - ${start})")

    echo "" && echo "Done. The runtime was $runtime seconds." && echo ""

}

bg &>>rnaseq_arbovirus_zkproject_log.txt &
