#!/bin/bash

bg() {

    start=$(date +%s.%N)

    source activate RNA-Seq

    RAWDIR="/home/laise/arboBA-RNAseq/ArbovirusFiocruzBA-83677594"
    RAWRUNDIR1="$RAWDIR"/FASTQ_Generation_2018-08-17_15_33_27Z-116577544
    RAWRUNDIR2="$RAWDIR"/FASTQ_Generation_2018-10-17_13_51_20Z-130391988
    RAWRUNDIR3="$RAWDIR"/FASTQ_Generation_2018-12-11_13_17_10Z-141945716
    RAWRUNDIR4="$RAWDIR"/FASTQ_Generation_2018-12-13_14_45_41Z-143225097
    RAWRUNDIR5="$RAWDIR"/FASTQ_Generation_2018-12-15_10_26_29Z-143716575
    THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

    mkdir "$RAWDIR"/ANALYSIS "$RAWDIR"/ANALYSIS/ALIGN "$RAWDIR"/ANALYSIS/INDEX "$RAWDIR"/ANALYSIS/QC_RUN{1..5} "$RAWDIR"/ANALYSIS/QC_RUNS_MERGED "$RAWDIR"/ANALYSIS/REFERENCE "$RAWDIR"/ANALYSIS/RUN{1..5} "$RAWDIR"/ANALYSIS/RUNS_MERGED "$RAWDIR"/ANALYSIS/TRIMMED

    find "$RAWRUNDIR1" -type f -name '*.fastq.gz' -exec cp -vat "$RAWDIR"/ANALYSIS/RUN1 {} +
    find "$RAWRUNDIR2" -type f -name '*.fastq.gz' -exec cp -vat "$RAWDIR"/ANALYSIS/RUN2 {} +
    find "$RAWRUNDIR3" -type f -name '*.fastq.gz' -exec cp -vat "$RAWDIR"/ANALYSIS/RUN3 {} +
    find "$RAWRUNDIR4" -type f -name '*.fastq.gz' -exec cp -vat "$RAWDIR"/ANALYSIS/RUN4 {} +
    find "$RAWRUNDIR5" -type f -name '*.fastq.gz' -exec cp -vat "$RAWDIR"/ANALYSIS/RUN5 {} +

    for i in $(find "$RAWDIR"/ANALYSIS -maxdepth 1 -mindepth 1 -type d -name "RUN*[1-5]" | while read o; do basename $o; done); do
        fastqc -t "$THREADS" "$RAWDIR"/ANALYSIS/"$i"/*.fastq.gz -o "$RAWDIR"/ANALYSIS/QC_"$i"
        multiqc -s -i "ArbovirusFiocruzBA-83677594 "$i" lanes" -ip --no-data-dir -n "$RAWDIR"/ANALYSIS/"$i"_lanes_multiqc_report "$RAWDIR"/ANALYSIS/QC_"$i"/*
    done

    for i in $(find "$RAWDIR"/ANALYSIS/RUN1 -type f -name "*.fastq.gz" | while read o; do basename $o | cut -d_ -f1,2; done | sort | uniq); do
        cat "$RAWDIR"/ANALYSIS/RUN1/"$i"_L00*_R1_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN1/"$i"_RUN1_R1.fastq.gz
        cat "$RAWDIR"/ANALYSIS/RUN1/"$i"_L00*_R2_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN1/"$i"_RUN1_R2.fastq.gz
    done

    for i in $(find "$RAWDIR"/ANALYSIS/RUN2 -type f -name "*.fastq.gz" | while read o; do basename $o | cut -d_ -f1,2; done | sort | uniq); do
        cat "$RAWDIR"/ANALYSIS/RUN2/"$i"_L00*_R1_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN2/"$i"_RUN2_R1.fastq.gz
        cat "$RAWDIR"/ANALYSIS/RUN2/"$i"_L00*_R2_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN2/"$i"_RUN2_R2.fastq.gz
    done

    for i in $(find "$RAWDIR"/ANALYSIS/RUN3 -type f -name "*.fastq.gz" | while read o; do basename $o | cut -d_ -f1,2; done | sort | uniq); do
        cat "$RAWDIR"/ANALYSIS/RUN3/"$i"_L00*_R1_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN3/"$i"_RUN3_R1.fastq.gz
        cat "$RAWDIR"/ANALYSIS/RUN3/"$i"_L00*_R2_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN3/"$i"_RUN3_R2.fastq.gz
    done

    for i in $(find "$RAWDIR"/ANALYSIS/RUN4 -type f -name "*.fastq.gz" | while read o; do basename $o | cut -d_ -f1,2; done | sort | uniq); do
        cat "$RAWDIR"/ANALYSIS/RUN4/"$i"_L00*_R1_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN4/"$i"_RUN4_R1.fastq.gz
        cat "$RAWDIR"/ANALYSIS/RUN4/"$i"_L00*_R2_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN4/"$i"_RUN4_R2.fastq.gz
    done

    for i in $(find "$RAWDIR"/ANALYSIS/RUN5 -type f -name "*.fastq.gz" | while read o; do basename $o | cut -d_ -f1,2; done | sort | uniq); do
        cat "$RAWDIR"/ANALYSIS/RUN5/"$i"_L00*_R1_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN5/"$i"_RUN5_R1.fastq.gz
        cat "$RAWDIR"/ANALYSIS/RUN5/"$i"_L00*_R2_001.fastq.gz > "$RAWDIR"/ANALYSIS/RUN5/"$i"_RUN5_R2.fastq.gz
    done

    rm -rf "$RAWDIR"/ANALYSIS/RUN{1..5}/*_001.fastq.gz

    for i in $(find "$RAWDIR"/ANALYSIS -maxdepth 1 -mindepth 1 -type d -name "RUN*" | while read o; do basename $o; done); do
        fastqc -t "$THREADS" "$RAWDIR"/ANALYSIS/"$i"/*RUN* -o "$RAWDIR"/ANALYSIS/QC_"$i"
        multiqc -s -i "$i" -b "ArbovirusFiocruzBA-83677594 "$i" merged lanes" -ip --no-data-dir -n "$RAWDIR"/ANALYSIS/"$i"_lanes_MERGED_multiqc_report "$RAWDIR"/ANALYSIS/QC_"$i"/*RUN*
    done

    mv "$RAWDIR"/ANALYSIS/RUN{1..5}/* "$RAWDIR"/ANALYSIS/RUNS_MERGED

    rm -rf "$RAWDIR"/ANALYSIS/RUN{1..5}

    for i in $(find "$RAWDIR"/ANALYSIS/RUNS_MERGED -type f -name "*.fastq.gz" | while read o; do basename $o | cut -d_ -f1,2; done | sort | uniq); do
        cat "$RAWDIR"/ANALYSIS/RUNS_MERGED/"$i"_RUN*_R1.fastq.gz > "$RAWDIR"/ANALYSIS/RUNS_MERGED/"$i"_MERGED_R1.fastq.gz
        cat "$RAWDIR"/ANALYSIS/RUNS_MERGED/"$i"_RUN*_R2.fastq.gz > "$RAWDIR"/ANALYSIS/RUNS_MERGED/"$i"_MERGED_R2.fastq.gz
    done

    rm -rf "$RAWDIR"/ANALYSIS/RUNS_MERGED/*RUN*

    for i in $(find "$RAWDIR"/ANALYSIS -maxdepth 1 -mindepth 1 -type d -name "RUNS_MERGED" | while read o; do basename $o; done); do
        fastqc -t "$THREADS" "$RAWDIR"/ANALYSIS/"$i"/*.fastq.gz -o "$RAWDIR"/ANALYSIS/QC_"$i"
        multiqc -s -i "ArbovirusFiocruzBA-83677594 "$i"" -ip --no-data-dir -n "$RAWDIR"/ANALYSIS/"$i"_multiqc_report "$RAWDIR"/ANALYSIS/QC_"$i"/*MERGED*
    done

    for i in $(find "$RAWDIR"/ANALYSIS/RUNS_MERGED -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1,2 | sort -u); do
        trimmomatic PE -threads "$THREADS" -phred33 "$RAWDIR"/ANALYSIS/RUNS_MERGED/"$i"_MERGED_R1.fastq.gz "$RAWDIR"/ANALYSIS/RUNS_MERGED/"$i"_MERGED_R2.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R1_PAIRED.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R1_UNPAIRED.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R2_PAIRED.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"R2_UNPAIRED.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36
    done

    wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -q -O "$RAWDIR"/ANALYSIS/REFERENCE/GRCh38.p13.genome.fa.gz
    wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz -q -O "$RAWDIR"/ANALYSIS/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz

    gunzip "$RAWDIR"/ANALYSIS/REFERENCE/GRCh38.p13.genome.fa.gz
    gunzip "$RAWDIR"/ANALYSIS/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz

    STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$RAWDIR"/ANALYSIS/INDEX --genomeFastaFiles "$RAWDIR"/ANALYSIS/REFERENCE/GRCh38.p13.genome.fa

    for i in $(find "$RAWDIR"/ANALYSIS/TRIMMED -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1,2 | sort -u); do
        STAR --runThreadN "$THREADS" --runMode alignReads --genomeDir "$RAWDIR"/ANALYSIS/INDEX --sjdbGTFfile "$RAWDIR"/ANALYSIS/REFERENCE/gencode.v38.chr_patch_hapl_scaff.annotation.gtf --readFilesIn "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R1_PAIRED.fastq.gz "$RAWDIR"/ANALYSIS/TRIMMED/"$i"_R2_PAIRED.fastq.gz --readFilesCommand zcat --outFileNamePrefix "$RAWDIR"/ANALYSIS/ALIGN/"$i"_ --outSAMtype BAM Unsorted --outReadsUnmapped Fastx
    done

    end=$(date +%s.%N)

    runtime=$(python -c "print(${end} - ${start})")

    echo "" && echo "Done. The runtime was $runtime seconds."

}

bg &>rnaseq_arbovirus_log_"$(date +'%Y-%m-%d')".txt &
