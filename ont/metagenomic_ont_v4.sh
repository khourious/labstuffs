#!/bin/bash

# Set the path of RAW fast5 files
RAWDIR=/mnt/c/Dropbox/kalabric/NGS_LIBRARY01_20210408

# Set the path of kraken2 database
KRAKEN2DB=/mnt/x/kraken2db/minikraken2_v2_8GB_201904_UPDATE

# Set the path of reference sequences
REFSEQS=/mnt/x/kalabric/REFSEQS

[ ! -d "$REFSEQS" ] && mkdir -p "$REFSEQS"

# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -q -O "$REFSEQS"/GRCh38.p13.genome.fa.gz

# source activate ont_metagenomic

# esearch -db nucleotide -query "NC_004162.2" | efetch -format fasta > "$REFSEQS"/CHIKV_NC_004162.2.fa
# esearch -db nucleotide -query "NC_001477.1" | efetch -format fasta > "$REFSEQS"/DENV1_NC_001477.1.fa
# esearch -db nucleotide -query "NC_001474.2" | efetch -format fasta > "$REFSEQS"/DENV2_NC_001474.2.fa
# esearch -db nucleotide -query "NC_001475.2" | efetch -format fasta > "$REFSEQS"/DENV3_NC_001475.2.fa
# esearch -db nucleotide -query "NC_002640.1" | efetch -format fasta > "$REFSEQS"/DENV4_NC_002640.1.fa
# esearch -db nucleotide -query "NC_035889.1" | efetch -format fasta > "$REFSEQS"/ZIKV_NC_035889.1.fa

# kraken2 db: viral
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20210517.tar.gz
# tar -vzxf k2_viral_20210517.tar.gz

# kraken2 db: archaea, human, plasmid, UniVec_Core, viral
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20210517.tar.gz
# tar -vzxf k2_minusb_20210517.tar.gz

# kraken2 db: archaea, bacteria, human, plasmid, UniVec_Core, viral
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz
# tar -vzxf k2_standard_8gb_20210517.tar.gz

# kraken2 db: archaea, bacteria, fungi, human, plasmid, protozoa, UniVec_Core, viral
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_8gb_20210517.tar.gz
# tar -vzxf k2_pluspf_8gb_20210517.tar.gz

# kraken2 db: archaea, bacteria, fungi, human, plant, plasmid, protozoa, UniVec_Core, viral
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_8gb_20210517.tar.gz
# tar -vzxf k2_pluspfp_8gb_20210517.tar.gz

MYSHELL=$(echo $SHELL | awk -F/ '{print $NF}')

# if [[ -z "$(which guppy_basecaller)" ]]; then
    # version=5.0.11
    # cd
    # curl https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_"$version"_linux64.tar.gz -o ont-guppy.tar.gz
    # tar -vzxf ont-guppy.tar.gz
    # rm -rf ont-guppy.tar.gz
    # echo 'export PATH=$HOME/ont-guppy/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
    # export PATH=$HOME/ont-guppy/bin:/usr/local/share/rsi/idl/bin:$PATH
# else
    # guppy_basecaller --version
# fi

# if [[ -z "$(which conda)" ]]; then
    # cd
    # wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    # bash Miniconda3-latest-Linux-x86_64.sh -bfp miniconda3
    # rm Miniconda3-latest-Linux-x86_64.sh
    # echo 'export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
    # export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH
    # conda install -y -c conda-forge mamba
    # mamba update -y -n base conda
    # mamba create -y -n ont_metagenomic -c conda-forge -c anaconda -c bioconda -c defaults entrez-direct ivar minimap2 nanopolish racon samtools
    # mamba create -y -n ont_qc -c aleg -c conda-forge -c anaconda -c bioconda -c defaults python=3.6 pycoqc
    # mamba create -y -n plot -c conda-forge -c anaconda -c bioconda -c defaults pysam numpy pandas seaborn
# else
    # if [[ -z "$(which mamba)" ]]; then
        # conda install -y -c conda-forge mamba
        # mamba update -y -n base conda
        # mamba create -y -n ont_metagenomic -c conda-forge -c anaconda -c bioconda -c defaults entrez-direct ivar minimap2 nanopolish racon samtools
        # mamba create -y -n ont_qc -c aleg -c conda-forge -c anaconda -c bioconda -c defaults python=3.6 pycoqc
        # mamba create -y -n plot -c conda-forge -c anaconda -c bioconda -c defaults pysam numpy pandas seaborn
    # else
        # mamba update -y -n base conda
        # mamba create -y -n ont_metagenomic -c conda-forge -c anaconda -c bioconda -c defaults entrez-direct ivar minimap2 nanopolish racon samtools
        # mamba create -y -n ont_qc -c aleg -c conda-forge -c anaconda -c bioconda -c defaults python=3.6 pycoqc
        # mamba create -y -n plot -c conda-forge -c anaconda -c bioconda -c defaults pysam numpy pandas seaborn
    # fi
# fi

# if [[ -z "$(which kraken2)" ]]; then
    # version=2.1.2
    # cd
    # wget https://github.com/DerrickWood/kraken2/archive/v$version.tar.gz
    # tar -vzxf v$version.tar.gz
    # rm -rf v$version.tar.gz
    # cd kraken2-$version
    # ./install_kraken2.sh $HOME/kraken2-$version
    # [ ! -d "$HOME/bin" ] && mkdir -p $HOME/bin
    # cp $HOME/kraken2-$version/kraken2{,-build,-inspect} $HOME/bin
# else
    # kraken2 --help
# fi

# if [[ -z "$(which fastcov.py)" ]]; then
    # cd
    # git clone https://github.com/RaverJay/fastcov
    # cd fastcov
    # echo 'export PATH=$HOME/fastcov:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
    # export PATH=$HOME/fastcov:/usr/local/share/rsi/idl/bin:$PATH
# else
    # source activate plot
    # fastcov.py --help
# fi

bg() {

    start=$(date +%s.%N)

    THREADS=$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')

    RUNNAME=$(basename "$RAWDIR")
    BASECALLDIR="$RAWDIR"/BASECALL
    DEMUXDIR="$RAWDIR"/DEMUX
    DEMUXCATDIR="$RAWDIR"/DEMUX_CAT
    READLEVELDIR="$RAWDIR"/LEVEL_READS
    CONTIGLEVELDIR="$RAWDIR"/LEVEL_CONTIGS
    ASSEMBLYDIR="$RAWDIR"/ASSEMBLY

    [ ! -d "$DEMUXCATDIR" ] && mkdir -vp "$DEMUXCATDIR"
    [ ! -d "$READLEVELDIR" ] && mkdir -vp "$READLEVELDIR"
    [ ! -d "$ASSEMBLYDIR" ] && mkdir -vp "$ASSEMBLYDIR"

    HUMANREFSEQ="$RAWDIR"/REFSEQS/GRCh38.p13.genome.fa.gz
    CHIKVREFSEQ="$RAWDIR"/REFSEQS/CHIKV_NC_004162.2.fa
    DENV1REFSEQ="$RAWDIR"/REFSEQS/DENV1_NC_001477.1.fa
    DENV2REFSEQ="$RAWDIR"/REFSEQS/DENV2_NC_001474.2.fa
    DENV3REFSEQ="$RAWDIR"/REFSEQS/DENV3_NC_001475.2.fa
    DENV4REFSEQ="$RAWDIR"/REFSEQS/DENV4_NC_001475.2.fa
    ZIKVREFSEQ="$RAWDIR"/REFSEQS/ZIKV_NC_035889.1.fa

    CONFIG=dna_r9.4.1_450bps_sup.cfg #dna_r9.4.1_450bps_hac.cfg
    ARRANGEMENTS="barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
    TRIMADAPTER=18

    guppy_basecaller -r -x auto --verbose_logs -c "$CONFIG" -i "$RAWDIR" -s "$BASECALLDIR" -q 1 --min_qscore 10 \
    --chunk_size 1000 --num_callers "$THREADS" --gpu_runners_per_device 1 --disable_pings

    guppy_barcoder -r --require_barcodes_both_ends --trim_barcodes -t "$THREADS" \
    -i "$BASECALLDIR" -s "$DEMUXDIR" --arrangements_files "$ARRANGEMENTS" \
    --num_extra_bases_trim "$TRIMADAPTER" -x auto

    source activate ont_qc

    pycoQC -q -f "$BASECALLDIR"/sequencing_summary.txt -b "$DEMUXDIR"/barcoding_summary.txt -o "$RAWDIR"/"$RUNNAME"_qc.html --report_title "$RUNNAME"

    for i in $(find "$DEMUXDIR" -mindepth 1 -type d -name "barcode*" -exec basename {} \; | sort); do
        [ -d "$DEMUXDIR"/"$i" ] && cat "$DEMUXDIR"/"$i"/*.fastq > "$DEMUXCATDIR"/"$i".fastq
    done

    source activate ont_metagenomic

    for i in $(find "$DEMUXCATDIR" -type f -name "*.fastq" -exec basename {} \; | awk -F. '{print $1}' | sort -u); do
        minimap2 -ax map-ont -t "$THREADS" "$HUMANREFSEQ" "$DEMUXCATDIR"/"$i".fastq | samtools sort -@ "$THREADS" -o "$READLEVELDIR"/"$i".sorted.bam -
        samtools index -@ "$THREADS" "$READLEVELDIR"/"$i".sorted.bam
        samtools view -@ "$THREADS" -bS -f 4 "$READLEVELDIR"/"$i".sorted.bam > "$READLEVELDIR"/"$i".sorted.filtered.sam
        samtools fastq -@ "$THREADS" -f 4 "$READLEVELDIR"/"$i".sorted.filtered.sam > "$READLEVELDIR"/"$i".sorted.filtered.fastq
        minimap2 -ax ava-ont -t "$THREADS" "$READLEVELDIR"/"$i".sorted.filtered.fastq "$READLEVELDIR"/"$i".sorted.filtered.fastq > "$READLEVELDIR"/"$i".sorted.filtered.overlap.sam
        racon -t "$THREADS" -f -u "$READLEVELDIR"/"$i".sorted.filtered.fastq "$READLEVELDIR"/"$i".sorted.filtered.overlap.sam "$READLEVELDIR"/"$i".sorted.filtered.fastq > "$READLEVELDIR"/"$i".sorted.filtered.corrected.fasta
    done

    for i in $(find "$READLEVELDIR" -type f -name "*.fasta" -exec basename {} \; | awk -F. '{print $1}' | sort -u); do
        kraken2 --db "$KRAKEN2DB" --threads "$THREADS" --report "$READLEVELDIR"/"$i"_report.txt --report-minimizer-data --output "$READLEVELDIR"/"$i"_output.txt "$READLEVELDIR"/"$i".sorted.filtered.corrected.fasta
        minimap2 -t "$THREADS" -ax map-ont "$CHIKVREFSEQ" "$READLEVELDIR"/"$i".sorted.filtered.corrected.fasta | samtools sort -@ "$THREADS" -o "$ASSEMBLYDIR"/"$i".chikv.sorted.bam -
        samtools view -@ "$THREADS" -h -F 4 -b "$ASSEMBLYDIR"/"$i".chikv.sorted.bam > "$ASSEMBLYDIR"/"$i".chikv.sorted.mapped.bam
        samtools index -@ "$THREADS" "$ASSEMBLYDIR"/"$i".chikv.sorted.mapped.bam
        samtools mpileup -A -B -Q 0 --reference "$CHIKVREFSEQ" "$ASSEMBLYDIR"/"$i".chikv.sorted.mapped.bam | ivar consensus -p "$ASSEMBLYDIR"/"$i".chikv -n N -i "$i"
        minimap2 -t "$THREADS" -ax map-ont "$DENV1REFSEQ" "$READLEVELDIR"/"$i".corrected.fasta | samtools sort -@ "$THREADS" -o "$ASSEMBLYDIR"/"$i".denv1.sorted.bam -
        samtools view -@ "$THREADS" -h -F 4 -b "$ASSEMBLYDIR"/"$i".denv1.sorted.bam > "$ASSEMBLYDIR"/"$i".denv1.sorted.mapped.bam
        samtools index -@ "$THREADS" "$ASSEMBLYDIR"/"$i".denv1.sorted.mapped.bam
        samtools mpileup -A -B -Q 0 --reference "$DENV1REFSEQ" "$ASSEMBLYDIR"/"$i".denv1.sorted.mapped.bam | ivar consensus -p "$ASSEMBLYDIR"/"$i".denv1 -n N -i "$i"
        minimap2 -t "$THREADS" -ax map-ont "$DENV2REFSEQ" "$READLEVELDIR"/"$i".corrected.fasta | samtools sort -@ "$THREADS" -o "$ASSEMBLYDIR"/"$i".denv2.sorted.bam -
        samtools view -@ "$THREADS" -h -F 4 -b "$ASSEMBLYDIR"/"$i".denv2.sorted.bam > "$ASSEMBLYDIR"/"$i".denv2.sorted.mapped.bam
        samtools index -@ "$THREADS" "$ASSEMBLYDIR"/"$i".denv2.sorted.mapped.bam
        samtools mpileup -A -B -Q 0 --reference "$DENV2REFSEQ" "$ASSEMBLYDIR"/"$i".denv2.sorted.mapped.bam | ivar consensus -p "$ASSEMBLYDIR"/"$i".denv2 -n N -i "$i"
        minimap2 -t "$THREADS" -ax map-ont "$DENV3REFSEQ" "$READLEVELDIR"/"$i".corrected.fasta | samtools sort -@ "$THREADS" -o "$ASSEMBLYDIR"/"$i".denv3.sorted.bam -
        samtools view -@ "$THREADS" -h -F 4 -b "$ASSEMBLYDIR"/"$i".denv3.sorted.bam > "$ASSEMBLYDIR"/"$i".denv3.sorted.mapped.bam
        samtools index -@ "$THREADS" "$ASSEMBLYDIR"/"$i".denv3.sorted.mapped.bam
        samtools mpileup -A -B -Q 0 --reference "$DENV3REFSEQ" "$ASSEMBLYDIR"/"$i".denv3.sorted.mapped.bam | ivar consensus -p "$ASSEMBLYDIR"/"$i".denv3 -n N -i "$i"
        minimap2 -t "$THREADS" -ax map-ont "$DENV4REFSEQ" "$READLEVELDIR"/"$i".corrected.fasta | samtools sort -@ "$THREADS" -o "$ASSEMBLYDIR"/"$i".denv4.sorted.bam -
        samtools view -@ "$THREADS" -h -F 4 -b "$ASSEMBLYDIR"/"$i".denv4.sorted.bam > "$ASSEMBLYDIR"/"$i".denv4.sorted.mapped.bam
        samtools index -@ "$THREADS" "$ASSEMBLYDIR"/"$i".denv4.sorted.mapped.bam
        samtools mpileup -A -B -Q 0 --reference "$DENV4REFSEQ" "$ASSEMBLYDIR"/"$i".denv4.sorted.mapped.bam | ivar consensus -p "$ASSEMBLYDIR"/"$i".denv4 -n N -i "$i"
        minimap2 -t "$THREADS" -ax map-ont "$ZIKVREFSEQ" "$READLEVELDIR"/"$i".corrected.fasta | samtools sort -@ "$THREADS" -o "$ASSEMBLYDIR"/"$i".zikv.sorted.bam -
        samtools view -@ "$THREADS" -h -F 4 -b "$ASSEMBLYDIR"/"$i".zikv.sorted.bam > "$ASSEMBLYDIR"/"$i".zikv.sorted.mapped.bam
        samtools index -@ "$THREADS" "$ASSEMBLYDIR"/"$i".zikv.sorted.mapped.bam
        samtools mpileup -A -B -Q 0 --reference "$ZIKVREFSEQ" "$ASSEMBLYDIR"/"$i".zikv.sorted.mapped.bam | ivar consensus -p "$ASSEMBLYDIR"/"$i".zikv -n N -i "$i"
    done

    source activate plot
	
    for i in $(find "$ASSEMBLYDIR" -type f -name "*.sorted.mapped.bam" -exec basename {} \; | awk -F. '{print $1}' | sort -u); do
        fastcov.py -l "$ASSEMBLYDIR"/barcode*.chikv.sorted.mapped.bam -o "$ASSEMBLYDIR"/chikv.coverage.pdf
        fastcov.py -l "$ASSEMBLYDIR"/barcode*.denv1.sorted.mapped.bam -o "$ASSEMBLYDIR"/denv1.coverage.pdf
        fastcov.py -l "$ASSEMBLYDIR"/barcode*.denv2.sorted.mapped.bam -o "$ASSEMBLYDIR"/denv2.coverage.pdf
        fastcov.py -l "$ASSEMBLYDIR"/barcode*.denv3.sorted.mapped.bam -o "$ASSEMBLYDIR"/denv3.coverage.pdf
        fastcov.py -l "$ASSEMBLYDIR"/barcode*.denv4.sorted.mapped.bam -o "$ASSEMBLYDIR"/denv4.coverage.pdf
        fastcov.py -l "$ASSEMBLYDIR"/barcode*.zikv.sorted.mapped.bam -o "$ASSEMBLYDIR"/zikv.coverage.pdf
    done
    
    end=$(date +%s.%N)

    runtime=$(python -c "print(${end} - ${start})")

    echo "" && echo "Done. The runtime was $runtime seconds." && echo ""

}

bg &>>metagenomic_ont_v4_log.txt &
