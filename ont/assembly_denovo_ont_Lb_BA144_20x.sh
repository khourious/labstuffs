#!/bin/bash

if [[ -z "$(which conda)" ]]; then
    cd
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -bfp miniconda3
    rm Miniconda3-latest-Linux-x86_64.sh
    MYSHELL=$(echo $SHELL | awk -F/ '{print $NF}')
    echo 'export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
    export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH
    conda install -y -c conda-forge mamba
    mamba update -y -n base conda
    mamba create -y -n flye -c conda-forge -c anaconda -c bioconda -c defaults
    mamba create -y -n medaka -c conda-forge -c anaconda -c bioconda -c defaults medaka bcftools minimap2 samtools
    mamba create -y -n nanopolish -c conda-forge -c anaconda -c bioconda -c defaults nanopolish samtools
    mamba create -y -n racon -c conda-forge -c anaconda -c bioconda -c defaults nanopolish racon
    mamba create -y -n ragtag -c conda-forge -c anaconda -c bioconda -c defaults
else
    if [[ -z "$(which mamba)" ]]; then
        conda install -y -c conda-forge mamba
        mamba update -y -n base conda
        mamba create -y -n flye -c conda-forge -c anaconda -c bioconda -c defaults
        mamba create -y -n medaka -c conda-forge -c anaconda -c bioconda -c defaults medaka bcftools minimap2 samtools
        mamba create -y -n nanopolish -c conda-forge -c anaconda -c bioconda -c defaults nanopolish samtools
        mamba create -y -n racon -c conda-forge -c anaconda -c bioconda -c defaults nanopolish racon
        mamba create -y -n ragtag -c conda-forge -c anaconda -c bioconda -c defaults
    else
        mamba update -y -n base conda
        mamba create -y -n flye -c conda-forge -c anaconda -c bioconda -c defaults
        mamba create -y -n medaka -c conda-forge -c anaconda -c bioconda -c defaults medaka bcftools minimap2 samtools
        mamba create -y -n nanopolish -c conda-forge -c anaconda -c bioconda -c defaults nanopolish samtools
        mamba create -y -n racon -c conda-forge -c anaconda -c bioconda -c defaults nanopolish racon
        mamba create -y -n ragtag -c conda-forge -c anaconda -c bioconda -c defaults
    fi
fi

THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

RAWDIR="/media/khourilab/Data/Lb_BA144/library15-2019-11-15"
HACDEMUXDIR="/media/khourilab/Data/Lb_BA144/Lb_BA144_HAC_Demux"
OUTPUTDIR="/media/khourilab/Data/Lb_BA144/assembly_denovo_20x"

SAMPLE="/media/khourilab/Data/Lb_BA144/Lb_BA144.fastq"
SAMPLEID="Lb_BA144"
REFSEQ="/media/khourilab/Data/Lb_BA144/TriTrypDB/TriTrypDB-50_LbraziliensisMHOMBR75M2904_2019_Genome.fasta"
PLOIDY="2"

mkdir flye
source activate flye
flye -t $THREADS --nano-raw $SAMPLE -g 33m -o $OUTPUTDIR/flye --asm-coverage 40
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye/ragtag $REFSEQ $OUTPUTDIR/flye/assembly.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye/ragtag $REFSEQ $OUTPUTDIR/flye/ragtag/assembly.corrected.fasta

mkdir flye.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye/assembly.fasta -o $OUTPUTDIR/flye.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.medaka/ragtag/consensus.corrected.fasta

mkdir flye.medaka.nanopolish
mkdir flye.medaka.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.medaka/consensus.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.medaka.nanopolish/consensus.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.medaka.nanopolish/consensus.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.medaka/consensus.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.medaka.nanopolish/consensus.sorted.bam -o $OUTPUTDIR/flye.medaka.nanopolish/vcf/consensus.{#}.vcf -g $OUTPUTDIR/flye.medaka/consensus.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.medaka.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.medaka.nanopolish/consensus.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.medaka/consensus.fasta -f $OUTPUTDIR/flye.medaka.nanopolish/consensus.vcflist.txt > $OUTPUTDIR/flye.medaka.nanopolish/consensus.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish/consensus.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish/ragtag/consensus.nanopolish.corrected.fasta

mkdir flye.medaka.nanopolish.racon1
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.medaka.nanopolish/consensus.nanopolish.fasta $OUTPUTDIR/flye.medaka.nanopolish/consensus.nanopolish.fasta > $OUTPUTDIR/flye.medaka.nanopolish.racon1/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.medaka.nanopolish/consensus.nanopolish.fasta $OUTPUTDIR/flye.medaka.nanopolish.racon1/assembly.sam $OUTPUTDIR/flye.medaka.nanopolish/consensus.nanopolish.fasta > $OUTPUTDIR/flye.medaka.nanopolish.racon1/assembly.racon1.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.nanopolish.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish.racon1/assembly.racon1.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.nanopolish.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish.racon1/ragtag/assembly.racon1.corrected.fasta

mkdir flye.medaka.nanopolish.racon2
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.medaka.nanopolish.racon1/assembly.racon1.fasta $OUTPUTDIR/flye.medaka.nanopolish.racon1/assembly.racon1.fasta > $OUTPUTDIR/flye.medaka.nanopolish.racon2/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.medaka.nanopolish.racon1/assembly.racon1.fasta $OUTPUTDIR/flye.medaka.nanopolish.racon2/assembly.sam $OUTPUTDIR/flye.medaka.nanopolish.racon1/assembly.racon1.fasta > $OUTPUTDIR/flye.medaka.nanopolish.racon2/assembly.racon2.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.nanopolish.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish.racon2/assembly.racon2.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.nanopolish.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish.racon2/ragtag/assembly.racon2.corrected.fasta

mkdir flye.medaka.nanopolish.racon3
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.medaka.nanopolish.racon2/assembly.racon2.fasta $OUTPUTDIR/flye.medaka.nanopolish.racon2/assembly.racon2.fasta > $OUTPUTDIR/flye.medaka.nanopolish.racon3/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.medaka.nanopolish.racon2/assembly.racon2.fasta $OUTPUTDIR/flye.medaka.nanopolish.racon3/assembly.sam $OUTPUTDIR/flye.medaka.nanopolish.racon2/assembly.racon2.fasta > $OUTPUTDIR/flye.medaka.nanopolish.racon3/assembly.racon3.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.nanopolish.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish.racon3/assembly.racon3.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.nanopolish.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish.racon3/ragtag/assembly.racon3.corrected.fasta

mkdir flye.medaka.nanopolish.racon4
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.medaka.nanopolish.racon3/assembly.racon3.fasta $OUTPUTDIR/flye.medaka.nanopolish.racon3/assembly.racon3.fasta > $OUTPUTDIR/flye.medaka.nanopolish.racon4/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.medaka.nanopolish.racon3/assembly.racon3.fasta $OUTPUTDIR/flye.medaka.nanopolish.racon4/assembly.sam $OUTPUTDIR/flye.medaka.nanopolish.racon3/assembly.racon3.fasta > $OUTPUTDIR/flye.medaka.nanopolish.racon4/assembly.racon4.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.nanopolish.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish.racon4/assembly.racon4.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.nanopolish.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.nanopolish.racon4/ragtag/assembly.racon4.corrected.fasta

mkdir flye.medaka.racon1
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.medaka/consensus.fasta $OUTPUTDIR/flye.medaka/consensus.fasta > $OUTPUTDIR/flye.medaka.racon1/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.medaka/consensus.fasta $OUTPUTDIR/flye.medaka.racon1/assembly.sam $OUTPUTDIR/flye.medaka/consensus.fasta > $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon1/ragtag/assembly.racon1.corrected.fasta

mkdir flye.medaka.racon2
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta > $OUTPUTDIR/flye.medaka.racon2/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta $OUTPUTDIR/flye.medaka.racon2/assembly.sam $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta > $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon2/ragtag/assembly.racon2.corrected.fasta

mkdir flye.medaka.racon3
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta > $OUTPUTDIR/flye.medaka.racon3/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta $OUTPUTDIR/flye.medaka.racon3/assembly.sam $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta > $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon3/ragtag/assembly.racon3.corrected.fasta

mkdir flye.medaka.racon4
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta > $OUTPUTDIR/flye.medaka.racon4/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta $OUTPUTDIR/flye.medaka.racon4/assembly.sam $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta > $OUTPUTDIR/flye.medaka.racon4/assembly.racon4.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon4/assembly.racon4.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon4/ragtag/assembly.racon4.corrected.fasta

mkdir flye.medaka.racon1.nanopolish
mkdir flye.medaka.racon1.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.medaka.racon1.nanopolish/assembly.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.medaka.racon1.nanopolish/assembly.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.medaka.racon1.nanopolish/assembly.sorted.bam -o $OUTPUTDIR/flye.medaka.racon1.nanopolish/vcf/assembly.{#}.vcf -g $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.medaka.racon1.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.medaka.racon1.nanopolish/assembly.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.medaka.racon1/assembly.racon1.fasta -f $OUTPUTDIR/flye.medaka.racon1.nanopolish/assembly.vcflist.txt > $OUTPUTDIR/flye.medaka.racon1.nanopolish/assembly.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.racon1.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon1.nanopolish/assembly.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.racon1.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon1.nanopolish/ragtag/assembly.nanopolish.corrected.fasta

mkdir flye.medaka.racon2.nanopolish
mkdir flye.medaka.racon2.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.medaka.racon2.nanopolish/assembly.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.medaka.racon2.nanopolish/assembly.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.medaka.racon2.nanopolish/assembly.sorted.bam -o $OUTPUTDIR/flye.medaka.racon2.nanopolish/vcf/assembly.{#}.vcf -g $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.medaka.racon2.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.medaka.racon2.nanopolish/assembly.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.medaka.racon2/assembly.racon2.fasta -f $OUTPUTDIR/flye.medaka.racon2.nanopolish/assembly.vcflist.txt > $OUTPUTDIR/flye.medaka.racon2.nanopolish/assembly.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.racon2.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon2.nanopolish/assembly.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.racon2.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon2.nanopolish/ragtag/assembly.nanopolish.corrected.fasta

mkdir flye.medaka.racon3.nanopolish
mkdir flye.medaka.racon3.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.medaka.racon3.nanopolish/assembly.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.medaka.racon3.nanopolish/assembly.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.medaka.racon3.nanopolish/assembly.sorted.bam -o $OUTPUTDIR/flye.medaka.racon3.nanopolish/vcf/assembly.{#}.vcf -g $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.medaka.racon3.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.medaka.racon3.nanopolish/assembly.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.medaka.racon3/assembly.racon3.fasta -f $OUTPUTDIR/flye.medaka.racon3.nanopolish/assembly.vcflist.txt > $OUTPUTDIR/flye.medaka.racon3.nanopolish/assembly.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.racon3.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon3.nanopolish/assembly.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.racon3.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon3.nanopolish/ragtag/assembly.nanopolish.corrected.fasta

mkdir flye.medaka.racon4.nanopolish
mkdir flye.medaka.racon4.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.medaka.racon4/assembly.racon4.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.medaka.racon4.nanopolish/assembly.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.medaka.racon4.nanopolish/assembly.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.medaka.racon4/assembly.racon4.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.medaka.racon4.nanopolish/assembly.sorted.bam -o $OUTPUTDIR/flye.medaka.racon4.nanopolish/vcf/assembly.{#}.vcf -g $OUTPUTDIR/flye.medaka.racon4/assembly.racon4.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.medaka.racon4.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.medaka.racon4.nanopolish/assembly.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.medaka.racon4/assembly.racon4.fasta -f $OUTPUTDIR/flye.medaka.racon4.nanopolish/assembly.vcflist.txt > $OUTPUTDIR/flye.medaka.racon4.nanopolish/assembly.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.medaka.racon4.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon4.nanopolish/assembly.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.medaka.racon4.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.medaka.racon4.nanopolish/ragtag/assembly.nanopolish.corrected.fasta

mkdir flye.nanopolish
mkdir flye.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye/assembly.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.nanopolish/assembly.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.nanopolish/assembly.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye/assembly.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.nanopolish/assembly.sorted.bam -o $OUTPUTDIR/flye.nanopolish/vcf/assembly.{#}.vcf -g $OUTPUTDIR/flye/assembly.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.nanopolish/assembly.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye/assembly.fasta -f $OUTPUTDIR/flye.nanopolish/assembly.vcflist.txt > $OUTPUTDIR/flye.nanopolish/assembly.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish/assembly.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish/ragtag/assembly.nanopolish.corrected.fasta

mkdir flye.nanopolish.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.nanopolish/assembly.nanopolish.fasta -o $OUTPUTDIR/flye.nanopolish.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka/ragtag/consensus.corrected.fasta

mkdir flye.nanopolish.medaka.racon1
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.nanopolish.medaka/consensus.fasta $OUTPUTDIR/flye.nanopolish.medaka/consensus.fasta > $OUTPUTDIR/flye.nanopolish.medaka.racon1/consensus.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.nanopolish.medaka/consensus.fasta $OUTPUTDIR/flye.nanopolish.medaka.racon1/consensus.sam $OUTPUTDIR/flye.nanopolish.medaka/consensus.fasta > $OUTPUTDIR/flye.nanopolish.medaka.racon1/consensus.racon1.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.medaka.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka.racon1/consensus.racon1.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.medaka.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka.racon1/ragtag/consensus.racon1.corrected.fasta

mkdir flye.nanopolish.medaka.racon2
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.nanopolish.medaka.racon1/consensus.racon1.fasta $OUTPUTDIR/flye.nanopolish.medaka.racon1/consensus.racon1.fasta > $OUTPUTDIR/flye.nanopolish.medaka.racon2/consensus.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.nanopolish.medaka.racon1/consensus.racon1.fasta $OUTPUTDIR/flye.nanopolish.medaka.racon2/consensus.sam $OUTPUTDIR/flye.nanopolish.medaka.racon1/consensus.racon1.fasta > $OUTPUTDIR/flye.nanopolish.medaka.racon2/consensus.racon2.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.medaka.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka.racon2/consensus.racon2.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.medaka.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka.racon2/ragtag/consensus.racon2.corrected.fasta

mkdir flye.nanopolish.medaka.racon3
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.nanopolish.medaka.racon2/consensus.racon2.fasta $OUTPUTDIR/flye.nanopolish.medaka.racon2/consensus.racon2.fasta > $OUTPUTDIR/flye.nanopolish.medaka.racon3/consensus.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.nanopolish.medaka.racon2/consensus.racon2.fasta $OUTPUTDIR/flye.nanopolish.medaka.racon3/consensus.sam $OUTPUTDIR/flye.nanopolish.medaka.racon2/consensus.racon2.fasta > $OUTPUTDIR/flye.nanopolish.medaka.racon3/consensus.racon3.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.medaka.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka.racon3/consensus.racon3.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.medaka.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka.racon3/ragtag/consensus.racon3.corrected.fasta

mkdir flye.nanopolish.medaka.racon4
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.nanopolish.medaka.racon3/consensus.racon3.fasta $OUTPUTDIR/flye.nanopolish.medaka.racon3/consensus.racon3.fasta > $OUTPUTDIR/flye.nanopolish.medaka.racon4/consensus.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.nanopolish.medaka.racon3/consensus.racon3.fasta $OUTPUTDIR/flye.nanopolish.medaka.racon4/consensus.sam $OUTPUTDIR/flye.nanopolish.medaka.racon3/consensus.racon3.fasta > $OUTPUTDIR/flye.nanopolish.medaka.racon4/consensus.racon4.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.medaka.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka.racon4/consensus.racon4.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.medaka.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.medaka.racon4/ragtag/consensus.racon4.corrected.fasta

mkdir flye.nanopolish.racon1
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.nanopolish/assembly.nanopolish.fasta $OUTPUTDIR/flye.nanopolish/assembly.nanopolish.fasta > $OUTPUTDIR/flye.nanopolish.racon1/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.nanopolish/assembly.nanopolish.fasta $OUTPUTDIR/flye.nanopolish.racon1/assembly.sam $OUTPUTDIR/flye.nanopolish/assembly.nanopolish.fasta > $OUTPUTDIR/flye.nanopolish.racon1/assembly.racon1.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon1/assembly.racon1.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon1/ragtag/assembly.racon1.corrected.fasta

mkdir flye.nanopolish.racon2
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.nanopolish.racon1/assembly.racon1.fasta $OUTPUTDIR/flye.nanopolish.racon1/assembly.racon1.fasta > $OUTPUTDIR/flye.nanopolish.racon2/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.nanopolish.racon1/assembly.racon1.fasta $OUTPUTDIR/flye.nanopolish.racon2/assembly.sam $OUTPUTDIR/flye.nanopolish.racon1/assembly.racon1.fasta > $OUTPUTDIR/flye.nanopolish.racon2/assembly.racon2.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon2/assembly.racon2.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon2/ragtag/assembly.racon2.corrected.fasta

mkdir flye.nanopolish.racon3
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.nanopolish.racon2/assembly.racon2.fasta $OUTPUTDIR/flye.nanopolish.racon2/assembly.racon2.fasta > $OUTPUTDIR/flye.nanopolish.racon3/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.nanopolish.racon2/assembly.racon2.fasta $OUTPUTDIR/flye.nanopolish.racon3/assembly.sam $OUTPUTDIR/flye.nanopolish.racon2/assembly.racon2.fasta > $OUTPUTDIR/flye.nanopolish.racon3/assembly.racon3.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon3/assembly.racon3.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon3/ragtag/assembly.racon3.corrected.fasta

mkdir flye.nanopolish.racon4
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.nanopolish.racon3/assembly.racon3.fasta $OUTPUTDIR/flye.nanopolish.racon3/assembly.racon3.fasta > $OUTPUTDIR/flye.nanopolish.racon4/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.nanopolish.racon3/assembly.racon3.fasta $OUTPUTDIR/flye.nanopolish.racon4/assembly.sam $OUTPUTDIR/flye.nanopolish.racon3/assembly.racon3.fasta > $OUTPUTDIR/flye.nanopolish.racon4/assembly.racon4.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon4/assembly.racon4.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon4/ragtag/assembly.racon4.corrected.fasta

mkdir flye.nanopolish.racon1.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.nanopolish.racon1/assembly.racon1.fasta -o $OUTPUTDIR/flye.nanopolish.racon1.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.racon1.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon1.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.racon1.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon1.medaka/ragtag/consensus.corrected.fasta

mkdir flye.nanopolish.racon2.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.nanopolish.racon2/assembly.racon2.fasta -o $OUTPUTDIR/flye.nanopolish.racon2.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.racon2.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon2.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.racon2.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon2.medaka/ragtag/consensus.corrected.fasta

mkdir flye.nanopolish.racon3.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.nanopolish.racon3/assembly.racon3.fasta -o $OUTPUTDIR/flye.nanopolish.racon3.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.racon3.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon3.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.racon3.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon3.medaka/ragtag/consensus.corrected.fasta

mkdir flye.nanopolish.racon4.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.nanopolish.racon4/assembly.racon4.fasta -o $OUTPUTDIR/flye.nanopolish.racon4.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.nanopolish.racon4.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon4.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.nanopolish.racon4.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.nanopolish.racon4.medaka/ragtag/consensus.corrected.fasta

mkdir flye.racon1
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye/assembly.fasta $OUTPUTDIR/flye/assembly.fasta > $OUTPUTDIR/flye.racon1/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye/assembly.fasta $OUTPUTDIR/flye.racon1/assembly.sam $OUTPUTDIR/flye/assembly.fasta > $OUTPUTDIR/flye.racon1/assembly.racon1.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.racon1/assembly.racon1.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon1/ragtag $REFSEQ $OUTPUTDIR/flye.racon1/ragtag/assembly.racon1.corrected.fasta

mkdir flye.racon2
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.racon1/assembly.racon1.fasta $OUTPUTDIR/flye.racon1/assembly.racon1.fasta > $OUTPUTDIR/flye.racon2/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.racon1/assembly.racon1.fasta $OUTPUTDIR/flye.racon2/assembly.sam $OUTPUTDIR/flye.racon1/assembly.racon1.fasta > $OUTPUTDIR/flye.racon2/assembly.racon2.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.racon2/assembly.racon2.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon2/ragtag $REFSEQ $OUTPUTDIR/flye.racon2/ragtag/assembly.racon2.corrected.fasta

mkdir flye.racon3
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.racon2/assembly.racon2.fasta $OUTPUTDIR/flye.racon2/assembly.racon2.fasta > $OUTPUTDIR/flye.racon3/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.racon2/assembly.racon2.fasta $OUTPUTDIR/flye.racon3/assembly.sam $OUTPUTDIR/flye.racon2/assembly.racon2.fasta > $OUTPUTDIR/flye.racon3/assembly.racon3.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.racon3/assembly.racon3.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon3/ragtag $REFSEQ $OUTPUTDIR/flye.racon3/ragtag/assembly.racon3.corrected.fasta

mkdir flye.racon4
source activate racon
minimap2 -t $THREADS -ax ava-ont $OUTPUTDIR/flye.racon3/assembly.racon3.fasta $OUTPUTDIR/flye.racon3/assembly.racon3.fasta > $OUTPUTDIR/flye.racon4/assembly.sam
racon -t $THREADS -f -u $OUTPUTDIR/flye.racon3/assembly.racon3.fasta $OUTPUTDIR/flye.racon4/assembly.sam $OUTPUTDIR/flye.racon3/assembly.racon3.fasta > $OUTPUTDIR/flye.racon4/assembly.racon4.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.racon4/assembly.racon4.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon4/ragtag $REFSEQ $OUTPUTDIR/flye.racon4/ragtag/assembly.racon4.corrected.fasta

mkdir flye.racon1.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.racon1/assembly.racon1.fasta -o $OUTPUTDIR/flye.racon1.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon1.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon1.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon1.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon1.medaka/ragtag/consensus.corrected.fasta

mkdir flye.racon2.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.racon2/assembly.racon2.fasta -o $OUTPUTDIR/flye.racon2.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon2.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon2.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon2.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon2.medaka/ragtag/consensus.corrected.fasta

mkdir flye.racon3.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.racon3/assembly.racon3.fasta -o $OUTPUTDIR/flye.racon3.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon3.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon3.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon3.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon3.medaka/ragtag/consensus.corrected.fasta

mkdir flye.racon4.medaka
source activate medaka
time medaka_consensus -t $THREADS -i $SAMPLE -d flye.racon4/assembly.racon4.fasta -o $OUTPUTDIR/flye.racon4.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon4.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon4.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon4.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon4.medaka/ragtag/consensus.corrected.fasta

mkdir flye.racon1.medaka.nanopolish
mkdir flye.racon1.medaka.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.racon1.medaka/consensus.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.racon1.medaka.nanopolish/consensus.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.racon1.medaka.nanopolish/consensus.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.racon1.medaka/consensus.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.racon1.medaka.nanopolish/consensus.sorted.bam -o $OUTPUTDIR/flye.racon1.medaka.nanopolish/vcf/consensus.{#}.vcf -g $OUTPUTDIR/flye.racon1.medaka/consensus.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.racon1.medaka.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.racon1.medaka.nanopolish/consensus.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.racon1.medaka/consensus.fasta -f $OUTPUTDIR/flye.racon1.medaka.nanopolish/consensus.vcflist.txt > $OUTPUTDIR/flye.racon1.medaka.nanopolish/consensus.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon1.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon1.medaka.nanopolish/consensus.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon1.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon1.medaka.nanopolish/ragtag/consensus.nanopolish.corrected.fasta

mkdir flye.racon2.medaka.nanopolish
mkdir flye.racon2.medaka.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.racon2.medaka/consensus.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.racon2.medaka.nanopolish/consensus.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.racon2.medaka.nanopolish/consensus.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.racon2.medaka/consensus.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.racon2.medaka.nanopolish/consensus.sorted.bam -o $OUTPUTDIR/flye.racon2.medaka.nanopolish/vcf/consensus.{#}.vcf -g $OUTPUTDIR/flye.racon2.medaka/consensus.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.racon2.medaka.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.racon2.medaka.nanopolish/consensus.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.racon2.medaka/consensus.fasta -f $OUTPUTDIR/flye.racon2.medaka.nanopolish/consensus.vcflist.txt > $OUTPUTDIR/flye.racon2.medaka.nanopolish/consensus.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon2.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon2.medaka.nanopolish/consensus.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon2.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon2.medaka.nanopolish/ragtag/consensus.nanopolish.corrected.fasta

mkdir flye.racon3.medaka.nanopolish
mkdir flye.racon3.medaka.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.racon3.medaka/consensus.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.racon3.medaka.nanopolish/consensus.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.racon3.medaka.nanopolish/consensus.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.racon3.medaka/consensus.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.racon3.medaka.nanopolish/consensus.sorted.bam -o $OUTPUTDIR/flye.racon3.medaka.nanopolish/vcf/consensus.{#}.vcf -g $OUTPUTDIR/flye.racon3.medaka/consensus.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.racon3.medaka.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.racon3.medaka.nanopolish/consensus.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.racon3.medaka/consensus.fasta -f $OUTPUTDIR/flye.racon3.medaka.nanopolish/consensus.vcflist.txt > $OUTPUTDIR/flye.racon3.medaka.nanopolish/consensus.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon3.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon3.medaka.nanopolish/consensus.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon3.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon3.medaka.nanopolish/ragtag/consensus.nanopolish.corrected.fasta

mkdir flye.racon4.medaka.nanopolish
mkdir flye.racon4.medaka.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.racon4.medaka/consensus.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.racon4.medaka.nanopolish/consensus.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.racon4.medaka.nanopolish/consensus.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.racon4.medaka/consensus.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.racon4.medaka.nanopolish/consensus.sorted.bam -o $OUTPUTDIR/flye.racon4.medaka.nanopolish/vcf/consensus.{#}.vcf -g $OUTPUTDIR/flye.racon4.medaka/consensus.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.racon4.medaka.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.racon4.medaka.nanopolish/consensus.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.racon4.medaka/consensus.fasta -f $OUTPUTDIR/flye.racon4.medaka.nanopolish/consensus.vcflist.txt > $OUTPUTDIR/flye.racon4.medaka.nanopolish/consensus.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon4.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon4.medaka.nanopolish/consensus.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon4.medaka.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon4.medaka.nanopolish/ragtag/consensus.nanopolish.corrected.fasta

mkdir flye.racon1.nanopolish
mkdir flye.racon1.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.racon1/assembly.racon1.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.racon1.nanopolish/assembly.racon1.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.racon1.nanopolish/assembly.racon1.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.racon1/assembly.racon1.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.racon1.nanopolish/assembly.racon1.sorted.bam -o $OUTPUTDIR/flye.racon1.nanopolish/vcf/assembly.{#}.vcf -g $OUTPUTDIR/flye.racon1/assembly.racon1.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.racon1.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.racon1.nanopolish/assembly.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.racon1/assembly.racon1.fasta -f $OUTPUTDIR/flye.racon1.nanopolish/assembly.vcflist.txt > $OUTPUTDIR/flye.racon1.nanopolish/assembly.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon1.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon1.nanopolish/assembly.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon1.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon1.nanopolish/ragtag/assembly.nanopolish.corrected.fasta

mkdir flye.racon2.nanopolish
mkdir flye.racon2.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.racon2/assembly.racon2.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.racon2.nanopolish/assembly.racon2.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.racon2.nanopolish/assembly.racon2.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.racon2/assembly.racon2.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.racon2.nanopolish/assembly.racon2.sorted.bam -o $OUTPUTDIR/flye.racon2.nanopolish/vcf/assembly.{#}.vcf -g $OUTPUTDIR/flye.racon2/assembly.racon2.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.racon2.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.racon2.nanopolish/assembly.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.racon2/assembly.racon2.fasta -f $OUTPUTDIR/flye.racon2.nanopolish/assembly.vcflist.txt > $OUTPUTDIR/flye.racon2.nanopolish/assembly.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon2.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon2.nanopolish/assembly.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon2.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon2.nanopolish/ragtag/assembly.nanopolish.corrected.fasta

mkdir flye.racon3.nanopolish
mkdir flye.racon3.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.racon3/assembly.racon3.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.racon3.nanopolish/assembly.racon3.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.racon3.nanopolish/assembly.racon3.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.racon3/assembly.racon3.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.racon3.nanopolish/assembly.racon3.sorted.bam -o $OUTPUTDIR/flye.racon3.nanopolish/vcf/assembly.{#}.vcf -g $OUTPUTDIR/flye.racon3/assembly.racon3.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.racon3.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.racon3.nanopolish/assembly.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.racon3/assembly.racon3.fasta -f $OUTPUTDIR/flye.racon3.nanopolish/assembly.vcflist.txt > $OUTPUTDIR/flye.racon3.nanopolish/assembly.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon3.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon3.nanopolish/assembly.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon3.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon3.nanopolish/ragtag/assembly.nanopolish.corrected.fasta

mkdir flye.racon4.nanopolish
mkdir flye.racon4.nanopolish/vcf
source activate nanopolish
minimap2 -t $THREADS -ax map-ont $OUTPUTDIR/flye.racon4/assembly.racon4.fasta $SAMPLE | samtools sort -@ $THREADS -o $OUTPUTDIR/flye.racon4.nanopolish/assembly.racon4.sorted.bam -
samtools index -@ $THREADS $OUTPUTDIR/flye.racon4.nanopolish/assembly.racon4.sorted.bam
nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE
nanopolish_makerange.py $OUTPUTDIR/flye.racon4/assembly.racon4.fasta --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $OUTPUTDIR/flye.racon4.nanopolish/assembly.racon4.sorted.bam -o $OUTPUTDIR/flye.racon4.nanopolish/vcf/assembly.{#}.vcf -g $OUTPUTDIR/flye.racon4/assembly.racon4.fasta -w {1} -p $PLOIDY -v
ls $OUTPUTDIR/flye.racon4.nanopolish/vcf/*.vcf > $OUTPUTDIR/flye.racon4.nanopolish/assembly.vcflist.txt
nanopolish vcf2fasta --skip-checks -g $OUTPUTDIR/flye.racon4/assembly.racon4.fasta -f $OUTPUTDIR/flye.racon4.nanopolish/assembly.vcflist.txt > $OUTPUTDIR/flye.racon4.nanopolish/assembly.nanopolish.fasta
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon4.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon4.nanopolish/assembly.nanopolish.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon4.nanopolish/ragtag $REFSEQ $OUTPUTDIR/flye.racon4.nanopolish/ragtag/assembly.nanopolish.corrected.fasta

mkdir flye.racon1.nanopolish.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.racon1.nanopolish/assembly.nanopolish.fasta -o $OUTPUTDIR/flye.racon1.nanopolish.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon1.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon1.nanopolish.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon1.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon1.nanopolish.medaka/ragtag/consensus.corrected.fasta

mkdir flye.racon2.nanopolish.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.racon2.nanopolish/assembly.nanopolish.fasta -o $OUTPUTDIR/flye.racon2.nanopolish.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon2.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon2.nanopolish.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon2.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon2.nanopolish.medaka/ragtag/consensus.corrected.fasta

mkdir flye.racon3.nanopolish.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.racon3.nanopolish/assembly.nanopolish.fasta -o $OUTPUTDIR/flye.racon3.nanopolish.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon3.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon3.nanopolish.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon3.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon3.nanopolish.medaka/ragtag/consensus.corrected.fasta

mkdir flye.racon4.nanopolish.medaka
source activate medaka
medaka_consensus -t $THREADS -i $SAMPLE -d flye.racon4.nanopolish/assembly.nanopolish.fasta -o $OUTPUTDIR/flye.racon4.nanopolish.medaka
source activate ragtag
ragtag.py correct -t $THREADS -u -o $OUTPUTDIR/flye.racon4.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon4.nanopolish.medaka/consensus.fasta
ragtag.py scaffold -t $THREADS -u -C -o $OUTPUTDIR/flye.racon4.nanopolish.medaka/ragtag $REFSEQ $OUTPUTDIR/flye.racon4.nanopolish.medaka/ragtag/consensus.corrected.fasta
