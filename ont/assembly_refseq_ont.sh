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
    mamba create -y -n minimap2 -c conda-forge -c anaconda -c bioconda -c defaults minimap2 samtools
    mamba create -y -n nanopolish -c conda-forge -c anaconda -c bioconda -c defaults nanopolish samtools
else
    if [[ -z "$(which mamba)" ]]; then
        conda install -y -c conda-forge mamba
        mamba update -y -n base conda
        mamba create -y -n minimap2 -c conda-forge -c anaconda -c bioconda -c defaults minimap2 samtools
        mamba create -y -n nanopolish -c conda-forge -c anaconda -c bioconda -c defaults nanopolish samtools
    else
        mamba update -y -n base conda
        mamba create -y -n minimap2 -c conda-forge -c anaconda -c bioconda -c defaults minimap2 samtools
        mamba create -y -n nanopolish -c conda-forge -c anaconda -c bioconda -c defaults nanopolish samtools
    fi
fi

THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

RAWDIR="/path/directory"
HACDEMUXDIR="/path/directory"

SAMPLE=""
SAMPLEID=""
REFSEQ=""
PLOIDY=""

source activate minimap2

minimap2 -t $THREADS -ax map-ont $REFSEQ $SAMPLE | samtools sort -@ $THREADS -o $SAMPLEID.sorted.bam -

samtools view -@ $THREADS -h -F 4 -b $SAMPLEID.sorted.bam > $SAMPLEID.sorted.mapped.bam

samtools view -@ $THREADS -bS -f 4 $SAMPLEID.sorted.bam > $SAMPLEID.sorted.unmapped.bam

samtools index -@ $THREADS $SAMPLEID.sorted.mapped.bam

source activate nanopolish

nanopolish index -d $RAWDIR -s $HACDEMUXDIR/sequencing_summary.txt $SAMPLE

mkdir nanopolish

nanopolish_makerange.py $REFSEQ --overlap-length -1 | parallel -P 1 nanopolish variants -t $THREADS -r $SAMPLE -b $SAMPLEID.sorted.mapped.bam -o nanopolish/$SAMPLEID.{#}.vcf -g $REFSEQ -w {1} -p $PLOIDY -v

ls nanopolish/*.vcf > $SAMPLEID.vcflist.txt

nanopolish vcf2fasta --skip-checks -g $REFSEQ -f $SAMPLEID.vcflist.txt > $SAMPLEID.consensus.fasta

source activate minimap2

samtools coverage $SAMPLEID.sorted.mapped.bam > $SAMPLEID.sorted.mapped.coverage

samtools coverage $SAMPLEID.sorted.mapped.bam -A > $SAMPLEID.sorted.mapped.histogram.coverage
