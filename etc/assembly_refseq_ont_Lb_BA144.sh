#!/bin/bash

# author: Laise de Moraes <laisepaixao@live.com>
# institution: Universidade Federal da Bahia, Brazil
# URL: https://github.com/lpmor22
# date: 024 JAN 2021

THREADS="12" #12

RAWDIR="/mnt/x/Lb_BA144/library15-2019-11-15"
HACDEMUXDIR="/mnt/x/Lb_BA144/Lb_BA144_HAC_Demux"

SAMPLE="/mnt/x/Lb_BA144/Lb_BA144.fastq"
SAMPLEID="Lb_BA144"

REFSEQ="/mnt/x/Lb_BA144/TriTrypDB/TriTrypDB-50_LbraziliensisMHOMBR75M2904_2019_Genome.fasta"
PLOIDY="2" #2

# conda create -y -n minimap2 -c conda-forge -c bioconda -c defaults minimap2 samtools
mamba create -y -n minimap2 -c conda-forge -c bioconda -c defaults minimap2 samtools
# conda create -y -n nanopolish -c conda-forge -c bioconda -c defaults nanopolish samtools
mamba create -y -n nanopolish -c conda-forge -c bioconda -c defaults nanopolish samtools
# conda create -y -n medaka -c conda-forge -c bioconda -c defaults medaka bcftools minimap2 samtools
mamba create -y -n medaka -c conda-forge -c bioconda -c defaults medaka bcftools minimap2 samtools

source activate minimap2
minimap2 -t $THREADS -ax map-ont $REFSEQ $SAMPLE | samtools sort -@ $THREADS -o $SAMPLEID.sorted.bam -
samtools view -@ $THREADS -h -F 4 -b $SAMPLEID.sorted.bam > $SAMPLEID.sorted.mapped.bam
# samtools view -@ $THREADS -bS -f 4 $SAMPLEID.sorted.bam > $SAMPLEID.sorted.unmapped.bam
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
