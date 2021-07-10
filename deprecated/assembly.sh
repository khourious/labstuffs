# create a reference asembly index
bwa index $reference.fasta

# align and sort a BAM file
bwa mem $reference.fasta $sample.fastq -t 12 | samtools view -@ 12 -bS - | samtools sort -@ 12 -o $sample.sorted.bam -

# create a BAM index file
samtools index $sample.sorted.bam -@ 12

# convert a BAM file to a SAM file
samtools view -h $sample.sorted.bam -o $sample.sorted.sam -@ 12

# view stats to a BAM file
samtools flagstat $sample.sorted.bam > $sample.sorted.flagstat.txt

# filtering out UNMAPPED reads in BAM files
samtools view -h -F 4 -b $sample.sorted.bam > $sample.sorted.mapped.bam -@ 12

# view stats to a BAM file
samtools flagstat $sample.sorted.mapped.bam > $sample.sorted.mapped.flagstat.txt

# conver a BAM file to a SAM file
samtools view -h $sample.sorted.mapped.bam -o $sample.sorted.mapped.sam -@ 12

# filtering out MAPPED reads in BAM files
samtools view -bS -f 4 $sample.sorted.bam > $sample.sorted.unmapped.bam -@ 12

# view stats to a BAM file
samtools flagstat $sample.sorted.unmapped.bam > $sample.sorted.unmapped.flagstat.txt

# convert a BAM file to a SAM file
samtools view -h $sample.sorted.unmapped.bam -o $sample.sorted.unmapped.sam -@ 12


# create a variant calling
bcftools mpileup -Ob -f $reference.fasta $sample.sorted.mapped.bam -t 12 | bcftools call -mv -Oz -o $sample.sorted.mapped.vcf.gz

# normalize indels
bcftools norm -f $reference.fasta $sample.sorted.mapped.vcf.gz -Oz -o $sample.sorted.mapped.norm.vcf.gz

# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 $sample.sorted.mapped.norm.vcf.gz -Oz -o $sample.sorted.mapped.flt-indels.vcf.gz

# index VCF file
bcftools index $sample.sorted.mapped.flt-indels.vcf.gz

# produce a consensus FASTA file
cat $reference.fasta | bcftools consensus $sample.sorted.mapped.flt-indels.vcf.gz > $sample.consensus.fasta

# get coverage and depth
samtools coverage $sample.sorted.mapped.bam > $sample.coverage
samtools depth $sample.sorted.mapped.bam > $sample.depth