sample="$1"

library="$2"

pathref="$3"

ref=`echo "$pathref" | sed -e 's/\/.*//g'`

echo -n "$sample""@" | tr '@' '\t' >> "$library".stats.txt

nb_reads_mapped.sh "$sample".trimmed.sorted.bam | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

depth.sh "$sample".trimmed.sorted.bam | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

bases_covered_x.sh "$sample".trimmed.sorted.bam 10 | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

bases_covered_x.sh "$sample".trimmed.sorted.bam 25 | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

coverage.sh "$ref".reference.fasta "$sample".consensus.fa | awk '{print $2}' >> "$library".stats.txt