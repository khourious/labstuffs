sample="$1"

library="$2"

pathref="$3"

ref=`echo "$pathref" | sed -e 's/\/.*//g'`

echo -n "$sample""@" | tr '@' '\t' >> "$library".stats.txt

$HOME/ONT/SCRIPTS/nb_reads_mapped.sh "$sample".sorted.bam | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

$HOME/ONT/SCRIPTS/depth.sh "$sample".sorted.bam | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

$HOME/ONT/SCRIPTS/bases_covered_x.sh "$sample".sorted.bam 10 | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

$HOME/ONT/SCRIPTS/bases_covered_x.sh "$sample".sorted.bam 25 | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

$HOME/ONT/SCRIPTS/coverage.sh $HOME/ONT/PRIMER_SCHEMES/"$pathref"/"$ref".reference.fasta "$sample".consensus.fasta | awk '{print $2}' >> "$library".stats.txt
