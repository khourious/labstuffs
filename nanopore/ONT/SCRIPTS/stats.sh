sample="$1"

library="$2"

pathref="$3"

ref=`echo "$pathref" | sed -e 's/\/.*//g'`

echo -n "$sample""@" | tr '@' '\t' >> "$library".stats.txt

$HOME/CADDE/SCRIPTS/nb_reads_mapped.sh "$sample".sorted.bam | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

$HOME/CADDE/SCRIPTS/depth.sh "$sample".sorted.bam | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

$HOME/CADDE/SCRIPTS/bases_covered_x.sh "$sample".sorted.bam 10 | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

$HOME/CADDE/SCRIPTS/bases_covered_x.sh "$sample".sorted.bam 25 | awk '{printf $1"@"}' | tr '@' '\t' >> "$library".stats.txt

$HOME/CADDE/SCRIPTS/coverage.sh $HOME/CADDE/PRIMER_SCHEMES/"$pathref"/"$ref".reference.fasta "$sample".consensus.fasta | awk '{print $2}' >> "$library".stats.txt
