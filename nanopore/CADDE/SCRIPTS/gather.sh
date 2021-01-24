csv="$1"

library=`echo "$csv" | sed -e 's/\.csv//g' -e 's/.*\///g'`

if [ -d CONSENSUS ] ; then rm -r CONSENSUS ; mkdir CONSENSUS ; else mkdir CONSENSUS ; fi

nb_ref=`cat "$csv" | awk -F"," '{print $3}' | sed '/^$/d' | sort | uniq | wc -l`

if [ "$nb_ref" -eq 1 ] 
	
	then 
                ref=`cat "$csv" | awk -F"," '{print $3}' | sed '/^$/d' | sort | uniq`
                min=`awk -F"\t" '{print $2,$3}' ~/Desktop/PRIMER_SCHEMES/"$ref".scheme.bed | tr '\n' ' ' | awk '{for (i=1;i<=(NF/2);i=i+2) {print $(i*2+1)-$(i*2)}}' | sort -n | awk 'NR==1{print}' | awk '{print $1-50}'`

                max=`awk -F"\t" '{print $2,$3}' ~/Desktop/PRIMER_SCHEMES/"$ref".scheme.bed | tr '\n' ' ' | awk '{for (i=1;i<=(NF/2);i=i+2) {print $(i*2+1)-$(i*2)}}' | sort -nr | awk 'NR==1{print}' | awk '{print $1+300}'`
	else
		min=100 ; max=2000
fi

artic gather --min-length "$min" --max-length "$max" --guppy --prefix "$library" BASECALLED_FILES

mv "$library"_all.fastq CONSENSUS/"$library".fastq ; rm *.fastq

mv "$library"_sequencing_summary.txt CONSENSUS/

seqtk seq -A CONSENSUS/"$library".fastq > CONSENSUS/"$library".fasta ; rm CONSENSUS/"$library".fastq

