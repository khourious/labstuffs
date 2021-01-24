location="../RAW_FILES/"	### path to library reads -- need full path with "/" at the end
library="$1"			### name of the library
ref="$2"			### reference sequence and primer scheme basename
barcode="$3"    		### barcode number : e.g. BC01
sample_name="$4"		### name of the sample corrspodnng to the barcode and library

cd CONSENSUS

if [ `awk 'NR==1{print $1}' "$library"_sequencing_summary.txt`  == "filename_fastq" ] 

	then cut -d$'\t' -f2- "$library"_sequencing_summary.txt | sed -e 's/^filename_fast5/filename/g' > tmp ; mv tmp "$library"_sequencing_summary.txt

fi 

grep ">" "$barcode"_"$library".fasta | awk '{print $1}' | sed -e 's/>//g' | sort -k1,1 > "$barcode"_"$library".txt
					
sed '1d' "$library"_sequencing_summary.txt | awk '{print $2,$1}' | sort -k1,1 > "$barcode"_"$library".link.txt
								
join -11 -11 "$barcode"_"$library".txt "$barcode"_"$library".link.txt | awk '{print $1"@""'$location'"$2}' | tr '@' '\t' > "$barcode"_"$library".fast5.txt


nanopolish index -d "$location" -s "$library"_sequencing_summary.txt "$barcode"_"$library".fasta

mv "$barcode"_"$library".fast5.txt "$barcode"_"$library".fasta.index.readdb

artic minion --normalise 200 --threads 10 --scheme-directory ~/Desktop/PRIMER_SCHEMES --read-file "$barcode"_"$library".fasta --nanopolish-read-file "$barcode"_"$library".fasta "$ref" "$sample_name"

cd ..

rm ~/Desktop/PRIMER_SCHEMES/"$ref".reference.fasta.*
