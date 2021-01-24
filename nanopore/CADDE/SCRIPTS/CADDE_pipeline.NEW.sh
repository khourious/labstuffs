csv="$1"	# csv file name corresponds to the library name with this format NameVirus(es)_LocationSequencing_libraryNumber_Date   e.g. YFV_SaoPaulo_library1_20190305.csv
		# csv file contains this format sample;barcode;virus_reference, NO HEADER !!	 e.g. sample01;BC01;YFV500

#librarypath = "$2"

library=`echo "$csv" | sed -e 's/\.csv//g' -e 's/.*\///g'`


cd ~/Desktop/LIBRARIES/"$library"

### Guppy commad just in case:
#/usr/local/guppy-gpu/ont-guppy/bin/guppy_basecaller -i RAW_FILES -s BASECALLED_FILES -c dna_r9.4.1_450bps_fast.cfg -q 0 -r --device cuda:00000000:01:00.0 --disable_pings

#rm final_summary.txt

source activate artic-workshop

#~/Desktop/SCRIPTS/clean_raw_files.sh

~/Desktop/SCRIPTS/gather.sh "$csv"

~/Desktop/SCRIPTS/porechop.parallel.sh "$library"

echo "Sample@Nb of reads mapped@Average depth coverage@Bases covered >10x@Bases covered >25x@Reference covered (%)" | tr '@' '\t' > CONSENSUS/"$library".stats.txt

for i in `cat "$csv"`

	do

		sample=`echo "$i" | awk -F"," '{print $1}' | sed '/^$/d'`
		barcode=`echo "$i"| awk -F"," '{print $2}' | sed '/^$/d'`
		ref=`echo "$i" | awk -F"," '{print $3}' | sed '/^$/d'`

		~/Desktop/SCRIPTS/consensus.sh "$library" "$ref" "$barcode" "$sample"

		~/Desktop/SCRIPTS/stats.sh "$sample" "$library" "$ref"

		sed -e '/>/ s/\.primertrimmed\.sorted\.bam//g' CONSENSUS/"$sample".consensus.fasta

        done > CONSENSUS/"$library".consensus.fasta

cp CONSENSUS/"$library".consensus.fasta /home/artic/Desktop/CADDE_RESULTS
cp CONSENSUS/"$library".stats.txt /home/artic/Desktop/CADDE_RESULTS
tar zvcf ~/Desktop/CADDE_RESULTS/"$library".tar.gz ~/Desktop/LIBRARIES/"$library"

echo "########################"
echo "####### Finish!! #######"
echo "########################"
