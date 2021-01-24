library="$1"

cd CONSENSUS

if [ -d "$library"_TMPDIR ] ; then rm -r "$library"_TMPDIR ; mkdir "$library"_TMPDIR ; else mkdir "$library"_TMPDIR ; fi

fastasplit -f "$library".fasta -o "$library"_TMPDIR -c 100 

ls "$library"_TMPDIR/"$library".fasta_chunk_00000* | parallel bash ~/Desktop/SCRIPTS/porechop.cmd {}

for i in `ls "$library"_TMPDIR/*_DIR/*.fasta | sed -e 's/.*\///g'| sort | uniq | sed -e 's/\.fasta//g'` ; do j=`echo "$i" | sed -e 's/NB/BC/g'` ; cat "$library"_TMPDIR/*/"$i".fasta > "$j"_"$library".fasta ; done

#rm none_"$library".fasta
rm -r "$library"_TMPDIR

cd ..
