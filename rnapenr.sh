#!/bin/bash

# mamba create -y -n rnapenr -c conda-forge -c anaconda -c bioconda -c defaults blast bwa entrez-direct exonerate fastp ivar mafft samtools seqkit seqtk spades

bg() {

    start=$(date +%s.%N)

    source activate rnapenr

    # DATABASE=$HOME/rnapenr/blastdb/multifasta.fasta
    # makeblastdb -in "$DATABASE" -dbtype nucl -out $HOME/rnapenr/blastdb/"$(basename "$DATABASE" | awk -F. '{print $1}')"
    # awk '/^>/{if(s) close(s); s="'"$HOME"'/rnapenr/refseq/" substr($0,2) ".fasta"; print > s; next}{if(s) print >> s}' "$DATABASE"

    THREADS=$(nproc)
    SAMPLE_SHEET=$1
    DEPTH=10
    INPUT_DIR=$HOME/BaseSpace
    OUTPUT_DIR=$HOME/rnapenr/assembly

    [[ ! -d "$OUTPUT_DIR" ]] && mkdir "$OUTPUT_DIR" && chmod 700 -R "$OUTPUT_DIR"

    # cat $SAMPLE_SHEET | tr -dc '[:print:]\n' | sed -e '1,18d' | awk -v SAMPLE_SHEET="$(basename "$SAMPLE_SHEET" .csv)" -F, '{print $1","SAMPLE_SHEET","$NF}' | sort > "$OUTPUT_DIR"/sample_sheet.csv

    for i in $(cat "$OUTPUT_DIR"/sample_sheet.csv); do
        LIBRARY=$(echo "$i" | awk -F, '{print $2}' | sort -u)
        bs download project --no-metadata --summary --extension=fastq.gz -o "$INPUT_DIR"/"$LIBRARY" -n "$LIBRARY"
    done

    for i in $(cat "$OUTPUT_DIR"/sample_sheet.csv); do
        SAMPLE_ID=$(echo "$i" | awk -F, '{print $1}')
        echo "ref_seq#target#sample_id#num_total_reads#num_mapp_reads#avg_depth#depth_10x#depth_100x#depth_1000x#ref_cov#ncount#ncount_perc" | tr '#' '\t' > "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
    done

    for i in $(cat "$OUTPUT_DIR"/sample_sheet.csv); do
        SAMPLE_ID=$(echo "$i" | awk -F, '{print $1}')
        LIBRARY=$(echo "$i" | awk -F, '{print $2}')
        PANEL=$(echo "$i" | awk -F, '{print $3}')
        if [[ $(find "$INPUT_DIR"/"$LIBRARY" -type f -name '*_L00*') ]]; then
            cp -v "$INPUT_DIR"/"$LIBRARY"/"$SAMPLE_ID"_*/"$SAMPLE_ID"_*_R1_001.fastq.gz "$OUTPUT_DIR"/"$SAMPLE_ID".R1.fastq.gz
            cp -v "$INPUT_DIR"/"$LIBRARY"/"$SAMPLE_ID"_*/"$SAMPLE_ID"_*_R2_001.fastq.gz "$OUTPUT_DIR"/"$SAMPLE_ID".R2.fastq.gz
            mkdir "$OUTPUT_DIR"/"$SAMPLE_ID"
            fastp --cut_front --cut_tail --qualified_quality_phred 20 -l 75 --thread "$THREADS" -f 0 -t 0 -F 0 -T 0 -i "$OUTPUT_DIR"/"$SAMPLE_ID".R1.fastq.gz -I "$OUTPUT_DIR"/"$SAMPLE_ID".R2.fastq.gz -o "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID".trimmed.R1.fastq.gz -O "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID".trimmed.R2.fastq.gz -h "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID".trimmed.html -j "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID".trimmed.json
            rnaspades.py -o "$OUTPUT_DIR"/"$SAMPLE_ID"/denovo -1 "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID".trimmed.R1.fastq.gz -2 "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID".trimmed.R2.fastq.gz -t "$THREADS" -m 15 -k 127
            if [ -f "$OUTPUT_DIR"/"$SAMPLE_ID"/denovo/transcripts.fasta ]; then
                blastn -query "$OUTPUT_DIR"/"$SAMPLE_ID"/denovo/transcripts.fasta -db $HOME/rnapenr/blastdb/"$PANEL" -out "$OUTPUT_DIR"/"$SAMPLE_ID"/denovo/results.txt -outfmt "6 qseqid sseqid pident length"
            else
                blastn -query "$OUTPUT_DIR"/"$SAMPLE_ID"/denovo/before_rr.fasta -db $HOME/rnapenr/blastdb/"$PANEL" -out "$OUTPUT_DIR"/"$SAMPLE_ID"/denovo/results.txt -outfmt "6 qseqid sseqid pident length"
            fi
        fi
        for j in $(awk -F"\t" '$3 > 98 {print $2}' "$OUTPUT_DIR"/"$SAMPLE_ID"/denovo/results.txt | sort -u); do
            bwa index $HOME/rnapenr/refseq/"$j".fasta
            bwa mem -t "$THREADS" $HOME/rnapenr/refseq/"$j".fasta "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID".trimmed.R1.fastq.gz "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID".trimmed.R2.fastq.gz -o "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$j".bam
            samtools sort -o "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$j".sorted.bam "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$j".bam
            samtools index "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$j".sorted.bam
            samtools mpileup -d 50000 --reference $HOME/rnapenr/refseq/"$j".fasta -a -B "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$j".sorted.bam | ivar variants -p "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$j" -q 30 -t 0.05 -r $HOME/rnapenr/refseq/"$j".fasta
            samtools mpileup -d 50000 --reference $HOME/rnapenr/refseq/"$j".fasta -a -B "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$j".sorted.bam | ivar consensus -p "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$j".depth"$DEPTH" -q 30 -t 0 -m 10 -n N
            sed -i -e 's/>.*/>'${SAMPLE_ID}'/g' "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$j".depth"$DEPTH".fa
        done
        rm -rf "$OUTPUT_DIR"/"$SAMPLE_ID".*.fastq.gz
        for k in $(awk -F"\t" '$3 > 98 {print $2}' "$OUTPUT_DIR"/"$SAMPLE_ID"/denovo/results.txt | sort -u); do
            echo "$k" | awk '{printf $0"#"}' | tr '#' '\t' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            grep "^$k" $HOME/rnapenr/blastdb/"$PANEL".tsv | awk -F"\t" '{print $2}' | awk '{printf $0"#"}' | tr '#' '\t' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            echo "$SAMPLE_ID" | awk '{printf $0"#"}' | tr '#' '\t' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            samtools view -c "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$k".sorted.bam | awk '{printf $0"#"}' | tr '#' '\t' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            samtools view -c -h -F 4 "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$k".sorted.bam | awk '{printf $0"#"}' | tr '#' '\t' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            AVG_DEPTH=$(samtools depth "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$k".sorted.bam | awk '{sum+=$3} END {print sum/NR}')
            if [[ "$AVG_DEPTH" == "" || "$AVG_DEPTH" == 0 ]]; then
                echo "0.00""#" | tr '#' '\t' | tr -d '\n' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            else
                echo "$AVG_DEPTH" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            fi
            paste <(samtools depth "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$k".sorted.bam | awk '{if ($3 > '"10"') {print $0}}' | wc -l) <(fastalength $HOME/rnapenr/refseq/"$k".fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            paste <(samtools depth "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$k".sorted.bam | awk '{if ($3 > '"100"') {print $0}}' | wc -l) <(fastalength $HOME/rnapenr/refseq/"$k".fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            paste <(samtools depth "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$k".sorted.bam | awk '{if ($3 > '"1000"') {print $0}}' | wc -l) <(fastalength $HOME/rnapenr/refseq/"$k".fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            N_COUNT=$(seqtk comp "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$k".depth10.fa | awk -F"\t" '{print $9}')
            N_COUNT_PER=$(paste <(seqtk comp "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$k".depth10.fa | awk -F"\t" '{print $9}') <(fastalength $HOME/rnapenr/refseq/"$k".fasta | awk '{print $1}')| awk -F"\t" '{printf("%0.2f\n", ($1/$2)*100)}')
            REF_SEQ_LENGHT=$(fastalength $HOME/rnapenr/refseq/"$k".fasta | awk -F" " '{print $1}')
            REF_SEQ_COV=$(paste <(fastalength $HOME/rnapenr/refseq/"$k".fasta | awk '{print $1}') <(seqtk comp "$OUTPUT_DIR"/"$SAMPLE_ID"/"$SAMPLE_ID"."$k".depth10.fa | awk -F"\t" '{print $9}') | awk -F"\t" '{printf("%0.2f\n", ($1-$2)/$1*100)}')
            if [[ "$N_COUNT" == 0 ]]; then
                echo "0.00" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
                echo "$REF_SEQ_LENGHT" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
                echo "100.00" | awk '{printf $0"\n"}' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            else
                echo "$REF_SEQ_COV" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
                echo "$N_COUNT" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
                echo "$N_COUNT_PER" | awk '{printf $0"\n"}' >> "$OUTPUT_DIR"/"$SAMPLE_ID".summary.tsv
            fi
        done
    done

    conda deactivate

    end=$(date +%s.%N)
    runtime=$(python -c "print(${end} - ${start})")
    echo "" && echo "Done. The runtime was "$runtime" seconds." && echo ""

}

bg $1 >> rnapenr.log.$(uname -n).$(date +'%Y-%m-%d').txt 2>&1 &

exit 0

