mkdir $1

cat $1_Illumina_ARTIC_V4/$1.R1.fq.gz $1_Illumina_FIOCRUZ-IOC_V2/$1.R1.fq.gz > $1/$1_merge.R1.fq.gz

cat $1_Illumina_ARTIC_V4/$1.R2.fq.gz $1_Illumina_FIOCRUZ-IOC_V2/$1.R2.fq.gz > $1/$1_merge.R2.fq.gz

cd $1

source activate igm-sars2_assembly

python $HOME/IGM_SARSCOV2/scripts/bwa_index.py -in /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta

python $HOME/IGM_SARSCOV2/scripts/bwa_mem.py -f /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta -pr $1_merge -p 12

python $HOME/IGM_SARSCOV2/scripts/ivar.py -f /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta -pr $1_merge -dp 10

python $HOME/IGM_SARSCOV2/scripts/get_mvs.py -f /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta -pr $1_merge -dp 10 -p 12 -di 10

source activate igm-sars2_summary

echo "biobanco_seq#num_total_reads#num_mapp_reads#avg_depth#depth_10x#depth_100x#depth_1000x#ref_cov#ncount#ncount_perc#pango_ver#pango_learn_ver#pango_lin#nextclade_ver#clade#nucl_substitutions#nucl_deletions#nucl_inserc#nucl_missing#aa_substitutions#aa_deletions" | tr '#' '\t' > $1_merge.summary.tmp

echo -n "#" | tr '#' '\n' >> $1_merge.summary.tmp
echo -n "$1_merge""#" | tr '#' '\t' >> $1_merge.summary.tmp

samtools view -c $1_merge.sorted.bam | awk '{printf $0"#"}' | tr '#' '\t' >> $1_merge.summary.tmp

samtools view -c -h -F 4 $1_merge.sorted.bam | awk '{printf $0"#"}' | tr '#' '\t' >> $1_merge.summary.tmp

AVGDEPTH=$(samtools depth $1_merge.sorted.bam | awk '{sum+=$3} END {print sum/NR}')

if [[ "$AVGDEPTH" == 0 ]]; then
    echo "0""#" | tr '#' '\t' | tr -d '\n' >> $1_merge.summary.tmp
else
    echo "$AVGDEPTH" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> $1_merge.summary.tmp
fi

paste <(samtools depth $1_merge.sorted.bam | awk '{if ($3 > '"10"') {print $0}}' | wc -l) <(fastalength /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> $1_merge.summary.tmp

paste <(samtools depth $1_merge.sorted.bam | awk '{if ($3 > '"100"') {print $0}}' | wc -l) <(fastalength /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> $1_merge.summary.tmp

paste <(samtools depth $1_merge.sorted.bam | awk '{if ($3 > '"1000"') {print $0}}' | wc -l) <(fastalength /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> $1_merge.summary.tmp

NCOUNT=$(seqtk comp $1_merge.depth10.fa | awk -F"\t" '{print $9}')

REVCOV=$(paste <(fastalength /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta | awk '{print $1}') <(seqtk comp $1_merge.depth10.fa | awk -F"\t" '{print $9}') | awk -F"\t" '{printf("%0.2f\n", ($1-$2)/$1*100)}')

NCOUNTPER=$(paste <(seqtk comp $1_merge.depth10.fa | awk -F"\t" '{print $9}') <(fastalength /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", ($1/$2)*100)}')

if [[ "$NCOUNT" == 0 ]]; then
    echo "0""#" | tr '#' '\t' | tr -d '\n' >> $1_merge.summary.tmp
else
    echo "$REVCOV" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> $1_merge.summary.tmp
fi

if [[ "$NCOUNT" == 0 ]]; then
    echo "29903""#" | tr '#' '\t' | tr -d '\n' >> $1_merge.summary.tmp
else
    echo "$NCOUNT" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> $1_merge.summary.tmp
fi

if [[ "$NCOUNT" == 0 ]]; then
    echo "100""#" | tr '#' '\t' | tr -d '\n' >> $1_merge.summary.tmp
else
    echo "$NCOUNTPER" | awk '{printf $0"#"}' | tr '#' '\t' | tr -d '\n' >> $1_merge.summary.tmp
fi

pangolin $1_merge.depth10.fa -t 12 --outfile $1_merge.lineage_report.csv

cat $1_merge.lineage_report.csv | sed -n 2p | awk -F, '{print $10"\t"$9"\t"$2}' | awk '{printf $0"#"}' | tr '#' '\t' >> $1_merge.summary.tmp

nextclade --version | awk '{printf $0"#"}' | tr '#' '\t' >> $1_merge.summary.tmp

nextclade -i $1_merge.depth10.fa -j 12 -t $1.nextclade.tsv

cat $1.nextclade.tsv | sed -n 2p | awk -F"\t" '{print $2"\t"$11"\t"$12"\t"$13"\t"$14"\t"$17"\t"$19}' | awk '{printf $0}' >> $1_merge.summary.tmp

sed '/^[[:space:]]*$/d' $1_merge.summary.tmp > $1_merge.summary.tsv

mafft --thread 12 --auto --quiet --keeplength --inputorder --6merpair --leavegappyregion --addfragments $1_merge.depth10.fa /home/lpmor22/IGM_SARSCOV2/ref_seq/MN908947.3.fasta | seqkit grep -vip MN908947.3 | sed '/>/!y/atcgn-/ATCGNN/' >> $1_merge.consensus.fasta

rm *.tmp