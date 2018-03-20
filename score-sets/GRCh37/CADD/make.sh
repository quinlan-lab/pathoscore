if [ ! -s whole_genome_SNVs.tsv.gz ]; then
    wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz
    wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi
fi

if [ ! -s InDels.tsv.gz ]; then
    wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels.tsv.gz
    wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels.tsv.gz.tbi
fi

if [ ! -s CADD_sorted.scores.txt.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/CADD_sorted.scores.txt.gz
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/CADD_sorted.scores.txt.gz.tbi
fi

# to ignore duplicate entries and prevent vcfanno errors

python filter.py whole_genome_SNVs.tsv.gz InDels.tsv.gz CADD_sorted.scores.txt.gz | uniq | bgzip -c > CADD_precomputed.txt.gz; tabix -b 2 -e 2 CADD_precomputed.txt.gz

#NOTE: use whole_genome_SNVs.tsv.gz, InDels.tsv.gz and CADD_precomputed.txt.gz, all column 6
