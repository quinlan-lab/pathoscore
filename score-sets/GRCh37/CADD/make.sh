if [ ! -s whole_genome_SNVs.tsv.gz ] | [ ! -s InDels.tsv.gz ]; then
    wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz
    wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi
    wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels.tsv.gz
    wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels.tsv.gz.tbi
fi
