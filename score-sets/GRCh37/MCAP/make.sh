if [ ! -s mcap.txt.gz ]; then
    wget http://bejerano.stanford.edu/mcap/downloads/dat/mcap_v1_0.txt.gz -O mcap.txt.gz
    gunzip mcap.txt.gz
    sed 's/^#grch37_/#/g' mcap.txt | bgzip -c > mcap.txt.gz
    tabix mcap.txt.gz -b 2 -e 2
    rm mcap.txt
fi
