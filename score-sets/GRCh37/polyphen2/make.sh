#HumDiv version (5th column)

if [ ! -s Polyphen2_HDIV.sorted.txt.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/Polyphen2_HDIV.sorted.txt.gz
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/Polyphen2_HDIV.sorted.txt.gz.tbi
fi

#HumVar version (5th column)

if [ ! -s Polyphen2_HVAR.sorted.txt.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/Polyphen2_HVAR.sorted.txt.gz
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/Polyphen2_HVAR.sorted.txt.gz.tbi
fi
