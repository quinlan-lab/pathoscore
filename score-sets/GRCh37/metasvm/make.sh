# 5th column
if [ ! -s MetaSVM.sorted.txt.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/MetaSVM.sorted.txt.gz
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/MetaSVM.sorted.txt.gz.tbi
fi
