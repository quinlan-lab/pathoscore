# 5th column
if [ ! -s SIFT_sorted.scores.txt.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/SIFT_sorted.scores.txt.gz
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/precomputed-metrics/SIFT_sorted.scores.txt.gz.tbi
fi
