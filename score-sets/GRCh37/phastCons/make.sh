# Use the 4th column (from vertebrate and from UCSC)
if [ ! -s phastcons_vertebrate.txt.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/phastcons/phastcons_vertebrate.txt.gz
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/phastcons/phastcons_vertebrate.txt.gz.tbi
fi
