# column 4 is the score
if [ ! -s ccrs.autosomes.v2.20180420.bed.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz
    wget https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz.tbi
fi
# x chrom, same column
if [ ! -s ccrs.xchrom.v2.20180420.bed.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz
    wget https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz.tbi
fi
