set -euo pipefail

mkdir -p gerp
cd gerp
if [ ! -s hg19.GERP_scores.tar.gz ]; then
    wget http://mendel.stanford.edu/SidowLab/downloads/gerp/hg19.GERP_scores.tar.gz
    tar xzvf hg19.GERP_scores.tar.gz
fi

if [ ! -s ../gerp_rs.txt.gz ]; then
    ( echo -e "#chrom\tstart\tgerp_rs"
    for chr in $(seq 1 22) X Y; do
        >&2 echo $chr
        awk -vchr=$chr 'BEGIN{OFS="\t"}{ print chr,NR,$2 }' chr$chr.maf.rates
    done ) | bgzip -@ 2 -c > ../gerp_rs.txt.gz
    tabix -b2 -e2 ../gerp_rs.txt.gz
fi
