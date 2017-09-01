mkdir -p linsight
cd linsight
wget  http://compgen.cshl.edu/%7Eyihuang/tracks/LINSIGHT.bw
set -euo pipefail
bigWigToBedGraph LINSIGHT.bw linsight.bedg
sed -s 's/^chr//' linsight.bedg \
| sort -k1,1 -k2,2n | bgzip -c > linsight.bed.gz
tabix -p bed linsight.bed.gz
rm linsight.bedg
