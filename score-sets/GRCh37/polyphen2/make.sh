wget ftp://genetics.bwh.harvard.edu/pph2/whess/polyphen-2.2.2-whess-2011_12.tab.tar.bz2
mkdir -p polyphen2
cd polyphen2
tar xjvf ../polyphen-2.2.2-whess-2011_12.tab.tar.bz2
set -euo pipefail
(echo -e "#chrom\tpos\tref\talt\thdiv\thvar";

for f in polyphen-2.2.2-whess-2011_12/*features.tab; do
    set -e
    s=$(dirname $f)/$(basename $f .features.tab).scores.tab
    paste $f $s | cut -f 1,2,58,64 
done | grep -v ^# \
    | perl -pe 's/:|\//\t/g' \
    | sed 's/^chr//; s/ //g' \
    | awk '$4 != ""' \
    | sort -k1,1 -k2,2n) | bgzip -c > polyphen2.txt.gz

tabix -b2 -e2 polyphen2.txt.gz
