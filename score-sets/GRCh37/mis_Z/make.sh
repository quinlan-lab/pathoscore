set -euo pipefail

name=nature19057-s2.zip # latest, most up-to-date release
if [[ ! -e "$name" ]]; then
    wget https://media.nature.com/original/nature-assets/nature/journal/v536/n7616/extref/nature19057-s2.zip
unzip nature19057-s2.zip
fi
python make.py nature19057-s2/nature19057-SI\ Table\ 13.xlsx
sort -k1,1 -k2,2n missensez.bed | bgzip -c > missensez.bed.gz; tabix missensez.bed.gz # value in 5th column
