set -euo pipefail

name=RVIS_Unpublished_ExACv2_March2017.txt # latest, most up-to-date release
if [[ ! -e "$name" ]]; then
    wget http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt
fi
if [[ ! -e ../Homo_sapiens.GRCh37.82.gff3.gz ]]; then
        wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
        mv Homo_sapiens.GRCh37.82.gff3.gz ..
fi
#creates flattened transcriptome
zgrep -P 'ensembl\tgene|ensembl_havana\tgene' ../Homo_sapiens.GRCh37.82.gff3.gz | grep -v "^#" | awk '{$4=$4-1; print $0}' OFS='\t' | tr -s ";" "\t" | cut -f 1,4,5,10 | tr -d "Name=" | sort -k1,1 -k2,2n > gene.bed
python make.py RVIS_Unpublished_ExACv2_March2017.txt gene.bed | bgzip -c > rvis.bed.gz; tabix rvis.bed.gz # use 5th column
