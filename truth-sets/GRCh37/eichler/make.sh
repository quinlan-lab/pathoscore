if [ ! -s eichler.benign.vcf ] | [ ! -s eichler.pathogenic.vcf ]; then
    wget -c https://www.nature.com/nature/journal/v515/n7526/extref/nature13908-s2.zip
    unzip nature13908-s2.zip
    python make.py nature13908-s2/Supplementary\ Table\ 2.xlsx
fi

fasta=/uufs/chpc.utah.edu/common/home/u6000771/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
if [[ ! -e ../Homo_sapiens.GRCh37.82.gff3.gz ]]; then
    wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
    mv Homo_sapiens.GRCh37.82.gff3.gz ..
fi

gff=../Homo_sapiens.GRCh37.82.gff3.gz

sortit() {
    set -euo pipefail
    fn=$1
    grep "^#" $fn > tt.tmp.vcf
    grep -v "^#" $fn | sort -k1,1 -k2,2n >> tt.tmp.vcf
    bash ../../../scripts/bcsq.sh $gff tt.tmp.vcf $fasta | bgzip -c > $fn.gz # also annotates with bcftools
    rm $fn tt.tmp.vcf
}

sortit eichler.benign.vcf
tabix -f eichler.benign.vcf.gz
sortit eichler.pathogenic.vcf
tabix -f eichler.pathogenic.vcf.gz
