set -euo pipefail
fname=148353-6.xlsx
if [[ ! -e $fname ]]; then
wget http://www.biorxiv.org/highwire/filestream/44323/field_highwire_adjunct_files/5/$fname
fi
python make.py $fname

fasta=/uufs/chpc.utah.edu/common/home/u6000771/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
fasta=/data/human/g1k_v37_decoy.fa
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
	bash ../../../scripts/bcsq.sh $gff tt.tmp.vcf $fasta \
        | python ../../score.py - \
        | bgzip -c > $fn.gz

    rm $fn tt.tmp.vcf
}

sortit samocha.benign.vcf
tabix -f samocha.benign.vcf.gz
sortit samocha.pathogenic.vcf
tabix -f samocha.pathogenic.vcf.gz
