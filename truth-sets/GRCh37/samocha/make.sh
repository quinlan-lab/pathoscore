set -euo pipefail
fname=148353-6.xlsx
if [[ ! -e $fname ]]; then
wget http://www.biorxiv.org/highwire/filestream/44323/field_highwire_adjunct_files/5/$fname
fi
python make.py $fname

gff=../Homo_sapiens.GRCh37.82.gff3.gz

if [[ ! -s $gff ]]; then
    wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
    mv Homo_sapiens.GRCh37.82.gff3.gz ..
fi

fasta=/uufs/chpc.utah.edu/common/home/u6000771/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa

if [ ! -s $fasta ]; then
    fasta=/data/human/g1k_v37_decoy.fa
fi

if [ ! -s $fasta ]; then
    fasta=../human_g1k_v37.fasta
fi

if [ ! -s $fasta ]; then
    wget -P .. ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
    gunzip ../human_g1k_v37.fasta.gz
fi

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
