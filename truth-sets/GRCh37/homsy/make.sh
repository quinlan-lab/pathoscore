set -euo pipefail

if [[ ! -e aac9396_SupportingFile_Other_seq1_v4.zip ]]; then
wget http://science.sciencemag.org/highwire/filestream/646355/field_highwire_adjunct_files/1/aac9396_SupportingFile_Other_seq1_v4.zip
unzip aac9396_SupportingFile_Other_seq1_v4.zip
fi

python make.py homsy_database_S02.xlsx homsy_database_S03.xlsx

fasta=/uufs/chpc.utah.edu/common/home/u6000771/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
if [ ! -s $fasta ]; then
    fasta=../human_g1k_v37.fasta
fi
if [ ! -s $fasta ]; then
    wget -P .. ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
    gunzip ../human_g1k_v37.fasta.gz
fi
if [[ ! -e ../Homo_sapiens.GRCh37.82.gff3.gz ]]; then
    wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
    mv Homo_sapiens.GRCh37.82.gff3.gz ..
fi

gff=../Homo_sapiens.GRCh37.82.gff3.gz

sortit() {
    set -euo pipefail
    fn=$1
    grep "^#" $fn > tt.tmp.vcf
    grep -v "^#" $fn | sort -k1,1 -k2,2n | grep -vw 46279769 >> tt.tmp.vcf
	bash ../../../scripts/bcsq.sh $gff tt.tmp.vcf $fasta \
		| bgzip -c > t.vcf.gz
	tabix t.vcf.gz
    python ../../score.py t.vcf.gz \
        | bgzip -c > $fn.gz

    rm $fn tt.tmp.vcf t.vcf.gz*
}

sortit homsy.benign.vcf
tabix -f homsy.benign.vcf.gz
sortit homsy.pathogenic.vcf
tabix -f homsy.pathogenic.vcf.gz
