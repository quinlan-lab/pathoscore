if [ ! -s wigler-autism.benign.vcf ] | [ ! -s wigler-autism.pathogenic.vcf ]; then
    wget -c https://www.nature.com/nature/journal/v515/n7526/extref/nature13908-s2.zip
    unzip nature13908-s2.zip
    python make.py nature13908-s2/Supplementary\ Table\ 2.xlsx
fi

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
		| bgzip -c > t.vcf.gz
	tabix t.vcf.gz
    python ../../score.py t.vcf.gz \
        | bgzip -c > $fn.gz

    rm $fn tt.tmp.vcf t.vcf.gz*
}

sortit wigler-autism.benign.vcf
tabix -f wigler-autism.benign.vcf.gz
sortit wigler-autism.pathogenic.vcf
tabix -f wigler-autism.pathogenic.vcf.gz
