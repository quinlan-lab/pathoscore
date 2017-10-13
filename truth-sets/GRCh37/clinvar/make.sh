date=20170905
set -euo pipefail
if [[ ! -e clinvar_${date}.vcf.gz ]]; then
    wget ftp://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2017/clinvar_${date}.vcf.gz
fi

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

tabix -f clinvar_${date}.vcf.gz
python make.py clinvar_${date}.vcf.gz $date
for vcf in clinvar-*.vcf; do
	bash ../../../scripts/bcsq.sh $gff $vcf $fasta \
		| python ../../score.py - \
		| bgzip -c > $vcf.gz
	rm $vcf
	tabix -f $vcf.gz
done
