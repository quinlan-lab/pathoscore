date=20170710
set -e
if [[ ! -e clinvar_${date}.vcf.gz ]]; then
wget ftp://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/vcf_2.0/clinvar_${date}.vcf.gz
fi
tabix -f clinvar_${date}.vcf.gz
python make.py clinvar_${date}.vcf.gz $date
for f in clinvar-*.vcf; do
	bgzip -f $f
	tabix -f $f.gz
done
