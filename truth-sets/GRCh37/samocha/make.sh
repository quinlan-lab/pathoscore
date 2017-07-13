set -euo pipefail
fname=148353-6.xlsx
if [[ ! -e $fname ]]; then
wget http://www.biorxiv.org/highwire/filestream/44323/field_highwire_adjunct_files/5/$fname
fi
python make.py $fname

sortit() {
    fn=$1
    grep "^#" $fn > tt.tmp
    grep -v "^#" $fn | sort -k1,1 -k2,2n >> tt.tmp
    bgzip -c tt.tmp > $fn.gz
    rm $fn tt.tmp
}

sortit samocha.benign.vcf
tabix -f samocha.benign.vcf.gz
sortit samocha.pathogenic.vcf
tabix -f samocha.pathogenic.vcf.gz
