set -euo pipefail
wget http://science.sciencemag.org/highwire/filestream/646355/field_highwire_adjunct_files/1/aac9396_SupportingFile_Other_seq1_v4.zip
unzip aac9396_SupportingFile_Other_seq1_v4.zip

python make.py homsy_database_S02.xlsx homsy_database_S03.xlsx

sortit() {
    fn=$1
    grep "^#" $fn > tt.tmp
    grep -v "^#" $fn | sort -k1,1 -k2,2n >> tt.tmp
    bgzip -c tt.tmp > $fn.gz
    rm $fn tt.tmp
}

sortit homsy.benign.vcf
tabix -f homsy.benign.vcf.gz
sortit homsy.pathogenic.vcf
tabix -f homsy.pathogenic.vcf.gz
