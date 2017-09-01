mkdir -p sift
cd sift
wget  ftp://ftp.jcvi.org/pub/data/sift/Human_db_37_ensembl_63/Human*sqlite.gz 
gunzip *gz
set -euo pipefail
(echo -e "#chrom\tpos\tref\talt\tsiftScore";
for f in `ls -v ./Human_CHR*.sqlite`; do
    set -e
    perl query.pl $f | sort --temporary-directory=. -k1,1 -k2,2n
done ) | bgzip -c > sift.txt.gz
tabix -b2 -e2 sift.txt.gz
rm Human*sqlite.gz
