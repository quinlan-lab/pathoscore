#mkdir -p fathmm
#cd fathmm
#wget  http://fathmm.biocompute.org.uk/database/fathmm-MKL_Current_zerobased.tab.gz 
#gunzip *gz
set -euo pipefail
(echo -e "#chrom\tpos\tref\talt\tnonCodingScore\tcodingScore" \
| perl -lane 'print "$F[0]\t$F[1]\t$F[3]\t$F[4]\t$F[5]\t$F[7]"' fathmm-MKL_Current_zerobased.tab \
| sort -k1,1 -k2,2n) | bgzip -c > fathmm.txt.gz
tabix -b2 -e2 fathmm.txt.gz
