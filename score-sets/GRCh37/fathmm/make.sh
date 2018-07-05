set -euo pipefail

mkdir -p fathmm
cd fathmm
if [ ! -s fathmm-MKL_Current_zerobased.tab.gz ]; then
    wget http://fathmm.biocompute.org.uk/database/fathmm-MKL_Current_zerobased.tab.gz 
fi
gunzip *gz
if [ ! -s fathmm.txt.gz ]; then
    (echo -e "#chrom\tpos\tref\talt\tnonCodingScore\tcodingScore";\
    perl -lane 'print "$F[0]\t$F[1]\t$F[3]\t$F[4]\t$F[5]\t$F[7]"' fathmm-MKL_Current_zerobased.tab) \
    | bgzip -c > fathmm.txt.gz
    tabix -b2 -e2 fathmm.txt.gz
fi
