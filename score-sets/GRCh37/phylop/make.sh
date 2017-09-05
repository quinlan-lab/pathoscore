rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals ./
cd placentalMammals
rm chrUn*.gz
rm chr*hap*.gz
rm chr*random*.gz
gunzip *gz
set -euo pipefail
for f in `ls -v ./chr*wigFix`; do
    /data2/qli/tool/bigwig/bin/convert2bed -i wig < $f > $f.bed 
    perl -lane 'next if @F<5;$F[0]=~ s/^chr//;print "$F[0]\t$F[1]\t$F[2]\t$F[4]";' $f.bed | sort --temporary-directory=. -k1,1 -k2,2n
    rm $f.bed 
done  | bgzip -c > phylop.bed.gz
tabix -b2 -e2 phylop.bed.gz
