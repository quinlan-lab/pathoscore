set -euo pipefail
mkdir -p fitcons
cd fitcons

wget http://compgen.cshl.edu/fitCons/0downloads/tracks/i6/scores/fc-i6-0.bw
wget http://compgen.cshl.edu/fitCons/0downloads/tracks/i6/scores/fc-i6-1.bw
wget http://compgen.cshl.edu/fitCons/0downloads/tracks/i6/scores/fc-i6-2.bw
wget http://compgen.cshl.edu/fitCons/0downloads/tracks/i6/scores/fc-i6-3.bw

for i in 0 1 2 3; do
    bigWigToBedGraph fc-i6-$i.bw fitcons.$i.bed
    sed -s 's/^chr//' fitcons.$i.bed
    rm fitcons.$i.bed
done | sort -k1,1 -k2,2n | bgzip -@ 4 -c > fitcons.bed.gz
tabix fitcons.bed.gz

