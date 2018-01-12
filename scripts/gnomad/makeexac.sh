original=ExAC.r1.sites.vep.vcf.gz
fileName="${original##*/}"
fileExt=${fileName#*.}
FILE=${fileName%*.$fileExt}
if [ ! -s $original ]; then
    wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz
    wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz.tbi
fi

fasta=../../truth-sets/GRCh37/human_g1k_v37.fasta

if [ ! -s $fasta ]; then
    wget -P ../../truth-sets/GRCh37/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
    gunzip ../../truth-sets/GRCh37/human_g1k_v37.fasta.gz
fi

gff=../../truth-sets/GRCh37/Homo_sapiens.GRCh37.82.gff3.gz

if [[ ! -s $gff ]]; then
    wget -P ../../truth-sets/GRCh37/ ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
fi

if [ ! -s ${FILE}_dec.vcf ]; then
    vt decompose $original -o ${FILE}_dec.vcf -s
fi

if [ ! -s ${FILE}.vcf ]; then
    vt normalize ${FILE}_dec.vcf -o ${FILE}.vcf -r $fasta
fi

sortit() {
    set -euo pipefail
    fn=$1
    grep "^#" $fn > tt.tmp.vcf
    grep -v "^#" $fn | sort -k1,1 -k2,2n -T . >> tt.tmp.vcf
	bash ../bcsq.sh $gff tt.tmp.vcf $fasta \
        | python ../../truth-sets/score.py - \
        | bgzip -c > $fn.gz

    rm $fn tt.tmp.vcf
}

sortit ExAC.vcf
tabix -f ExAC.vcf.gz
