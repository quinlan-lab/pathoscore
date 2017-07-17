set -euo pipefail
gff=$1
vcf=$2
fasta=$3
bcftools csq -f $fasta -g $gff --local-csq --samples - $vcf
