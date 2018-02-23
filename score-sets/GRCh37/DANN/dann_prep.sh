#!/bin/bash

mkdir -p DANN
cd DANN
dann_file="DANN_whole_genome_SNVs.tsv.bgz"
if [ ! -f $dann_file ]; then
	wget https://cbcl.ics.uci.edu/public_data/DANN/data/DANN_whole_genome_SNVs.tsv.bgz*
fi

echo "##INFO=<ID=DANN_SCORE, Number=1,Type=Float>" > DANN_whole_genome_SNVs.vcf
echo -e "#CHROM\tPOS\tREF\tALT\tINFO" >> DANN_whole_genome_SNVs.vcf
zcat DANN_whole_genome_SNVs.tsv.bgz | awk -v OFS="" '{print $1,"\t", $2, "\t", $3, "\t",  $4, "\tDANN_SCORE=", $5}' >> DANN_whole_genome_SNVs.vcf
bgzip DANN_whole_genome_SNVs.vcf
tabix -p vcf DANN_whole_genome_SNVs.vcf.gz
