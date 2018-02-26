#!/bin/bash

# can be used as is with pathoscore.py with --scores score-sets/GRCh37/DANN/DANN_whole_genome_SNVs.tsv.bgz:DANN:5:self
dann_file="DANN_whole_genome_SNVs.tsv.bgz"
if [ ! -f $dann_file ]; then
	wget https://cbcl.ics.uci.edu/public_data/DANN/data/DANN_whole_genome_SNVs.tsv.bgz
fi

dann_index="DANN_whole_genome_SNVs.tsv.bgz.tbi"
if [ ! -f $dann_index ]; then
	wget https://cbcl.ics.uci.edu/public_data/DANN/data/DANN_whole_genome_SNVs.tsv.bgz.tbi
fi

