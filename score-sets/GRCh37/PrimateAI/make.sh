# use 11th column
if [ ! -s PrimateAI_scores_v0.2.tsv.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/primateai/PrimateAI_scores_v0.2.tsv.gz
fi
zgrep -v "^#" PrimateAI_scores_v0.2.tsv.gz | sed '1d' | sed '1s/^/#/' | sort -k1,1 -k2,2n > PrimateAI_scores_v0.2.tsv
bgzip -f PrimateAI_scores_v0.2.tsv
tabix -b 2 -e 2 PrimateAI_scores_v0.2.tsv.gz
