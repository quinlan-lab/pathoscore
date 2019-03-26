# 12th column
if [ ! -s GRCh37.74.zip ]; then
    wget http://sift.bii.a-star.edu.sg/sift4g/public//Homo_sapiens/GRCh37.74.zip
    unzip GRCh37.74.zip
    rm GRCh37.74.zip
fi
cd GRCh37.74
rm *.regions *_stats.txt HG* GL*
for file in *gz; do
    filename="${file%.*}"
    gunzip $file -c > $filename
    sed -i '1d' $filename
    sed -i -e "1,$ s/^/$filename\t/" $filename
    cat $filename >> sift.txt
done
mv sift.txt ..
cd -
cat <(printf "#Chrom\tPosition\tRef_allele\tNew_allele\tTranscript_id\tGene_id\tGene_name\tRegion\tRef_amino_acid\tNew_amino_acid\tPosition_of_amino_acid_substitution\tSIFT_score\tSIFT_median_sequence_info\tNum_seqs_at_position\tdbSNP_id\n") sift.txt | bgzip -c > sift.txt.gz; tabix -b 2 -e 2 sift.txt.gz
