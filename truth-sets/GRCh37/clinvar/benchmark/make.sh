# clinvar ALL uses random variants from the whole gnomAD set to add to clinvar benigns
if [ ! -s clinvar-benign_randomALL.vcf.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/clinvar-benchmark/clinvar-benign_randomALL.vcf.gz
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/clinvar-benchmark/clinvar-benign_randomALL.vcf.gz.tbi
fi
# clinvar singleton uses random variants from the whole gnomAD set to add to clinvar benigns
if [ ! -s clinvar-benign_randomsingleton.vcf.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/clinvar-benchmark/clinvar-benign_randomsingleton.vcf.gz
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/clinvar-benchmark/clinvar-benign_randomsingleton.vcf.gz.tbi
fi
# clinvar AF < 0.005 uses random variants from the gnomAD set < 0.005 AF to add to clinvar benigns
if [ ! -s clinvar-benign_randomAF0.005.vcf.gz ]; then
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/clinvar-benchmark/clinvar-benign_randomAF0.005.vcf.gz
    wget https://s3.us-east-2.amazonaws.com/pathoscore-data/clinvar-benchmark/clinvar-benign_randomAF0.005.vcf.gz.tbi
fi
