if [ ! -s revel.txt.gz ]; then
    wget https://rothsj06.u.hpc.mssm.edu/revel/revel_all_chromosomes.csv.zip
    unzip revel_all_chromosomes.csv.zip
    tr -s ',' '\t' < revel_all_chromosomes.csv | sed '1 s/^/#/' | sed '1 s/hg19_//g' | bgzip -c > revel.txt.gz
    tabix revel.txt.gz -b 2 -e 2
    rm revel_all_chromosomes.csv.zip revel_all_chromosomes.csv
fi
