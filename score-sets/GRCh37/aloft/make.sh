if [ ! -s aloft.txt.gz ]; then
    wget http://files.gersteinlab.org/aloft/hg19/chr{{1..22},{X..Y}}.vcf.vat.aloft.lof.predict.gz
    python make.py
    bgzip aloft.txt
    tabix aloft.txt.gz -b 2 -e 2
    rm *predict.gz
fi
