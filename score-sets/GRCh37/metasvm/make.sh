# 5th column
if ! compgen -G "dbNSFP2.9.3*" > /dev/null ; then # last version on hg19, 3 and above are hg38
    fileid="0B60wROKy6OqceTNZRkZnaERWREk"
    filename="dbNSFP2.9.3.zip"
    curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
    curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
    rm cookie
    unzip dbNSFP2.9.zip
    rm dbNSFP2.9.zip
    rm try*; rm search*
fi
rm metasvm.txt
for chr in $(seq 1 22) X Y; do
    python make.py dbNSFP2.9.3_variant.chr$chr >> metasvm.txt
done 
sort -k1,1 -k2,2n metasvm.txt | bgzip -c > metasvm.txt.gz; tabix -b 2 -e 2 metasvm.txt.gz
