wget -nd -r --accept "*illumina.vcf.gz*"  https://genomics.scripps.edu/browser/files/wellderly/vcf/ --no-check-certificate

wget wget https://s3.amazonaws.com/gemini-annotations/gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz
wget wget https://s3.amazonaws.com/gemini-annotations/gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz.tbi

wget -O vcfanno https://github.com/brentp/vcfanno/releases/download/v0.2.8/vcfanno_linux64
chmod +x vcfanno

wget -O gargs https://github.com/brentp/gargs/releases/download/v0.3.8/gargs_linux
chmod +x gargs

cat > conf.toml << EOL
[[annotation]]
file="gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz"
fields = ["AC"]
ops=["flag"]
names=["gnomad_ac"]
EOL

ls chr*.illumina.vcf.gz \
    | cut -d"." -f1 \
    | ./gargs -p 24 'vcfanno conf.toml <(zcat {}.illumina.vcf.gz | grep -v "^####INFO" | cut -f1-7 | bcftools view) | bcftools view -Oz -i "gnomad_ac=0" > {}.vcf.gz'

ls chr*.vcf.gz \
    | grep -v illumina \
    | ./gargs -p 24 "tabix {}"

wget -O gsort https://github.com/brentp/gsort/releases/download/v0.0.6/gsort_linux_amd64
chmod +x gsort

bcftools view -h gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz \
    | grep contig \
    | sed -e "s/##contig=<ID=//" \
    | sed -e "s/>//" \
    | sed -e "s/,length=/\t/" \
    > genome.txt

bcftools concat $( ls chr*.vcf.gz | grep -v illumina ) \
    | ./gsort /dev/stdin genome.txt \
    | bcftools view -Oz > wellderly.vcf.gz

wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
bedtools intersect -u -a wellderly.vcf.gz -b Homo_sapiens.GRCh37.82.gff3.gz -header | bcftools view -Oz > wellderly.coding.vcf.gz
