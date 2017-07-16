date=2017_07_14
set -e
# grabs files
if [[ ! -e berg_ad.tsv ]]; then
    wget https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/berg_ad.tsv
fi
if [[ ! -s Homo_sapiens.GRCh37.75.gtf.gz ]]; then
    wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz # final release for GRCh37
fi
#only grabs protein coding sequence, CDS or stop_codon (separate from CDS for some reason) | then only grabs the autosomal and sex chromosomes | then turns 1-based GFF format into 0-based half-open, BED format | cuts out necessary columns, sorts by gene
zgrep -P 'protein_coding\tstop_codon|protein_coding\tCDS' Homo_sapiens.GRCh37.75.gtf.gz | grep -P "^1|^2|^3|^4|^5|^6|^7|^8|^9|^10|^11|^12|^13|^14|^15|^16|^17|^18|^19|^20|^21|^22|^X|^Y" | awk '{$4=$4-1; print $0}' OFS='\t' | cut -f 1,4,5,12,16 > codingtranscriptome.bed
# necessary edit for easy awk command
sed 's/\"//g' codingtranscriptome.bed | sed 's/;//g' > Homo_sapiens37.bed
# matches gene names to bed coordinates
awk 'FNR==NR{genes[$1]; next} {for (gene in genes) if (gene == $5) print $0, genes[gene]}' FS='\t' OFS='\t' berg_ad.tsv Homo_sapiens37.bed | sort -k1,1 -k2,2n | cut -f -3 > unflattened_ad_genes.bed
# creates flattened representation of protein-coding exome covering AD genes
bedtools merge -i unflattened_ad_genes.bed > ad_genes.bed

# bedtools complement so we can use the EXCLUDE option
if [[ ! -s hg19.chrom.sizes ]]; then
    wget https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
fi
bedtools complement -i <(cut -f -3 ad_genes.bed | sort -k1,1 -k2,2n) -g <(sed 's/^chr//g' hg19.chrom.sizes | sort -k1,1) > ad_gene_complement.bed
# sorted, bgzipped, and tabixed
sort -k1,1 -k2,2n ad_gene_complement.bed | bgzip -c > ad_gene_complement.bed.gz
tabix -f ad_gene_complement.bed.gz
