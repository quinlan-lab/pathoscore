set -euo pipefail

#Dang et al., manually curated HI set of genes
if [ ! -s ejhg2008111x1.xls ]; then
    wget https://media.nature.com/original/nature-assets/ejhg/journal/v16/n11/extref/ejhg2008111x1.xls
fi

python make.py ejhg2008111x1.xls

if [[ ! -s ../Homo_sapiens.GRCh37.75.gtf.gz ]]; then
    wget -P .. ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz # final release for GRCh37
#only grabs protein coding sequence, CDS or stop_codon (separate from CDS for some reason) | then only grabs the autosomal and sex chromosomes | then turns 1-based GFF format into 0-based half-open, BED format | cuts out necessary columns, sorts by gene

    zgrep 'protein_coding\tstop_codon\|protein_coding\tCDS' ../Homo_sapiens.GRCh37.75.gtf.gz | grep "^1\|^2\|^3\|^4\|^5\|^6\|^7\|^8\|^9\|^10\|^11\|^12\|^13\|^14\|^15\|^16\|^17\|^18\|^19\|^20\|^21\|^22\|^X\|^Y" | awk '{$4=$4-1; print $0}' OFS='\t' | cut -f 1,4,5,12,16 > ../codingtranscriptome.bed
# necessary edit for easy awk command
    sed 's/\"//g' ../codingtranscriptome.bed | sed 's/;//g' > ../Homo_sapiens37.bed
fi
# matches gene names to bed coordinates
awk 'FNR==NR{genes[$1]; next} {for (gene in genes) if (gene == $5) print $0, genes[gene]}' FS='\t' OFS='\t' dang_hi.tsv ../Homo_sapiens37.bed | sort -k1,1 -k2,2n | cut -f -3 > unflattened_hi_genes.bed
# creates flattened representation of protein-coding exome covering AD genes
bedtools merge -i unflattened_hi_genes.bed > hi_genes.bed

sort -k1,1 -k2,2n hi_genes.bed | bgzip -c > hi_genes.bed.gz
tabix hi_genes.bed.gz

# bedtools complement so we can use the EXCLUDE option
if [[ ! -s hg19.chrom.sizes ]]; then
    wget https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
fi
bedtools complement -i <(cut -f -3 hi_genes.bed | sort -k1,1 -k2,2n) -g <(sed 's/^chr//g' hg19.chrom.sizes | sort -k1,1) > hi_gene_complement.bed
# sorted, bgzipped, and tabixed
sort -k1,1 -k2,2n hi_gene_complement.bed | bgzip -c > hi_gene_complement.bed.gz
tabix -f hi_gene_complement.bed.gz

# ClinGen HI set from MacArthur lab
if [ ! -s clingen_level3_genes_2015_02_27.tsv ]; then
    wget https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/clingen_level3_genes_2015_02_27.tsv
fi

awk 'FNR==NR{genes[$1]; next} {for (gene in genes) if (gene == $5) print $0, genes[gene]}' FS='\t' OFS='\t' clingen_level3_genes_2015_02_27.tsv ../Homo_sapiens37.bed | sort -k1,1 -k2,2n | cut -f -3 > unflattened_clingen_hi_genes.bed
# creates flattened representation of protein-coding exome covering AD genes
bedtools merge -i unflattened_clingen_hi_genes.bed > clingen_hi_genes.bed

sort -k1,1 -k2,2n clingen_hi_genes.bed | bgzip -c > clingen_hi_genes.bed.gz
tabix clingen_hi_genes.bed.gz

bedtools complement -i <(cut -f -3 clingen_hi_genes.bed | sort -k1,1 -k2,2n) -g <(sed 's/^chr//g' hg19.chrom.sizes | sort -k1,1) > clingen_hi_gene_complement.bed
# sorted, bgzipped, and tabixed
sort -k1,1 -k2,2n clingen_hi_gene_complement.bed | bgzip -c > clingen_hi_gene_complement.bed.gz
tabix -f clingen_hi_gene_complement.bed.gz
