import sys
import toolshed as ts
from collections import defaultdict

genes = defaultdict(float)
for d in ts.reader(sys.argv[1]):
    genes[d['CCDSr20']] = float(d['%RVIS[pop_maf_0.05%(any)]']) # gene, value at 0.05% MAF
for line in open(sys.argv[2], "r"): # gene list
    fields = line.strip().split("\t") 
    chrom = fields[0]; start = fields[1]; end = fields[2]; gene = fields[3]
    RVIS=genes[gene]
    if RVIS == 0.0: # RVIS does not have this score, is dict default
        continue 
    print("{chrom}\t{start}\t{end}\t{gene}\t{RVIS}".format(chrom=chrom, start=start,
        end=end, gene=gene, RVIS=RVIS))
