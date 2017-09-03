from __future__ import print_function
import gzip

chroms = sorted(map(str,range(1,23)+['X','Y']))
out = open('aloft.txt', 'w')
print("#chrom\tpos\tref\talt\tDominant\tLoftol\tRecessive", file = out)
inds=[0,1,3,4,11,12,13]
for i in chroms:
    with gzip.open('chr'+i+'.vcf.vat.aloft.lof.predict.gz', 'rb') as f:
        ct = 0
        for line in f:
            if ct == 0:
                ct += 1
                continue 
            print("\t".join([line.strip().split("chr")[1].split("\t")[i].strip() for i in inds]), file = out)
