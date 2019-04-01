from __future__ import print_function
from cyvcf2 import VCF
import sys
from collections import Counter

vcf = VCF(sys.argv[1]) # clinvar/samocha files go here

"""
Simply grabs the first gene entry for each BCSQ info field.  Genes will often repeat multiple times in the same BCSQ string, and though they may be different in rare cases, this will give a nonetheless accurate count of the 10 most frequent genes in each variant set.
"""

c = Counter()
for variant in vcf:
    gene=str(variant.INFO.get("BCSQ").split(",")[0].split("|")[1])
    c.update({gene: 1})
topten = c.most_common(10)
out=sys.argv[1].split(".")[0]+"genes.txt"
f=open(out, 'w')
for i in topten:
    print (i[0], file=f)
f.close()
