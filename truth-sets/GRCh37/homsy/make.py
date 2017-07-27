from __future__ import print_function
import sys

import pyexcel as pe

path_book = pe.get_book(file_name=sys.argv[1])
benign_book = pe.get_book(file_name=sys.argv[2])

bfh = open('homsy.benign.vcf', 'w')
pfh = open('homsy.pathogenic.vcf', 'w')
header = """##fileformat=VCFv4.1
##source=pathoscore
##reference=GRCh37
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

print(header, file=bfh)
print(header, file=pfh)

sheet = path_book['Database S2']
for i, record in enumerate(sheet):
    if i == 0:
        continue
    if i == 1:
        keys = map(str, record)
        continue
    
    record = dict(zip(keys, record))
    if record['Variant Class'] == "Synonymous" or record['Variant Class'] == 'Splice site': continue
    print("{chrom}\t{pos}\t.\t{ref}\t{alt}\t10\tPASS\t.".format(\
            chrom=record['CHROM'], \
            pos=record['POS'], \
            ref=record['REF'], \
            alt=record['ALT']), \
            file=pfh)

sheet = benign_book['Database S3']
for i, record in enumerate(sheet):
    if i == 0:
        continue
    if i == 1:
        keys = map(str, record)
        continue
    
    record = dict(zip(keys, record))
    if record['Variant Class'] == "Synonymous" or record['Variant Class'] == 'Splice site': continue
    print("{chrom}\t{pos}\t.\t{ref}\t{alt}\t10\tPASS\t.".format(\
            chrom=record['CHROM'], \
            pos=record['POS'], \
            ref=record['REF'], \
            alt=record['ALT']), \
            file=bfh)

bfh.close()
pfh.close()
print("DONE")
