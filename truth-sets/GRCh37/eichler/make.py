from __future__ import print_function
import sys

import pyexcel as pe

path_book = pe.get_book(file_name=sys.argv[1])

bfh = open('eichler.benign.vcf', 'w')
pfh = open('eichler.pathogenic.vcf', 'w')
header = """##fileformat=VCFv4.1
##source=pathoscore
##reference=GRCh37
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

print(header, file=bfh)
print(header, file=pfh)

sheet = path_book['Supplement-T2-eventsTable']
for i, record in enumerate(sheet):
    #if i == 0:
    #    continue
    if i == 0:
        keys = map(str, record)
        continue
    
    record = dict(zip(keys, record))
    chrom, pos, ref, alt = record['vcfVariant'].split(":")
    if not any(["strong" in i or "valid" in i for i in [record['CSHL'] + record['YALE'] + record['UW']]]): continue #if strong GT confidence in one pipeline or validated, keep
    if record['familyDescription'] in ['pM', 'pF']: # if in probands only, pathogenic
        print("{chrom}\t{pos}\t.\t{ref}\t{alt}\t10\tPASS\t.".format(\
                chrom=chrom, \
                pos=pos, \
                ref=ref, \
                alt=alt), \
                file=pfh)
    else: #if not in proband only, probably not pathogenic...
        print("{chrom}\t{pos}\t.\t{ref}\t{alt}\t10\tPASS\t.".format(\
                chrom=chrom, \
                pos=pos, \
                ref=ref, \
                alt=alt), \
                file=bfh)

bfh.close()
pfh.close()
print("DONE")
