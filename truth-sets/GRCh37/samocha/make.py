from __future__ import print_function
import sys

import pyexcel as pe

book = pe.get_book(file_name=sys.argv[1])

bfh = open('samocha.benign.vcf', 'w')
pfh = open('samocha.pathogenic.vcf', 'w')
header = """##fileformat=VCFv4.1
##source=pathoscore
##reference=GRCh37
##INFO=<ID=MPC,Number=1,Type=Float,Description="MPC score">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

print(header, file=bfh)
print(header, file=pfh)

sheet = book['Table_S8']
for i, record in enumerate(sheet):
    if i == 0:
        keys = map(str, record)
        continue
    record = dict(zip(keys, record))
    fh = bfh if (record['dataset'] == 'control') else pfh
    if record['MPC'] == "NA":
        print("{chrom}\t{pos}\t.\t{ref}\t{alt}\t10\tPASS\t.".format(chrom=record['chrom'], pos=record['pos'],
          ref=record['ref'], alt=record['alt']), file=fh)
    else:
        print("{chrom}\t{pos}\t.\t{ref}\t{alt}\t10\tPASS\tMPC={MPC}".format(chrom=record['chrom'], pos=record['pos'],
          ref=record['ref'], alt=record['alt'], MPC="%.4f" % float(record['MPC'])), file=fh)


bfh.close()
pfh.close()
print("DONE")
