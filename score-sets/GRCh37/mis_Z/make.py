from __future__ import print_function
import sys
import pyexcel as pe

book = pe.get_book(file_name=sys.argv[1])

fh = open('missensez.bed', 'w')

sheet = book['Gene Constraint']
for i, record in enumerate(sheet):
    if i == 0:
        keys = map(str, record)
        continue
    record = dict(zip(keys, record))
    print("{chrom}\t{start}\t{end}\t{gene}\t{mis_z}".format(chrom=record['chr'], start=record['cds_start'],
          end=record['cds_end'], gene=record['gene'], mis_z="%.4f" % float(record['mis_z'])), file=fh)
fh.close()
print("DONE")
