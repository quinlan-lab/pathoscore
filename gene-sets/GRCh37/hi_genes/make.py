from __future__ import print_function
import sys

import pyexcel as pe

book = pe.get_book(file_name=sys.argv[1])

genelist = open('dang_hi.tsv', 'w')

sheet = book['hap_with_blocks']
for i, record in enumerate(sheet):
    if i == 0:
        keys = map(str, record)
        continue
    record = dict(zip(keys, record))
    print("{gene}".format(gene=record['Gene Symbol']), file=genelist)

genelist.close()
print("DONE")
