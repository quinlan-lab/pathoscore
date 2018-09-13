from __future__ import print_function
import sys

import pyexcel as pe

book = pe.get_book(file_name=sys.argv[1])

genelist = open('heart_genes.tsv', 'w')

sheet = book['Suppl Table 20']
for i, record in enumerate(sheet):
    if i == 0:
        continue
    print("{gene}".format(gene=record[0]), file=genelist)

genelist.close()
print("DONE")
