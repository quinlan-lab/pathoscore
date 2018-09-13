from __future__ import print_function
import sys

import pyexcel as pe

book = pe.get_book(file_name=sys.argv[1])

genelist = open('eiee_genes.tsv', 'w')

sheet = book['Supp Table 4']
for i, record in enumerate(sheet):
    if i < 14: # list is in first column but has a bunch of information above it, list doesn't start till row 15
        continue
    print("{gene}".format(gene=record[0]), file=genelist)

genelist.close()
print("DONE")
