import gzip
import sys
import tabix
caddsnv = tabix.open(sys.argv[1])
caddindel = tabix.open(sys.argv[2])
caddprecomputed = gzip.open(sys.argv[3], 'rb')
ct = 0
for i, record in enumerate(caddprecomputed):
    if i == 0:
        print record,
        continue
    prefields = record.strip().split("\t")
    region = prefields[0]+":"+prefields[1]+"-"+prefields[1]
    for r in caddsnv.querys(region):
        if r[:4] == prefields[:4]:
            ct += 1
    for r in caddindel.querys(region):
        if r[:4] == prefields[:4]:
            ct += 1
    if ct == 0:
        print "\t".join(prefields[:-1]) + "\t" + "." + "\t" + prefields[-1]
    ct = 0
