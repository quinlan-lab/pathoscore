import sys
import toolshed as ts

vals = []
for d in ts.reader(sys.argv[1]):
    vals.append((d['chr'], int(d['tx_start']), int(d['tx_end']), d['gene'], d['pLI']))
vals.sort()
for v in vals:
    print("{chrom}\t{start}\t{end}\t{gene}\t{pLI}".format(chrom=v[0], start=v[1],
        end=v[2], gene=v[3], pLI=v[4]))
