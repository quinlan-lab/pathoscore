from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys
import pandas as pd

plt.rcParams["figure.figsize"]=(14,5)
with open(sys.argv[1], 'r') as clinicaltable:
    df = pd.read_csv(clinicaltable, sep='\t', index_col=0)
    df = df.T
    print(df)
    # header=clinicaltable.readline().strip().split("\t")
    # print (header)
    # for line in clinicaltable:
        # print(line.strip().split("\t"))
        

ax = df.plot.bar(rot=0,legend=False,width=.8,align='center', fontsize=18)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=18)
ax.set_ylabel("Clinical Utility", fontsize=18)
plt.tight_layout()
plt.savefig("".join(sys.argv[1].split(".cu.tsv",)[0])+"cubarplot.pdf", bbox_inches="tight")

#fig, ax = plt.subplots()
