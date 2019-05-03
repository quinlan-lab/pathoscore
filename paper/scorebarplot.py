from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys
import pandas as pd

plt.rcParams["figure.figsize"]=(14,5)
with open(sys.argv[1], 'r') as clinvartable, open(sys.argv[2], 'r') as samochatable, open(sys.argv[3], 'r') as homsytable, open(sys.argv[4], 'r') as wiglertable:
    c = pd.read_csv(clinvartable, sep=",", index_col=0)
    s = pd.read_csv(samochatable, sep=",", index_col=0)
    h = pd.read_csv(homsytable, sep=",", index_col=0)
    w = pd.read_csv(wiglertable, sep=",", index_col=0)
    auc = pd.DataFrame([c["AUC"],s["AUC"],h["AUC"],w["AUC"]]).T
    j = pd.DataFrame([c["J"],s["J"],h["J"],w["J"]]).T
    f1 = pd.DataFrame([c["F1@J"],s["F1@J"],h["F1@J"],w["F1@J"]]).T
    columns = ["ClinVar","Samocha","Homsy","Wigler"]
    auc.columns, j.columns, f1.columns = columns, columns, columns

ax = auc.plot.bar(rot=0,legend=False,width=.8,align='center', fontsize = 16)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 18)
ax.set_ylabel("AUC", fontsize = 18)
ax.set_ylim(0.4, 1)
ax.set_xlabel("method", fontsize = 16)
plt.tight_layout()
plt.savefig("".join(sys.argv[1].rsplit("/",1)[0])+"/aucbarplot.pdf", bbox_inches="tight")

ax = j.plot.bar(rot=0,legend=False,width=.8,align='center', fontsize = 16)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 18)
ax.set_ylabel("Peak J", fontsize = 18)
ax.set_xlabel("method", fontsize = 16)
plt.tight_layout()
plt.savefig("".join(sys.argv[1].rsplit("/",1)[0])+"/jbarplot.pdf", bbox_inches="tight")

ax = f1.plot.bar(rot=0,legend=False,width=.8,align='center', fontsize = 16)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 18)
ax.set_ylabel("F1 @ Peak J", fontsize = 18)
ax.set_xlabel("method", fontsize = 16)
plt.tight_layout()
plt.savefig("".join(sys.argv[1].rsplit("/",1)[0])+"/f1barplot.pdf", bbox_inches="tight")
