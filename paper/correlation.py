from __future__ import print_function
metrics = ["CCR", "CADD", "DANN", "GERP", "MCAP", "MPC", "MTR", "REVEL", "RVIS", "VVP", "aloft_het", "aloft_lof", "aloft_rec", "fathmm_non", "fathmm_coding", "fitCons", "MetaSVM", "missense_z", "pLI", "phastCons", "polyphen2_hvar", "SIFT", "mis_badness", "Grantham", "BLOSUM"]
from cyvcf2 import VCF
import pandas as pd
from argparse import ArgumentParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.backends.backend_pdf
from matplotlib import pyplot as plt
from matplotlib import ticker
#import seaborn as sns
#sns.set_style('white')
import sys
if sys.version_info.major > 2:
    xrange = range

def key_val(arg):
    return arg.split(",")
parser = ArgumentParser(description = "create correlation matrix for metrics across different truth sets")
parser.add_argument("-f", "--file", help = "input vcfs", type=key_val, nargs = "+")
parser.add_argument("-n", "--name", help = "evaluation set names/titles", nargs = "+")
parser.add_argument("-o", "--output", help = "output file")
args = parser.parse_args()

corr = {}
for var, name in zip(args.file, args.name):
    patho, benign = var
    vcfpatho = VCF(patho)
    vcfbenign = VCF(benign)
    p = pd.DataFrame([variant.INFO.get(key) for key in metrics] for variant in vcfpatho)
    b = pd.DataFrame([variant.INFO.get(key) for key in metrics] for variant in vcfbenign)
    p.columns = metrics
    b.columns = metrics
    d = p.append(b)
    corr[name] = d.corr(method = 'kendall')
    print(corr)
print(d)

plt.rcParams["figure.figsize"]=(5,5)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
print(len(args.file))
#fig, ax = plt.subplots(len(args.file), 1)
for i, name in enumerate(args.name):
    fig, ax = plt.subplots(1)
    g = ax.matshow(corr[name], interpolation = 'nearest', cmap = plt.cm.Blues)
    fig.colorbar(g, ax=ax)
    ax.xaxis.set_major_locator(ticker.LinearLocator(numticks = len(metrics)))
    ax.yaxis.set_major_locator(ticker.LinearLocator(numticks = len(metrics)))
    ax.set_xticklabels(metrics)
    ax.set_yticklabels(metrics)
    plt.setp(ax.get_xticklabels(), fontsize = 10, rotation = 'vertical')
    plt.title(name, y=1.35)
pdf = matplotlib.backends.backend_pdf.PdfPages(args.output)
for fig in xrange(1, plt.gcf().number+1): ## will open an empty extra figure :(
    pdf.savefig(fig, bbox_inches = 'tight')
pdf.close()
