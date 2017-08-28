from __future__ import print_function
import os
import sys
import math
import toolshed as ts
from cyvcf2 import VCF
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from scipy.stats import binom_test
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.preprocessing import minmax_scale
sns.set_style('white')

__version__ = "0.1.2"

WIDTH = 7

cmd = "vcfanno -lua {lua} -p {p} {conf} {query_vcf} | bgzip -c > {out_vcf}"

def infos(path):
    infos = []
    for x in ts.nopen(path):
        if x[1] != "#": break
        if not "INFO" in x: continue
        infos.append(x.split("ID=")[1].split(",")[0])
    return infos

def isfunctional(csqs):
    for csq in csqs.split(","):
        eff = csq.split("|", 2)[0]
        if eff in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon', 'rare_amino_acid',
                     'missense', 'protein_altering', 'frameshift'):
            return True
    return False

def evaluate(vcfs, fields, inverse_fields, include=None, functional=False):
    scored = {}
    unscored = {}
    for f in fields + inverse_fields:
        scored[f] = [[], []]
        unscored[f] = [0, 0]

    fields = [(f, False) for f in fields] + [(f, True) for f in inverse_fields]
    common_pathogenic = 0
    scorable = [[0, 0], [0, 0]]

    include_skipped = 0
    functional_skipped = 0

    for i, vcf in enumerate(vcfs):
        for v in VCF(vcf):
            if v.REF == v.ALT[0]:
                continue
            is_pathogenic = i == 0
            if include and v.INFO.get(include) is not None:
                include_skipped += 1
                continue

            is_indel = 1 - int(len(v.REF) == 1 and len(v.ALT[0]) == 1)

            if functional:
                csq = v.INFO.get("BCSQ")
                if csq is None or not isfunctional(csq):
                    functional_skipped += 1
                    continue

            if is_pathogenic and v.INFO.get('_exclude'):
                common_pathogenic += 1
                continue
            scorable[is_pathogenic][is_indel] += 1

            for f, invert in fields:
                score = v.INFO.get(f)
                if score is None or score == "NA":
                    unscored[f][is_pathogenic] += 1
                    continue
                try:
                    score = float(score)
                except: # handle multiple alts by recording both
                    iscores = map(float, score.split(","))
                    if invert:
                        score = min(iscores)
                    else:
                        score = max(iscores)

                if math.isnan(score):
                    unscored[f][is_pathogenic] += 1
                    continue
                if invert:
                    score = -score

                scored[f][is_pathogenic].append(score)

    methods = [f for f, _ in fields]
    for f in methods:
        for i in (0, 1):
            arr = np.array(scored[f][i], dtype=float)
            if np.any(np.isinf(arr)):
                imax = np.max(arr[~np.isinf(arr)])
                arr[np.isinf(arr)] = imax
                scored[f][i] = list(arr)

    print("unscored:", unscored)
    print("scored:", {k: {'benign': len(v[0]), 'pathogenic': len(v[1])} for k, v in scored.items()})
    print("pathogenics excluded (via '_exclude' flag): %d" % common_pathogenic)
    if include:
        print("variants skipped for lack of include: %d" % include_skipped)
    if functional:
        print("variants skipped as not functional: %d" % functional_skipped)

    print("scorable sites: benign (snp/indel): (%d/%d), pathogenic: (%d/%d)" % tuple(scorable[0] + scorable[1]))
    return methods, scored, unscored, scorable

def plot(score_methods, scored, unscored, scorable, prefix, title=None, suffix="png"):

    bar_colors = sns.color_palette()[:2]
    bar_colors = [bar_colors[0], tuple(x * 0.85 for x in bar_colors[0]), (0.9, 0.9, 0.9), (0.8, 0.8, 0.8)]

    if len(score_methods) <= 10:
        sns.set_palette(sns.color_palette("Vega10", 10))
    else:
        sns.set_palette(sns.color_palette("Vega20", 20))
    plt.figure(figsize=(WIDTH, 6))

    jcurves = {}
    for f in score_methods:
        if len(scored[f][0]) == 0:
            print("skipping %s because no negatives" % f, file=sys.stderr)
            continue
        if len(scored[f][1]) == 0:
            print("skipping %s because no positives" % f, file=sys.stderr)
            continue

        scores = scored[f][0] + scored[f][1]
        truth = ([0] * len(scored[f][0])) + ([1] * len(scored[f][1]))
        fpr, tpr, thresh = roc_curve(truth, scores, pos_label=1, drop_intermediate=True)
        auc_score = auc(fpr, tpr)

        ji = np.argmax(tpr - fpr)
        J = tpr[ji] - fpr[ji]
        score_at_maxJ = thresh[ji]
        jcurves[f] = tpr + (1 - fpr) - 1, score_at_maxJ, thresh
        plt.plot(fpr, tpr, label="%s (%.2f, %.2f)" % (f, auc_score, J))
        plt.plot([fpr[ji]], [tpr[ji]], 'ko')

    plt.plot([0, 1], [0, 1], linestyle='--', color='#777777', zorder=-1)

    # order is scored path, benign then unscored path, benign
    score_counts = [
        np.array([len(scored[key][1]) for key in score_methods]),
        np.array([len(scored[key][0]) for key in score_methods]),
        np.array([unscored[key][1] for key in score_methods]),
        np.array([unscored[key][0] for key in score_methods]),
        ]
    labels = ['scored pathogenic', 'scored benign', 'unscored pathogenic', 'unscored benign']
    pct_variants_scored = 100.0 *(score_counts[0] + score_counts[1]).astype(float) / np.array(score_counts).sum(axis=0)
    sns.despine()

    plt.xlim(-0.004, 1)
    plt.ylim(0, 1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    legend = plt.legend(loc="lower right", title="%s (AUC, J index)" % "method", handletextpad=1)
    if title:
        plt.title(title)
    plt.savefig(prefix + ".roc." + suffix)
    plt.close()

    for f in score_methods:
      jc, cutoff, thresh = jcurves[f]
      idx = np.argmax(jc)
      xs = minmax_scale(thresh)
      J = jc[idx]
      plt.plot(xs, jc, label="%s (%.2f @ %.2f)" % (f, J, cutoff))
      plt.plot([xs[idx]], [jc[idx]], 'ko')

    sns.despine()
    plt.ylabel('J-score')
    plt.xlabel('Normalized score')
    leg = plt.legend(title="method (J-index @ score)", bbox_to_anchor=(1, 1))
    plt.savefig(prefix + ".J." + suffix, bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.close()

    # histogram of J scores.
    Js = []
    for f in score_methods:
      jc, cutoff, thresh = jcurves[f]
      idx = np.argmax(jc)
      Js.append(jc[idx])
    inds = 0.1 + np.array(list(range(len(score_methods))))
    width = 0.72
    bars = plt.bar(inds, Js)
    colors = sns.color_palette()
    ymax, _ = plt.ylim()
    ax = plt.gca()

    for i, bar in enumerate(bars):
        bar.set_color(colors[i])
        label = "%.2f" % Js[i]
        height = 0.005 * ymax + bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height, label, ha='center', va='bottom', zorder=10)

    plt.xticks(np.array(inds), score_methods, rotation=30, ha='right')
    sns.despine()
    plt.ylabel('J-index')
    plt.savefig(prefix + ".Jbar." + suffix,  bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.close()

    # histogram of scored/unscored by pathogenic/benign

    inds = 0.1 + np.array(list(range(len(score_methods))))
    width = 0.72

    plt.figure(figsize=(WIDTH, 4))
    bottom = np.zeros_like(score_counts[0])
    shapes = []
    for i, sc in enumerate(score_counts):
        shapes.append(plt.bar(inds, sc, width, bottom=bottom, color=bar_colors[i], label=labels[i])[0])
        bottom += sc
    plt.xticks(np.array(inds), score_methods, rotation=30, ha='right')
    sns.despine()
    plt.ylabel('Variants')
    #ph = [plt.plot([],marker="", ls="")[0]]*2
    leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.20), ncol=2)
    #plt.tight_layout(h_pad=0.2)
    plt.savefig(prefix + ".stats." + suffix,  bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.close()

    # get the red and blue colors for path, benign
    sns.set_palette(sns.color_palette("Set1", 10))
    fig, axes = plt.subplots(len(score_methods), figsize=(WIDTH,
        2*len(score_methods)))
    try:
        axes[0]
    except:
        axes = (axes,)
    for i, f in enumerate(score_methods):
        ax = axes[i]

        vals = np.array(scored[f][1])
        step_plot(vals, ax, label="pathogenic", alpha=0.85)

        vals = np.array(scored[f][0])
        step_plot(vals, ax, label="benign", alpha=0.85)

        ax.set_xlabel(f)
        rng = vals.max() - vals.min()
        ax.set_xlim(vals.min() - 0.01 * rng, vals.max() + 0.01 * rng)
        ax.set_ylabel("Frequency")

    axes[0].legend(loc='upper left')
    if title:
        plt.suptitle(title)
    plt.tight_layout()
    sns.despine()
    plt.savefig(prefix + ".step." + suffix)
    write_html(prefix, scorable, title, suffix)

def write_html(prefix, scorable, title=None, suffix="png"):
    import datetime
    fh = open(prefix + ".overview.html", "w")
    fh.write("""<html>
<title>pathoscore summary {title}</title>
<body>

<pre>
date: {date}
created with <b><a href="https://github.com/quinlan-lab/pathoscore">pathoscore</a></b> version: {version}
{title}
</pre>

<i>pathoscore evaluates variant pathogenicity tools and scores.</i>

<p>
In this evaluation, there were <b>{pathogenic} pathogenic</b> ({pathogenic_pct_indel:.1f}% indels)
and <b>{benign} benign</b> ({benign_pct_indel:.1f}% indels) variants that could be scored.
</p>


<h3>Distribution of variants scored</h3>
<img src="{prefix}.stats.{suffix}"/>

<h3>Receiver Operating Characteristic Curve</h3>
<img src="{prefix}.roc.{suffix}"/>

<h3>Youden's J Statistic</h3>
<img src="{prefix}.Jbar.{suffix}"/>

<h3>Youden's J Statistic (distribution)</h3>
<img src="{prefix}.J.{suffix}"/>

<h3>Step Plot of Scores</h3>
<img src="{prefix}.step.{suffix}"/>

<pre>
invocation: {invocation}
</pre>
</body>
</html>""".format(prefix=prefix.split(os.path.sep)[-1], date=datetime.date.today(),
                  title=("for " + title) if title else "",
                  invocation=" ".join(sys.argv),
                  pathogenic=sum(scorable[1]),
                  pathogenic_pct_indel=100.0*scorable[1][1] / float(sum(scorable[1])),
                  benign=sum(scorable[0]),
                  benign_pct_indel=100.0*scorable[0][1] / float(sum(scorable[0])),
                  suffix=suffix,
                  version=__version__))
    fh.close()
    print("wrote overview to %s" % fh.name, file=sys.stderr)


def annotate(args):
    scores = [x.split(":") for x in args.scores]
    assert all(len(x) == 4 for x in scores), "scores must be specified as quartets of path:dest:source:op"

    fh = open("%s.conf" % args.prefix, "w")
    lua_fields = ['"%s"' % i for i in infos(args.query_vcf)]
    names = []
    for path, name, field, op in scores:
        names.append(name)
        lua_fields.append('"%s"' % name)
        if not field.isdigit():
            field = '"%s"' % field
            col = "fields"
        else:
            col = "columns"
        fh.write("""[[annotation]]
file="{path}"
names=["{name}"]
{col}=[{field}]
ops=["{op}"]
\n""".format(**locals()))

    for exclude in (args.exclude or []):
        field = """fields=["AF"]""" if exclude.endswith(".vcf.gz") else """columns=[1]"""
        fh.write("""
[[annotation]]
file="{path}"
names=["_exclude"]
{field}
ops=["flag"]
\n""".format(path=exclude, field=field))


    if args.conf:
        fh.write("\n")
        fh.write(open(args.conf).read())
    fh.close()

    if not args.lua:
        args.lua = """<(echo "")"""

    out = args.prefix + ".vcf.gz"

    fcmd = cmd.format(p=args.procs, conf=fh.name, query_vcf=args.query_vcf, out_vcf=out, lua=args.lua)
    print(fcmd)
    for d in ts.nopen("|" + fcmd):
        print(d)

def step_plot(vals, ax, **kwargs):
    p, p_edges = np.histogram(vals, bins=kwargs.pop('bins', 50), range=[vals.min(), vals.max()], normed=True)
    if p[0] > 2 * p[1:].max():
        ax.set_yscale('log', basey=2)
        ax.get_yaxis().get_major_formatter().labelOnlyBase = False
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    p = list(p)
    p.append(p[-1])
    ax.plot(p_edges, p, ls='steps', lw=1.9, **kwargs)


def add_eval_args(p):
    p.add_argument("query_vcf", nargs="+", help="vcf(s) to annotate if 2 are specified it must be pathogenic and then benign")
    p.add_argument("--score-columns", "-s", action="append", help="info fields on which to base evaluation.",
                     default=[])
    p.add_argument("--inverse-score-columns", "-i", action="append", default=[],
            help="like score columns but lower score is more constrained")
    p.add_argument("--include", help="only evaluate variants that have this Flag in the INFO field. (Useful for specifying include regions)")
    p.add_argument("--functional", action="store_true", default=False,
    help="only evaluate variants that are missense or loss-of-function per annotation from bcftools csq. Default is to evaluate all variants in the input")
    p.add_argument("--prefix", default="pathoscore", help="prefix for output files")
    p.add_argument("--title", help="optional title for figure")
    p.add_argument("--suffix", help="plot type", choices=("png", "svg", "pdf",
    "eps"), default="svg")


if __name__ == "__main__":
    from argparse import ArgumentParser

    p = ArgumentParser()
    subps = p.add_subparsers(help="sub-command", dest="command")

    ### annotation ###
    pan = subps.add_parser("annotate")

    pan.add_argument("--procs", "-p", default=3, help="number of processors to use for vcfanno")
    pan.add_argument("--exclude", default=[], action="append", help="optional exclude vcf or bed to filter supposed pathogenic variants")
    pan.add_argument("--prefix", default="pathoscore", help="prefix for output files")
    pan.add_argument("--conf", help="optional vcfanno conf file that will also be used for annotation")
    pan.add_argument("--lua", help="optional path to lua file if it's needed by the --conf argument")
    pan.add_argument("--scores", "-s", default=[], action="append", help="format of path:name:field:op e.g. some.bed:myscore:4:self or cadd.vcf:cadd:PHRED:concat that give the path of the annotation file, the name in the output, and the column in the input respectively. may be specified multiple times. op is one of those specified here: https://github.com/brentp/vcfanno#operations")
    pan.add_argument("query_vcf", help="vcf to annotate")

    ### evaluation ###
    pev = subps.add_parser("evaluate")
    add_eval_args(pev)

    a = p.parse_args()

    if a.command == "annotate":
        annotate(a)
    elif a.command == "evaluate":
        methods, scored, unscored, scorable = evaluate(a.query_vcf,
                a.score_columns, a.inverse_score_columns, include=a.include,
                functional=a.functional)
        plot(methods, scored, unscored, scorable, a.prefix, a.title, a.suffix)


