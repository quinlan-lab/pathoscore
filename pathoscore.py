from __future__ import print_function
import os
import sys
import datetime
import string
import math
import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.3f')
import itertools

import toolshed as ts
from cyvcf2 import VCF
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score, f1_score
from scipy.stats import binom_test
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from collections import OrderedDict, defaultdict
import seaborn as sns
from sklearn.preprocessing import minmax_scale
sns.set_style('white')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16)

__version__ = "0.1.3"

WIDTH = 7

cmd = "vcfanno -lua {lua} -p {p} {conf} {query_vcf} | bgzip -c > {out_vcf}"

def get_genes(csqs):
    genes = set()
    for csq in csqs.split(","):
        genes.add(csq.split("|", 2)[1])
    return genes

def clinical_utility(scoredbygene, unscoredbygene, jindices, prefix, goi):
    genes = set(itertools.chain.from_iterable(map(str,scoredbygene[key].keys()) for key in scoredbygene.keys()))
    genes = genes.intersection(goi)
    cu = defaultdict(lambda: dict.fromkeys(genes, 0))
    tp, tn, fp, fn = 0, 0, 0, 0
    for method in scoredbygene:
        for gene in scoredbygene[method]:
            if gene not in genes:
                continue
            scored = len(scoredbygene[method][gene][0])+len(scoredbygene[method][gene][1])
            unscored = unscoredbygene[method][gene][0]+unscoredbygene[method][gene][1]
            j = jindices[method]
            fracvars = scored/float(scored+unscored)
            for score in scoredbygene[method][gene][1]:
                if score >= j:
                    tp += 1
                else:
                    fn += 1
            for score in scoredbygene[method][gene][0]:
                if score <= j:
                    tn += 1
                else:
                    fp += 1
            acc = (tp + tn)/float(tp + tn + fp + fn)
            cu[method][gene] = acc * fracvars
            tp, tn, fp, fn = 0, 0, 0, 0
    
    culist = []
    header = [{'title': 'Genes'}]+[{'title': i} for i in map(str,cu.keys())]
    for gene in genes:
        culist.append([gene] + ["{num:.3f}".format(num=cu[method][gene]) for method in cu])
    return culist, header

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
        for c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon', 'rare_amino_acid',
                     'missense', 'protein_altering', 'frameshift', 'inframe_insertion', 'inframe_deletion'):
            if c in eff or (('splice_donor' in eff or 'splice_acceptor' in eff) and 'coding_sequence' in eff): 
                return True
    return False

def evaluate(vcfs, fields, inverse_fields, include=None, functional=False, goi=set(['BRCA2','BRCA1','SCN1A','LDLR','MLH1','MSH2','DMD','ATM','FBN1','CFTR'])):
    scored = {}
    unscored = {}
    scoredbygene = {}
    unscoredbygene = {}
    for f in fields + inverse_fields:
        scored[f] = [[], []]
        unscored[f] = [0, 0]
        scoredbygene[f] = defaultdict(lambda:[[], []]) # left is scored benign, right is scored pathogenic
        unscoredbygene[f] = defaultdict(lambda:[0, 0]) # left is unscored benign, right is unscored pathogenic

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

            csq = v.INFO.get("BCSQ")
            if functional:
                if csq is None or not isfunctional(csq):
                    functional_skipped += 1
                    continue
            genes = get_genes(csq)

            if v.INFO.get('_exclude'): # excludes benigns as well if desired
                if is_pathogenic:
                    common_pathogenic += 1
                continue
            scorable[is_pathogenic][is_indel] += 1

            for f, invert in fields:
                score = v.INFO.get(f)
                if score is None or score == "NA":
                    unscored[f][is_pathogenic] += 1
                    for gene in genes:
                        unscoredbygene[f][gene][is_pathogenic] += 1
                    continue
                try:
                    score = float(score)
                except: # handle multiple alts by recording both
                    try:
                        iscores = map(float, score.split(","))
                        if invert:
                            score = min(iscores)
                        else:
                            score = max(iscores)
                    except ValueError:
                        score = float('nan')

                if math.isnan(score):
                    unscored[f][is_pathogenic] += 1
                    for gene in genes:
                        unscoredbygene[f][gene][is_pathogenic] += 1
                    continue
                if invert:
                    score = -score

                scored[f][is_pathogenic].append(score)
                for gene in genes:
                    scoredbygene[f][gene][is_pathogenic].append(score)

    methods = [f for f, _ in fields]
    for f in methods:
        for i in (0, 1):
            arr = np.array(scored[f][i], dtype=float)
            if np.any(np.isinf(arr)):
                imax = np.max(arr[~np.isinf(arr)])
                arr[np.isinf(arr)] = imax
                scored[f][i] = list(arr)
            for gene in genes:
                arr = scoredbygene[f][gene][i]
                if np.any(np.isinf(arr)):
                    imax = np.max(arr[~np.isinf(arr)])
                    arr[np.isinf(arr)] = imax
                    scoredbygene[f][gene][i] = list(arr)

    print("unscored:", unscored)
    print("scored:", {k: {'benign': len(v[0]), 'pathogenic': len(v[1])} for k, v in scored.items()})
    print("pathogenics excluded (via '_exclude' flag): %d" % common_pathogenic)
    if include:
        print("variants skipped for lack of include: %d" % include_skipped)
    if functional:
        print("variants skipped as not functional: %d" % functional_skipped)

    print("scorable sites: benign (snp/indel): (%d/%d), pathogenic: (%d/%d)" % tuple(scorable[0] + scorable[1]))
    return methods, scored, unscored, scorable, scoredbygene, unscoredbygene

def get_se(A, B, C, D):
    L = float(A * B) / (A + B)**3
    R = float(C * D) / (C + D)**3
    return np.sqrt(L + R)

def color_to_rgb(c):
    return "rgb(%d, %d, %d)" % (255*c[0], 255*c[1], 255*c[2])

def step_traces_to_json(st):
    score_layout = {
        "height": 250,
        "title": "HERE",
        "yaxis": {"title": "Frequency"},
    }

    for method, tr in st.items():
        div = "score_step_%s" % method
        js = json.dumps(tr)
        score_layout["title"] = method
        yield "Plotly.newPlot('%s', %s, %s)" % (div, js, json.dumps(score_layout))

def plot(score_methods, scored, unscored, scorable, prefix, title=None, suffix="png", goi=set(['BRCA2','BRCA1','SCN1A','LDLR','MLH1','MSH2','DMD','ATM','FBN1','CFTR'])):

    pr_base=(scorable[1][0]+scorable[1][1])/float(scorable[0][0]+scorable[0][1]+scorable[1][0]+scorable[1][1])

    bar_colors = sns.color_palette()[:2]
    bar_colors = [bar_colors[0], tuple(x * 0.85 for x in bar_colors[0]), (0.9, 0.9, 0.9), (0.8, 0.8, 0.8)]

    if len(score_methods) <= 10:
        try:
            sns.set_palette(sns.color_palette("Vega10", 10))
        except ValueError:
            sns.set_palette(sns.color_palette("tab10", 10))
    else:
        try:
            sns.set_palette(sns.color_palette("Vega20", len(score_methods)))
        except ValueError:
            sns.set_palette(sns.color_palette("tab20", len(score_methods)))
    colors = sns.color_palette()
    fig, ax = plt.subplots(figsize=(WIDTH, 6))
    fig2, ax2 = plt.subplots(figsize=(WIDTH, 6))

    jcurves = {}
    roc_traces = []
    roc_traces.append({
       'x': [0, 1],
       'y': [0, 1],
       'mode': 'lines',
       'showlegend': False,
       'hoverinfo': 'none',
       'line': {'color': 'rgb(200, 200, 200)', 'width': 3, 'dash': 'dash'}
    })
    pr_traces = []
    pr_traces.append({
    'x': [0, 1],
    'y': [pr_base, pr_base],
    'mode': 'lines',
    'showlegend': False,
    'hoverinfo': 'none',
    'line': {'color': 'rgb(200, 200, 200)', 'width': 3, 'dash': 'dash'}
    })
    jbar_trace = [{
        'x': score_methods,
        'y': [],
        'type': 'bar',
        'text': score_methods,
        'marker': {'color': [color_to_rgb(c) for c, m in zip(colors, score_methods)]},
        'error_y': {
            'type': 'data',
            'array': [],
            'visible': True,
        }
    }]
    jdist_traces = []
    jindices = {}
    output = OrderedDict((k, []) for k in ('method', 'J', 'score@J', 'se(J)', 'TPR@J', 'FPR@J', 'AUC', 'TP@J', 'FP@J', 'TN@J', 'FN@J', 'F1@J'))
    for i, f in enumerate(score_methods):
        if len(scored[f][0]) == 0:
            print("skipping %s because no negatives" % f, file=sys.stderr)
            continue
        if len(scored[f][1]) == 0:
            print("skipping %s because no positives" % f, file=sys.stderr)
            continue


        scores = scored[f][0] + scored[f][1]
        truth = ([0] * len(scored[f][0])) + ([1] * len(scored[f][1]))
        fpr, tpr, thresh = roc_curve(truth, scores, pos_label=1, drop_intermediate=True)
        precision, recall, thresholds = precision_recall_curve(truth, scores, pos_label=1)
        auc_score = auc(fpr, tpr)
        ap = average_precision_score(truth, scores, average = 'macro')

        ji = np.argmax(tpr - fpr)
        J = tpr[ji] - fpr[ji]
        S = score_at_maxJ = thresh[ji]
        jindices[f] = S
        jcurves[f] = tpr + (1 - fpr) - 1, score_at_maxJ, thresh
        jbar_trace[0]['y'].append(round(J, 3))

        y_pred = [1 if y > S else 0 for y in scores]
        f1 = f1_score(truth, y_pred, average = 'binary')

        # naming from Youden
        # A: true+
        # B: false-
        # C: false+
        # D: true-
        A = sum(s >= S and t == 1 for s, t in zip(scores, truth))
        B = sum(s < S and t == 1 for s, t in zip(scores, truth))
        C = sum(s >= S and t == 0 for s, t in zip(scores, truth))
        D = sum(s < S and t == 0 for s, t in zip(scores, truth))
        jse = get_se(A, B, C, D)
        jcurves[f] = tpr + (1 - fpr) - 1, score_at_maxJ, thresh, jse
        jbar_trace[0]['error_y']['array'].append(round(jse, 3))


        output['method'].append(f)
        output['J'].append(J)
        output['score@J'].append(S)
        output['se(J)'].append(jse)
        output['TPR@J'].append(tpr[ji])
        output['FPR@J'].append(fpr[ji])
        output['AUC'].append(auc_score)
        output['TP@J'].append(A)
        output['FP@J'].append(C)
        output['TN@J'].append(D)
        output['FN@J'].append(B)
        output['F1@J'].append(f1)
        label = "%s (AUC: %.2f, Peak J-score: %.2f)" % (f, auc_score, J)
        label2 = "%s (%.2f, %.2f)" % (f, auc_score, J)
        roc_traces.append({
           'x': list(np.round(fpr, 3)),
           'y': list(np.round(tpr, 3)),
           'text': ["<b>%s</b> FPR: %.2f, TPR: %.2f at score: %.2f" % (f, ff, t, s) for ff,t,s in zip(fpr, tpr, thresh)],
           'mode': 'lines',
           'hoverinfo': 'text',
           'line': {'color': color_to_rgb(colors[i])},
           'name': label
        })
        roc_traces.append({
          'x': [float(fpr[ji])],
          'y': [float(tpr[ji])],
          'marker': {'color': color_to_rgb(colors[i]), 'size': 10},
          'text': ["<b>J-index</b>: %.2f => FPR: %.2f, TPR: %.2f at score: %.2f" % (J, fpr[ji], tpr[ji], S)],
          'mode': 'markers',
           'hoverinfo': 'text',
          'showlegend': False,
          'type': 'scatter',
        })
        label3 = "%s (F1 score @ Peak J: %.2f)" % (f, f1)
        label4 = "%s (%.2f)" % (f, f1)
        pr_traces.append({
           'x': list(np.round(recall, 3)),
           'y': list(np.round(precision, 3)),
           'text': ["<b>%s</b> Recall: %.2f, Precision: %.2f at score: %.2f" % (f, ff, t, s) for ff,t,s in zip(recall, precision, thresholds)],
           'mode': 'lines',
           'hoverinfo': 'text',
           'line': {'color': color_to_rgb(colors[i])},
           'name': label3
        })
        ind = list(recall).index(tpr[ji])
        pr_traces.append({
          'x': [float(recall[ind])],
          'y': [float(precision[ind])],
          'marker': {'color': color_to_rgb(colors[i]), 'size': 10},
          'text': ["<b>F1 score</b>: %.2f => Recall: %.2f, Precision: %.2f at score: %.2f" % (f1, recall[ind], precision[ind], S)],
          'mode': 'markers',
           'hoverinfo': 'text',
          'showlegend': False,
          'type': 'scatter',
        })
        ax.plot(fpr, tpr, label=label2)
        ax.plot([fpr[ji]], [tpr[ji]], 'ko')
        ax2.plot(recall, precision, label=label4)
        ax2.plot(recall[ind], precision[ind], 'ko')

    ax.plot([0, 1], [0, 1], linestyle='--', color='#777777', zorder=-1)
    ax2.axhline(pr_base, linestyle='--', color='#777777', zorder=-1)
    #baseline for precision recall is dependent on P:N ratio

    df = pd.DataFrame(output)
    df.to_csv(prefix + ".csv", index=False, float_format="%.4f")
    sys.stderr.write("pathoscore: wrote csv file to %s.csv\n" % prefix)

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

    ax.set_xlim(-0.004, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False Positive Rate", fontsize=16)
    ax.set_ylabel("True Positive Rate", fontsize=16)
    ax2.set_xlim(-0.004, 1)
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("Recall", fontsize=16)
    ax2.set_ylabel("Precision", fontsize=16)
    legend = ax.legend(loc="lower right", title="%s (AUC, J index)" % "method", handletextpad=1, fontsize=12)
    legend = ax2.legend(loc="lower right", title="%s (F1 Score @ Peak J)" % "method", handletextpad=1, fontsize=12)
    if title:
        plt.title(title)
    fig.savefig(prefix + ".roc." + suffix)
    fig2.savefig(prefix + ".pr." + suffix)
    plt.close()

    fig, ax = plt.subplots(figsize=(WIDTH, 6))
    for i, f in enumerate(score_methods):
      jc, cutoff, thresh, se = jcurves[f]
      idx = np.argmax(jc)
      xs = minmax_scale(thresh)
      J = jc[idx]
      label = "%s (Peak J: %.2f @ score: %.2f)" % (f, J, cutoff)
      ax.plot(xs, jc, label=label)
      ax.plot([xs[idx]], [jc[idx]], 'ko')
      jdist_traces.append({
           'x': list(np.round(xs, 3)),
           'y': list(np.round(jc, 3)),
           'text': ["<b>%s</b> score: %.2f, J: %.2f" % (f,tr, jj ) for tr,jj  in zip(thresh, jc)],
           'mode': 'lines',
           'hoverinfo': 'text',
           'line': {'color': color_to_rgb(colors[i])},
           'name': label,
      })
      jdist_traces.append({
           'x': [round(xs[idx], 3)],
           'y': [round(jc[idx], 3)],
           'text': ["<b>%s</b> score: %.2f, J: %.2f" % (f,xs[idx], jc[idx])],
           'mode': 'markers',
           'hoverinfo': 'text',
           'showlegend': False,
           'type': 'scatter',
           'marker': {'color': color_to_rgb(colors[i]), 'size': 10},
           'name': label,
      })

    sns.despine()
    ax.set_ylabel('J-score', fontsize=16)
    ax.set_xlabel('Normalized score', fontsize=16)
    leg = ax.legend(title="method (J-index @ score)", bbox_to_anchor=(1, 1), fontsize=14)
    plt.savefig(prefix + ".J." + suffix, bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.close()

    # histogram of J scores.
    Js, errors = [], []
    for f in score_methods:
      jc, cutoff, thresh, se = jcurves[f]
      errors.append(se)
      idx = np.argmax(jc)
      Js.append(jc[idx])
    inds = 0.1 + np.array(list(range(len(score_methods))))
    width = 3.72
    bars = plt.bar(inds, Js, yerr=errors, error_kw=dict(capsize=6))
    ymax, _ = plt.ylim()
    ax = plt.gca()

    for i, bar in enumerate(bars):
        bar.set_color(colors[i])
        label = "%.2f" % Js[i]
        height = 0.005 * ymax + bar.get_height()
        #ax.text(bar.get_x() + bar.get_width()/2, height, label, ha='center', va='bottom', zorder=10)

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
    pb_colors = sns.color_palette()
    fig, axes = plt.subplots(len(score_methods), figsize=(WIDTH,
        2*len(score_methods)))
    try:
        axes[0]
    except:
        axes = (axes,)
    step_traces = OrderedDict()
    for i, f in enumerate(score_methods):
        step_traces[f] = [None, None]
        ax = axes[i]

        vals = np.array(scored[f][1])
        px, py = step_plot(vals, ax, label="pathogenic", alpha=0.85)
        step_traces[f][0] = {
           'x': list(np.round(px, 3)),
           'y': list(np.round(py, 3)),
           'mode': 'lines',
           'line': {'color': color_to_rgb(pb_colors[0]), 'shape': 'hv'},
           'name': "pathogenic",
        }

        vals = np.array(scored[f][0])
        bx, by = step_plot(vals, ax, label="benign", alpha=0.85)
        step_traces[f][1] = {
           'x': list(np.round(bx, 3)),
           'y': list(np.round(by, 3)),
           'mode': 'lines',
           'line': {'color': color_to_rgb(pb_colors[1]), 'shape': 'hv'},
           'name': "benign",
        }

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

    score_step_divs = "\n".join(['<div id="score_step_%s"></div>' % s for s in score_methods])


    return jindices, score_methods, score_counts, roc_traces, pr_traces, jbar_trace, jdist_traces, score_step_divs, step_traces

def serialize(arr):
    return "[%s]" % ",".join([("%.3f" % v).rstrip("0").rstrip(".") for v in arr])

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

<h3>Precision-Recall Curve</h3>
<img src="{prefix}.pr.{suffix}"/>

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
    for path, name, field, op in scores:
        fields = []
        name = name.split(",")
        op = op.split(",")
        lua_fields.append('"%s"' % name)
        field = field.split(",")
        for f in field:
            if not f.isdigit():
                f = '"%s"' % f
                col = "fields"
            else:
                col = "columns"
            fields.append(int(f))
        fh.write("""[[annotation]]
file="{path}"
names={name}
{col}={fields}
ops={op}
\n""".format(**locals()))
        print (name, op, fields, col)

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
    p, p_edges = np.histogram(vals, bins=kwargs.pop('bins', 50), range=[vals.min(), vals.max()])
    p = p.astype(float) / p.sum()
    p = list(p)
    p.append(p[-1])
    ax.plot(p_edges, p, ls='steps', lw=1.9, **kwargs)
    return p_edges, p

def plotly_html(score_methods, score_counts, roc_traces, pr_traces, jbar_trace, jdist_traces, score_step_divs, step_traces, scorable, prefix, cu=[], header=[]):
    tmpl = string.Template(open(os.path.join(os.path.dirname(__file__), "tmpl.html")).read())
    with open(prefix + ".html", "w") as html:
        html.write(tmpl.substitute(methods=score_methods,
                        scored_pathogenic=serialize(score_counts[0]),
                        scored_benign=serialize(score_counts[1]),
                        unscored_pathogenic=serialize(score_counts[2]),
                        unscored_benign=serialize(score_counts[3]),
                        roc_data=json.dumps(roc_traces),
                        pr_data=json.dumps(pr_traces),
                        Jbar_data=json.dumps(jbar_trace),
                        Jdist_data=json.dumps(jdist_traces),
                        score_step_divs=score_step_divs,
                        plotly_score_steps="\n".join(step_traces_to_json(step_traces)),
                        command=" ".join(sys.argv),
                        version=__version__,
                        date=(datetime.date.today()),
                        n_benign=sum(scorable[0]),
                        n_pathogenic=sum(scorable[1]),
                        path_indel_pct="%.1f" % (100.0*scorable[1][1] /
                            float(sum(scorable[1]))),
                        benign_indel_pct="%.1f" % (100.0*scorable[0][1] /
                            float(sum(scorable[0]))),
                        culist=cu,
                        header=header
                  ))

def add_eval_args(p):
    p.add_argument("query_vcf", nargs="+", help="vcf(s) to annotate if 2 are specified it must be pathogenic and then benign")
    p.add_argument("--score-columns", "-s", action="append", help="info fields on which to base evaluation.",
                     default=[])
    p.add_argument("--inverse-score-columns", "-i", action="append", default=[],
            help="like score columns but lower score is more constrained")
    p.add_argument("--include", help="only evaluate variants that have this Flag in the INFO field. (Useful for specifying include regions)")
    p.add_argument("--functional", action="store_true", default=False,
    help="only evaluate variants that are missense or loss-of-function per annotation from bcftools csq. Default is to evaluate all variants in the input")
    p.add_argument("--goi", help="file containing genes of interest in 1 column separated by newlines")
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
        if a.goi:
            goi = set()
            with open(a.goi, "r") as f:
                for line in f:
                    goi.add(line.strip())
        else:
            goi = set(['BRCA2','BRCA1','SCN1A','LDLR','MLH1','MSH2','DMD','ATM','FBN1','CFTR'])
        methods, scored, unscored, scorable, scoredbygene, unscoredbygene = evaluate(a.query_vcf,
                a.score_columns, a.inverse_score_columns, include=a.include,
                functional=a.functional, goi=goi)
        jindices, score_methods, score_counts, roc_traces, pr_traces, jbar_trace, jdist_traces, score_step_divs, step_traces = plot(methods, scored, unscored, scorable, a.prefix, a.title, a.suffix, goi)
        cu, header = clinical_utility(scoredbygene, unscoredbygene, jindices, a.prefix, goi)
        print ("\t".join([i["title"] for i in header]), file=open(a.prefix+".cu.tsv","w"))
        print ("\n".join(["\t".join(j) for j in cu]), file=open(a.prefix+".cu.tsv","a"))
        plotly_html(score_methods, score_counts, roc_traces, pr_traces, jbar_trace, jdist_traces, score_step_divs, step_traces, scorable, a.prefix, cu, header)


