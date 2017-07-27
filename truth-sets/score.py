import sys
import os
import toolshed as ts
from cyvcf2 import VCF, Writer

HERE = os.path.dirname(__file__)

fname = os.path.join(HERE, "score-lookup.txt")

lookup = {(d['from_short'], d['to_short']): d for d in ts.reader(fname)}
maxg, maxmb, minb = 0, 0, 0
for v in lookup.values():
    if v['Grantham']:
        v['Grantham'] = int(v['Grantham'])
    else:
        del v['Grantham']
    if v['mis_badness']:
        v['mis_badness'] = float(v['mis_badness'])
    else:
        del v['mis_badness']
    v['BLOSUM'] = int(v['BLOSUM'])
    maxg = max(v.get('Grantham', maxg), maxg)
    maxmb = max(v.get('mis_badness', maxmb), maxmb)
    minb = min(v['BLOSUM'], minb)

froms = [x[0] for x in lookup.keys()]
for f in froms:
    lookup[(f, "*")] = dict(Grantham=maxg, mis_badness=maxmb, BLOSUM=minb)


def change(csq):
    """
    >>> change("missense|TMEM240|ENST00000378733|protein_coding|-|170P>170L|1470752G>A")
    ('170P', '170L', True)
    >>> change('synonymous|AGRN|ENST00000379370|protein_coding|+|268A|976629C>T')
    ('268A', '268A', False)
    """
    parts = csq.split("|")
    if len(parts) < 7:
        return 0, 0, False
    if not ">" in parts[-2]:
        return (parts[-2], parts[-2], False)

    aa_from, aa_to = parts[-2].split(">")
    # 533PLGQGYPYQYPGPPPCFPPAYQDPGFSYGSGSTGSQQSEGSKSSGSTRSSRRAPGREKERRAAGAGGSGSESDHTAPSGVGSSWRERPAGQLSRGSSPRSQASATAPGLPPPHPTTKAYTVVGGPPGGPPVRELAAVPPELTGSRQSFQKAMGNPCEFFVDIM*>533LWVRATPTSTRDPHPASRLPTRTRALAMAAAAPGVSRVKGAKAVGPPGAAAGPRAVRRSVGRRELGAVAVNRITRHRVGWGAAGESVRPASSAVAAAHAVRPRLPPRGSPRPTPRPRPIQWWGGHPGDPLSGSWLPSPRN*
    return aa_from, aa_to, True

def aasplit(aa):
    """
    >>> aasplit("533PLGQGYPYQYPGPPPCFPPAYQDPGFSYGSGSTGSQQSEGSKSSGSTRSSRRAPGREKERRAAGAGGSGSESDHTAPSGVGSSWRERPAGQLSRGSSPRSQASATAPGLPPPHPTTKAYTVVGGPPGGPPVRELAAVPPELTGSRQSFQKAMGNPCEFFVDIM")
    (533, 'PLGQGYPYQYPGPPPCFPPAYQDPGFSYGSGSTGSQQSEGSKSSGSTRSSRRAPGREKERRAAGAGGSGSESDHTAPSGVGSSWRERPAGQLSRGSSPRSQASATAPGLPPPHPTTKAYTVVGGPPGGPPVRELAAVPPELTGSRQSFQKAMGNPCEFFVDIM')
    """
    for i, v in enumerate(aa):
        if not v.isdigit():
            break
    return int(aa[:i]), aa[i:]

def score(csq, smax=None):
    """
    >>> exp = {'Grantham': 98, 'mis_badness': 0.3769446754, 'BLOSUM': -3}
    >>> score("missense|TMEM240|ENST00000378733|protein_coding|-|170P>170L|1470752G>A") == exp
    True

    >>> score("frameshift|DVL1|ENST00000378888|protein_coding|-|527HPAAPWPLGQGYPYQYPGPPPCFPPAYQDPGFSYGSGSTGSQQSEGSKSSGSTRSSRRAPGREKERRAAGAGGSGSESDHTAPSGVGSSWRERPAGQLSRGSSPRSQASATAPGLPPPHPTTKAYTVVGGPPGGPPVRELAAVPPELTGSRQSFQKAMGNPCEFFVDIM*>527PGLWVRATPTSTRDPHPASRLPTRTRALAMAAAAPGVSRVKGAKAVGPPGAAAGPRAVRRSVGRRELGAVAVNRITRHRVGWGAAGESVRPASSAVAAAHAVRPRLPPRGSPRPTPRPRPIQWWGGHPGDPLSGSWLPSPRN*|1273478GGGGGCAGCCGGGT>G") == {'BLOSUM': -3, 'Grantham': 215, 'mis_badness': 1.0}
    True
    """

    func_lookup = {'BLOSUM': min}
    f, t, ok = change(csq)
    if not ok:
        return smax

    fpos, faa = aasplit(f)
    tpos, taa = aasplit(t)

    if smax is None:
        smax = {'BLOSUM': 0, 'Grantham': 0, 'mis_badness': 0}
    keys = smax.keys()
    for f in set(faa):
        for t in set(taa):
            ftscore = lookup.get((f, t))
            if ftscore is None: continue
            for k in keys:
                try:
                    v = func_lookup.get(k, max)(ftscore[k], smax[k])
                except KeyError:
                    continue
                if v is None: continue
                if v != smax[k]:
                    smax[k] = ftscore[k]
    return smax

if __name__ == "__main__":
    import doctest
    doctest.testmod()

    vcf = VCF(sys.argv[1])
    vcf.add_info_to_header({'ID': 'Grantham', 'Description': 'grantham score',
                            'Type':'Integer', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'mis_badness', 'Description': 'missense badness score',
                            'Type':'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'BLOSUM', 'Description': 'blosum score',
                            'Type':'Integer', 'Number': '1'})


    wtr = Writer("-", vcf)

    for v in vcf:
        csqs = v.INFO.get("BCSQ", "").split(",")
        smax = None
        for csq in csqs:
            smax = score(csq, smax)
        if not smax:
            continue

        for k in ('Grantham', 'mis_badness', 'BLOSUM'):
            v.INFO[k.lower()] = smax[k]

        wtr.write_record(v)
    wtr.close()
