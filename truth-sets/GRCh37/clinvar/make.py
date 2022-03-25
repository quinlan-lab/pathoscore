import sys
import re
from cyvcf2 import VCF, Writer

date = sys.argv[2]

vcf = VCF(sys.argv[1])

fhs = {'benign': Writer('clinvar-benign.%s.vcf' % date, vcf),
       'benign-likely_benign': Writer('clinvar-benign-likely_benign.%s.vcf' % date, vcf),
       'pathogenic': Writer('clinvar-pathogenic.%s.vcf' % date, vcf),
       'expert': Writer('clinvar-pathogenic-expert.%s.vcf' % date, vcf),
       'likely_pathogenic-pathogenic': Writer('clinvar-pathogenic-likely_pathogenic.%s.vcf' % date, vcf),
       }
wtr = fhs['pathogenic']

fhs['likely_benign'] = fhs['benign-likely_benign']
fhs['likely_pathogenic'] = fhs['likely_pathogenic-pathogenic']

# benigns also get written just to the joint file
lookup = {}
lookup['benign'] = 'benign-likely_benign'
lookup['pathogenic'] = 'likely_pathogenic-pathogenic'

exclude = "criteria_provided,_conflicting_interpretations no_assertion_criteria_provided no_assertion_provided no_interpretation_for_the_single_variant".split()

for v in vcf:
    clnsig = v.INFO.get('CLNSIG')
    if clnsig is None: continue
    if len(v.ALT) == 0 or v.REF == v.ALT[0]: continue

    # exclude things with questionable review status
    crs = v.INFO.get('CLNREVSTAT', '').lower()
    if crs in exclude or any(c in exclude for c in crs.split(",")): continue
    # we filter this out as well, because it indicates a potential false positive variant
    if v.INFO.get('SSR', ''): continue
    info = v.INFO
    for k, _ in info:
        if k == 'CLNSIG': continue
        del info[k]

    # make a key that matches the fhs keys
    key = "-".join(sorted(set(re.split("/|,", clnsig.lower()))))
    if not key in fhs:
        continue
    fhs[key].write_record(v)
    if key in lookup:
        okey = lookup[key]
        fhs[okey].write_record(v)

    if 'expert' in crs or 'practice' in crs:
        fhs['expert'].write_record(v)

for fh in fhs.values():
    fh.close()
