pathoscore
==========

pathoscore evaluates variant pathogenicity tools and scores.

evaluating scores is hard because logic can be circular and benign and pathogenic sets are
hard to curate and evaluate.

`pathoscore` is software and datasets that facilitate applying evaluating pathogenicity scores.

The sections below describe the tools

Annotate
--------

Annotate a vcf with some scores (which can be bed or vcf)

```
python pathoscore.py annotate \
    --scores exac-ccrs.bed.gz:exac_ccr:14:max \
    --scores mpc.regions.clean.sorted.bed.gz:mpc_regions:5:max \
    --exclude "dataset ~= 'control'" \
    --population /data/gemini_install/data/gemini_data/ExAC.r0.3.sites.vep.tidy.vcf.gz \
    --conf combined-score.conf \
    testing-denovos.vcf.gz
```

The individual flags are described here:

### scores


The `scores` format is `path:name:column:op` where:

+ name becomes the new name in the INFO field.
+ column indicates the column number (or INFO name) to pull from the scores VCF.
+ op is a `vcfanno` operation.

### pathogenic

pathogenic is a `lua` expression that will be run on the dataset to indicate which variants are pathogenic.
It could also be something like `pathogenic ~= nil`.

### exclude

is a population VCF that is use to filter would-be pathogenic variants (as we know that common variants
can't be pathogenic).

### conf

an optional [vcfanno](https://github.com/brentp/vcfanno) conf file so users can specify exactly
how to annotate if they feel comfortable doing so.

This can also be used to specify vcfanno `[[postannotation]]` blocks, for example, to combine scores.

An example `conf` to combine 2 scores looks like:

```
[[postannotation]]
name="combined"
op="lua:exac_ccr+33\*cadd"
fields=["exac_ccr", "cadd"]
type="Float"
```

Evaluate
--------

```
python pathoscore.py evaluate \
    -s MPC \
    -s exac_ccr \
    -i mpc_regions \
    -s combined \
    pathogenic.vcf.gz \
    benign.vcf.gz
```

This will take the output(s) from `annotate` and create ROC curves and score distribution plots.
It assumes that the first VCF contains pathogenic variants and the 2nd contains benign variants.
It uses the columns specified via `-s` and `-i` as the scores.

`-i` indicates that lower scores are more constrained where as 

`-s` is for fields where higher scores are more constrained.


Install
-------

Download a [vcfanno binary](https://github.com/brentp/vcfanno/releases) for your system and make it available as
`vcfanno` on your `$PATH`

Then run:
```
pip install -r requirements.txt
```

Then you should be able to run the evaluation scripts.

Truth Sets
----------

Part of `pathoscore` is to provide curated truth sets that can be used for evaluation.

These are kept in `truth-sets/`. Each set has a benign and/or a pathogenic set. 

Pull-requests for recipes that add new truth sets are welcomed. These should include a `make.sh`
script that, when run will pull from the original data source and make a benign and/or pathogenic
vcf that is bgzipped and tabixed and made as small as possible (see the clinvar example for how
to remove unneeded fields from the INFO field).

Currently we have:

### clinvar

+ clinvar pathogenics are either `Pathogenic` or `Likely-Pathogenic` and variants with uncertainty are removed.
+ clinvar benigns are either `Benign` or `Likely-Benign` and variants with uncertainty are removed.

### samocha

These are from [Kaitlin Samocha's paper](http://www.biorxiv.org/content/early/2017/06/12/148353) on mis-sense contraint.

+ benigns are labelled as `control` in her source file
+ pathogenics are anything other than control.


