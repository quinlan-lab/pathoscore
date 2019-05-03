##############
# ANNOTATION #
##############
bash annotate.sh

##############
# EVALUATION #
##############

# Table 3
bash scangenes.sh # gets top 10 ClinVar/Samocha Genes

# Figures 1, 2, 5, 6; Table 3 (default list is top 10 ClinVar genes according to Yandell lab; we used top 10 in our ClinVar patho file, and top 10 list of Samocha set was manually curated via neurodev set as well)
bash evaluate.sh

# Table 1 (or results for cutoff section)
# to get this one, just used the results of evaluate.sh and looked by hand through ClinVar. no real script.
# for MetaSVM data...
python checkcutoff.py pathogenic.vcf.gz benign.vcf.gz

# Table 2
# grabbed them from the results of evaluate.sh by hand.  no script, parsing the HTML was too difficult, may want to use plot.ly python library in future releases (used unfiltered sets for all scores)

# Bar plots (Figures 4 and 5)
bash scorebarplots.sh # Figure 4
bash cubarplots.sh # Figure 5

###############
# CORRELATION #
###############
# Figure 5
bash correlation.sh
