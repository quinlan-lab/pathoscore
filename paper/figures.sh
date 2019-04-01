##############
# ANNOTATION #
##############
bash annotate.sh

##############
# EVALUATION #
##############

# Table 3
bash scangenes.sh # gets top 10 ClinVar/Samocha Genes

# Figures 1, 2, 3, 4; Table 3 (default list is top 10 ClinVar genes according to Yandell lab; we used top 10 in our ClinVar patho file, and top 10 list of Samocha set was manually curated via neurodev set as well)
bash evaluate.sh

# Table 1
# to get this one, just used the results of evaluate.sh and looked by hand through ClinVar. no real script.

# Table 2
# grabbed them from the results of evaluate.sh by hand.  no script, parsing the HTML was too difficult, may want to use plot.ly python library in future releases (used unfiltered sets for all scores)


###############
# CORRELATION #
###############
# Figure 5
bash correlation.sh
