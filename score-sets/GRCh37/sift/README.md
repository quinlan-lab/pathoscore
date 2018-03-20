SIFT uses the degree of protein sequence conservation to predict the impact of a missense variant. It uses alignment based method to give prediction scores in the range [0,1]. Values <= 0.05 are predicted to be deleterious, while those > 0.05 are predicted to be tolerated. The columns are "chrom  position ref alt siftScore".  
Paper at: http://genome.cshlp.org/content/11/5/863.full.pdf+html 

These scores are precomputed and derived from CADD's annotation collection and scored for the ClinVar benchmark and Samocha truth sets.
