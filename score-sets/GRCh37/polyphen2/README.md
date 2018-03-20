From the PolyPhen site:
"Two pairs of datasets were used to train and test PolyPhen-2 prediction models. The first pair, HumDiv, was compiled from all damaging alleles with known effects on the molecular function causing human Mendelian diseases, present in the UniProtKB database, together with differences between human proteins and their closely related mammalian homologs, assumed to be non-damaging. The second pair, HumVar, consisted of all human disease-causing mutations from UniProtKB, together with common human nsSNPs (MAF>1%) without annotated involvement in disease, which were treated as non-damaging.

Diagnostics of Mendelian diseases requires distinguishing mutations with drastic effects from all the remaining human variation, including abundant mildly deleterious alleles. Thus, the HumVar-trained model should be used for this task. In contrast, the HumDiv-trained model should be used for evaluating rare alleles at loci potentially involved in complex phenotypes, dense mapping of regions identified by genome-wide association studies, and analysis of natural selection from sequence data, where even mildly deleterious alleles must be treated as damaging."

These files are precomputed on the ClinVar benchmark and the Samocha truth sets.  There is some issue with the precomputed values derived from the Polyphen site, such as multiple values for the same allele at the same position for both model versions.

Paper for Polyphen-2: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4480630/
