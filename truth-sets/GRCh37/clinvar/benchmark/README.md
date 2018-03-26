These sets are for benchmark testing metrics more fairly -- although some already utilize gnomAD in their metric creation.  The idea is that because many ClinVar benigns are variants found at 5% or greater frequency in the 1000G, ESP6500, ExAC, and gnomAD datasets, incorporating more variants from gnomAD to obtain a better idea of how a metric calls true negatives is the next best estimate. Thus we created a version of ClinVar benigns that incorporates gnomAD variants to match the count of ClinVar pathogenics.

Info on the three sets of benigns:

+ ClinVar ALL uses random variants from the whole gnomAD set to add to the set of ClinVar benigns so that the pathogenic count matches the benign count of variants.
+ ClinVar singleton uses random variants from the whole gnomAD set to add to the set of ClinVar benigns so that the pathogenic count matches the benign count of variants.
+ ClinVar AF < 0.005 uses random variants from the gnomAD set < 0.005 AF to add to the set of ClinVar benigns so that the pathogenic count matches the benign count of variants.
