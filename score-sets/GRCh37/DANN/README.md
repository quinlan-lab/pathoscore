#DANN: a deep learning approach for annotating the pathogenicity of genetic variants.

DANN uses the same feature set and training data as CADD to train a deep neural network (DNN). DNNs can capture non-linear relationships among features and are better suited than SVMs for problems with a large number of samples and features. We exploit Compute Unified Device Architecture-compatible graphics processing units and deep learning techniques such as dropout and momentum training to accelerate the DNN training. DANN achieves about a 19% relative reduction in the error rate and about a 14% relative increase in the area under the curve (AUC) metric over CADD's SVM methodology.

Reference:https://www.ncbi.nlm.nih.gov/pubmed/25338716
