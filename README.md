## TB Gene Signature Pipeline

This repository contains a collection of scripts and notebooks working towards the development of a gene signature model predictive of stage of TB disease progression.


### Notebooks

The `notebooks/` directory contains examples of what the various components of an end-to-end pipeline might look like.

* `notebooks/machine-learning-pipeline.ipynb`: This notebook takes as input group phenotypes and normalized expression data from the ACS cohort. It first "screens" relevant features by either subsetting to genes from the network analysis or by performing a bootstrap stability analysis as described here: https://www.jmlr.org/papers/volume18/17-514/17-514.pdf. After pre-selecting relevant features, the notebook shows how sklearn's `Pipeline` object can be used in a nested cross-validation setup to 1) optimize hyperparameters in a set of different regressor models, and 2) select the best-performing regressor model among the optimized models.

* `notebooks/machine-learning-pipeline-smote.ipynb`: This notebook contains the same content as `notebooks/machine-learning-pipeline.ipynb`, but with an extra step at the beginning of the sklearn pipeline to oversample underrepresented outcome classes using the SMOTE algorithm (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-106).

* `notebooks/network-analysis-salmon.ipynb`: This notebook attempts to tie together the various components of the network analysis project, which are coded as functions in the `tb_gene_signature_pipeline/network_analysis.py` module. This version of the notebook uses bulk RNA-seq count matrices generated by Salmon.

* `notebooks/network-analysis-rsem.ipynb`: This notebook contains the same content as `notebooks/network-analysis-salmon.ipynb`, but with bulk RNA-seq count matrices generated by RSEM.

### tb_gene_signature_pipeline

* `tb_gene_signature_pipeline/network_analysis.py`: This script contains functions related to the network analysis.

* `tb_gene_signature_pipeline/plotting.py`: This script contains various plotting functions useful in plotting output from the network analysis.

