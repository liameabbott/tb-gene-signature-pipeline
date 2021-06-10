## TB Gene Signature Pipeline

This repository contains a collection of scripts and notebooks working towards the development of a gene signature model predictive of stage of TB disease progression.


### Notebooks

The `notebooks/` directory contains examples of what the various components of an end-to-end pipeline might look like.

* `notebooks/machine-learning-pipeline.ipynb`: This notebook takes as input group phenotypes and normalized expression data from the ACS cohort. It first "screens" relevant features by either subsetting to genes from the network analysis or by performing a bootstrap stability analysis as described here: https://www.jmlr.org/papers/volume18/17-514/17-514.pdf. After pre-selecting relevant features, the notebook shows how sklearn's `Pipeline` object can be used in a nested cross-validation setup to 1) optimize hyperparameters in a set of different regressor models, and 2) select the best-performing regressor model among the optimized models.

