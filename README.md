MetMSLine
=========

[![DOI](https://zenodo.org/badge/21719/WMBEdmands/MetMSLine.svg)](https://zenodo.org/badge/latestdoi/21719/WMBEdmands/MetMSLine)

R functions for automation of biomarker discovery based on processing downstream of peak picking softwares.

Overview
===============

The workflow consists of 5 stages:

1. pre-processing. Performs all multiparametric preprocessing steps for large-scale high-resolution LC-MS metabolomic datasets.

2. PCA outlier removal and cluster identification. Principal Component Analysis (PCA) of pre-processed LC-MS data with iterative automatic outlier removal based on a user defined Hotelling’s T2 ellipse expansion
and PCA scores cluster identification (using PAM clustering and regression).

3. Objective univariate statistical analysis based on covariate type.
Multiparametric, automatic regression analysis, biomarker discovery, hierarchical clustering analysis and cluster ion and isotope identification for high resolution LC-MS data with multiple continuous Y variables. 

Installation
==============
install directly from github using the devtools package. First install devtools,
instructions can be found here: https://github.com/hadley/devtools

```
library(devtools)

install_github('WMBEdmands/MetMSLine')
```

Getting started
===============

After MetMSLine installation has completed read the package vignette *"MetMSLineBasics"*
for some tips on getting started and also the vignette *MetMSLine* for further
examples of the package functions. 
Just type ```vignette('MetMSLineBasics')``` or ```vignette('MetMSLine')``` to view
the pdf versions of the vignettes.
The package also includes some example data to illustrate the workflow. 


licence
=============
The MetMSLine package is licenced under the GPLv3 (http://www.gnu.org/licenses/gpl.html).
