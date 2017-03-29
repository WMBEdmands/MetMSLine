MetMSLine
=========

[![DOI](https://zenodo.org/badge/21719/WMBEdmands/MetMSLine.svg)](https://zenodo.org/badge/latestdoi/21719/WMBEdmands/MetMSLine)

R functions for automation of biomarker discovery based on processing downstream of peak picking softwares.

Overview
===============

The workflow consists of 5 stages:
 
1. pre-processing. Performs all muMetMSLine
=========

[![DOI](https://zenodo.org/badge/14743752.svg)](https://zenodo.org/badge/latestdoi/14743752) latest stable release v1.1.0 (archived on zenodo repository). 

R functions for automation of biomarker discovery based on processing downstream of peak picking softwares.

![workflow_illustr](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4341062/bin/btu705f1p.jpg)

The MetMSLine workflow scripts and associated functions have now migrated to package form. The first 3 data processing scripts as discussed in the publication (*PreProc.QC.LSC.R*, *Auto.PCA.R*, *Auto.MV.Regress.R*) are incorporated within this package and explained in the two accompanying vignettes. The 4th and 5th scripts (*Auto.MS.MS.match.R*, *DBAnnotate.R*) which are concerned with metabolite annotation have now been largely supplanted by the [compMS2Miner](https://github.com/WMBEdmands/compMS2Miner) package. The original scripts discussed in the publication can still be found in the [MetMSLine_Scripts](https://github.com/WMBEdmands/MetMSLine_Scripts) repository.

MetMSLine combined with the [compMS2Miner](https://github.com/WMBEdmands/compMS2Miner) package are intended to facilitate autonomous metabolomic data analysis. These two packages combined with the [xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html) and [CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html) R packages a complete and largely autonomous metabolomic workflow can be achieved.

If you find MetMSLine useful please cite us:

**MetMSLine: an automated and fully integrated pipeline for rapid processing of high-resolution LC-MS metabolomic datasets.**
*William Matthew Bell Edmands, Dinesh Kumar Barupal, Augustin Scalbert*
Bioinformatics 2015; 31 (5): 788-790.
[DOI: 10.1093/bioinformatics/btu705](https://doi.org/10.1093/bioinformatics/btu705) 

Overview
===============

The workflow consists of 4 stages:

1. pre-processing. Performs all multiparametric preprocessing steps for large-scale high-resolution LC-MS metabolomic datasets.

2. PCA outlier removal and cluster identification. Principal Component Analysis (PCA) of pre-processed LC-MS data with iterative automatic outlier removal based on a user defined Hotelling’s T2 ellipse expansion
and PCA scores cluster identification (using PAM clustering and regression).

3. Objective univariate statistical analysis based on covariate type.
Multiparametric, automatic regression/statistical analysis, biomarker discovery for high resolution LC-MS data with multiple Y variables. 

4. Concerned with unknown LC-MS feature annotation has now been largely supplanted by the [compMS2Miner](https://github.com/WMBEdmands/compMS2Miner) package.
...

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
ltiparametric preprocessing steps for large-scale high-resolution LC-MS metabolomic datasets.

2. PCA outlier removal and cluster identification. Principal Component Analysis (PCA) of pre-processed LC-MS data with iterative automatic outlier removal based on a user defined Hotelling’s T2 ellipse expansion
and PCA scores cluster identification (using PAM clustering and regression).

3. Objective univariate statistical analysis based on covariate type.
Multiparametric, automatic regression analysis, biomarker discovery, hierarchical clustering analysis and cluster ion and isotope identification for high resolution LC-MS data with multiple continuous Y variables. 

4. Biomarker identification...

5. ...

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
