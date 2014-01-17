MetMSLine
=========

R functions for automation of biomarker discovery based on processing downstream of XCMS based peak picking

PreProc.QC.LSR.R
===============

Performs all multiparametric preprocessing steps for large-scale high-resolution LC-MS metabolomic datasets.

```
>PreProc.QC.RLSC(X="XCMS_output.tsv",CCQC=10,SGroups=3,QCInt=5,f=1/5,
RSD=30,a=1,scatter.plots=TRUE,wd="D:\\R_data_processing\\STUDY NAME\\",XCMS_dir="D:\\R_data_processing\\STUDY NAME\\XCMS\\")
```

Arguments: -
--------------------

1.	**X** - XCMS diffreport name in inverted commas saved as either a .tsv, .txt or .csv (tab separated values, text or comma seperated values) file
2.	**CCQC** - number of column conditioning QCs (default is 30).
3.	**QCInt** - QC injection interval (default is 10 or after every 9th sample). 
4.	**SGroups** � number of study groups in the sample ie. mzXML file subdirectories prior to XCMS peak-picking (defaults to 3 a good subdirectory structure is to have QCs, samples and blanks in 3 seperate subdirectories).
5.	**f** - the smoother span (proportion of points in the plot which influence the smooth at each value)(default is 1/5). 
6.	**a** - alpha for generalized log transform (default is 1).
7.	**RSD** - percentage relative standard deviation cutoff for reproducible features assessed from the pooled QC samples following signal smoothing (default is 30% RSD).
8.	**wd** - name of study parent working directory in inverted commas (default="D:\\R_data_processing\\STUDY NAME\\").
9.	**XCMS_dir** - working directory where XCMS output file is located in inverted commas (default="D:\\R_data_processing\\STUDY NAME\\XCMS\\").


Auto.PCA.R
==========

Principle Components Analysis of pre-processed LC-MS data with iterative (up to 2 rounds if necessary) automatic outlier removal based on a user defined Hotelling�s T2 ellipse expansion.

```
>Auto.PCA(X="QC_Corrected.csv",Yvar="Y.csv",QCInt=5,PCA.method="svd",scaling="uv",out.tol=1.2,wd="D:\\R_data_processing\\STUDY NAME\\")
```

Arguments: -
----------------------

1.	**X** - QC corrected output data from PreProc.QC.LSC.R (default = "QC_Corrected.csv")
2.	**Yvar** - Y variable matrix in the layout columns = Y variable names, rows = X matrix sample names (default = "Y.csv")
3.	**QCInt** - QC samples interval for your experiment (default = 5 ie. every 4th sample)
4.	scaling - scaling method for PCA calculation (default = "uv", other options include "none", "vector" and "uv")
5.	**PCA.method** � PCA calculation method (default = "svd", other options include �ppca�, �nipals�, �rnipals�, �bpca� and �svdImpute�)
6.	**out.tol** - outlier removal tolerance based on an proportional expansion of the Hotellings T2 ellipse beyond which to remove outliers (default =1.2)
7.	**wd** - name of study parent working directory in inverted commas (default="D:\\R_data_processing\\STUDY NAME\\"). The function will automatically go to the previous PreProc.QC.RLSC.results folder within the study so the user does not need to state this explicitly.


Auto.MV.Regress.R
================

Multiparametric, automatic regression analysis, biomarker discovery, hierarchical clustering analysis and cluster ion and isotope identification for high resolution LC-MS data with multiple continuous Y variables. 

```
> Auto.MV.Regress( X="PCA.outliers.removed.csv", Yvar="Y.outliers.removed.csv", box.prop=1/5, wd="D:\\R_data_processing\\STUDY NAME\\", p.adjust.methods="none", heatmap=TRUE, hclust.method="complete", dist.m.method="euclidean", Yunits="Y_units", Corr.thresh=0.3, pvalue=0.01, HMDBtol=0.005, Clust.RT.tol=2, mode="negative", Clust.ppm=10, non.zero=2) 
```

Arguments: -
--------------------

1.	**X** � Mass spectral output dataset from the previous Auto.PCA.R function, with outliers removed (default = "PCA.outliers.removed.csv").
2.	**Yvar** - Y variable matrix output from the previous Auto.PCA.R function, with outliers removed (default = "Y.outliers.removed.csv").
3.	**box.prop** � proportion of total sample number to use for high versus low box plots and score calculation (default = 1/5 ie. This will take quintiles of both highest and lowest Y variables).
4.	**wd** - name of study parent working directory in inverted commas (default="D:\\R_data_processing\\STUDY NAME\\"). The function will automatically go to the previous Auto.PCA.results folder within the study so the user does not need to state this explicitly.
5.	**p.adjust.methods** � p value multiple testing correction method, defaults to �BH� Benjamini-Hochberg or false discovery rate calculation as it is less stringent than other methods for biomarker discovery however other methods are optional (�"holm", "hochberg", "hommel", "bonferroni", "BY","fdr", "none").  
6.	**heatmap** � TRUE / FALSE condition for hierarchical clustering analysis of all potential biomarkers of all Y-variables supplied above the user defined threshold and both X-Y and X-X heatmap creation (default is TRUE)  
7.	**hclust.method** � hierarchical clustering agglomeration method (default = �complete� however other options include "ward", "single", "average", "mcquitty", "median" or "centroid").
8.	**dist.m.method** � distance to model measure (default = "euclidean", however other options include "maximum", "manhattan", "canberra", "binary" or "minkowski").
9.	**Yunits** � the names of the units of your continuous Y variables. This unit name will appear on the saved scatter and box and whisker plots of potential biomarkers identified.
10.	**Corr.thresh** � pearson product moment correlation threshold for biomarker selection. (default =0.3, although this can be set higher with greater risk of false negatives).
11.	**pvalue** � p value calculated from an F hypothesis test statistic and multiple testing method corrected threshold (default=0.01). 
12.	**HMDBtol** � mass accuracy tolerance for Human Metabolome Database (HMDB) search, a link to the sites search engine is provided in the results table (default=0.005).
13.	**Clust.RT.tol** � time in seconds to identify cluster ions in biomarkers above threshold following hierarchical clustering (default= 2).
14.	**mode** � mass spectrometer ionisation mode either positive or negative (NB. Must be lower case).
15.	**Clust.ppm** � mass accuracy tolerance for cluster ion and isotope identification (default = 10). 
16.	**non.zero** � minimum percentage of non-zero datapoints in a Y-variable for it to be utilised for biomarker discovery, this allows Auto.MV.regress to filter out continuous Y variables with low sample number (default= 2% eg. if number samples = 500 then 2% = 10 samples minimum).


Auto.MS.MS.match.R
==================

Automated biomarker � experimental MS/MS fragmentation spectra matching function based upon retention time and mass accuracy tolerance windows. The function also annotates possible common neutral losses and fragments based upon high mass accuracy MS/MS fragmentation data and outputs figures and results tables. The function was designed primarily for data-dependent MS/MS fragmentation method data although would in principle work equally for targeted MS/MS (.mzXML) data files.

```
> Auto.MS.MS.match(MSfeatures="Features_Above_threshold.csv",  mode="negative", wd="D:\\R_data_processing\\STUDY NAME\\",  mzXML.dir="D:\\R_data_processing\\STUDY NAME\\MS_MS_mzXML\\", TICfilter=5000,Precursor.ppm=10,Frag.ppm=20,ret=5, Parent.tol=0.1, Fragment.tol=0.5) 
```

Arguments: -
-------------------

1.	**MSfeatures** � Significant features identified by the previous Auto.MV.Regress.R function or the output of the DBAnnotate.R function, however it is flexible enough to work with other data files (.csv) originating from XCMS. (default = �Features_Above_threshold.csv ").
2.	**mode** � mass spectrometer ionisation mode either positive or negative (NB. Must be lower case).
3.	**wd** - the address of your parent study directory (�STUDY NAME�). The function will automatically go to the previous Auto.MV.Regress.results folder so the user does not need to state this explicitly, therefore do not move the result from this directory
4.	**mzXML.dir** � the directory address of the location of the MS/MS .mzXML files with which to match against potential biomarkers identified by the Auto.MV.Regress.R function. If the Auto.Dir.file.sorter.R has been used or the standardised subdirectory layout has been followed then the address would be in the text string form �D:\\R_data_processing\\STUDY NAME\\MS_MS_mzXML\\".
5.	**TICfilter** � total ion current filter, this filter is applied to ensure that high quality MS/MS spectra with a total ion current above the value are used to match to the unknown biomarkers (default = 5000).
6.	**Precursor.ppm** � mass accuracy in parts per million with which to match potential biomarkers to precursor ions of MS/MS data (default =10)
7.	**Frag.ppm** - mass accuracy in parts per million with which to match fragments and neutral losses in MS/MS data matches to common fragment masses, adducts and conjugates. The Frag.ppm argument should be set higher as with data-dependent MS/MS algorithms there can be a loss of mass accuracy compared to MS mode (default =20).
8.	**ret** � retention time window � in seconds with which to match potential biomarker mass spectral features to MS/MS data (default = �2 seconds).
9.	**Parent.tol** � delta mass accuracy window for parent/precursor ion for matching of unknown fragmentation spectra to the Human Metabolome Database (HMDB) experimental MS/MS database. (default = �0.1 m/z). 
10.	**Fragment.tol** - delta mass accuracy window for top six fragment ions (ranked by relative intensity) for matching of unknown fragmentation spectra to the Human Metabolome Database (HMDB) experimental MS/MS database. (default = �0.5 m/z). A broad tolerance is necessary as experimental MS/MS spectra is mainly acquired using nominal mass with less mass accuracy.


DBAnnotate.R
===========

Automated putative biomarker annotation function based upon high mass accuracy tolerance matching. A user supplied list of biologically plausible/selected metabolites is used to generate all iterations of common Phase II conjugate modifications and electrospray adducts in a �shotgun� approach. The combination of this �shotgun� method and the use of a narrowed list of biologically plausible metabolites represent a balance between the �false negative� and �false positive� rates of biomarker annotation and metabolite prediction. Following completion the results table must be visually inspected to decide on plausibility of subsequent biomarker identifications. 

```
>   DBAnnotate( X="Features_above_threshold.csv", database="metabolite_DB.csv", mode="negative", conjugates=c("Gluc","Sulf","DiGluc","DiSulf","GlucSulf","NAcCys"), MassAcc=10, wd="D:\\R_data_processing\\STUDY NAME\\",unknowns.dir=�"D:\\R_data_processing\\STUDY NAME\\Auto.MV.Regress.results\\")
```

Arguments: -
-------------------

1.	**X** � Significant features identified by the previous Auto.MV.Regress.R or the Auto.MS.MS.match.R function, however it is flexible enough to work with other data files (.csv) originating from XCMS. (default = �Features_Above_threshold.csv").
2.	**database** � an experimental-specific semi-targeted list of biologically plausible/expected metabolites (.csv), the only essential requirement for the database table is that the name of the compound is in the first column table and the monoisotopic/theoretical masses are in the second table column. Additional information columns adjacent to the first two columns will also be incorporated in the result table output by the function and not lost (default = �metabolites_DB.csv�).
3.	**mode** � mass spectrometer ionisation mode either positive or negative (NB. Must be lower case).
4.	**Conjugates** � a list (not exhaustive) of potential phase II modifications of metabolites in human biofluids. Included are glucuronides (�Gluc�), sulfates (�Sulf�), diglucuronides (�DiGluc�), disulfates (�DiSulf�), mixed glucuronide sulfates (�GlucSulf�) and N-acetyl cysteine conjugates (�NAcCys�). By default all conjugates are considered however to remove conjugates simply delete entries from the list (eg. conjugates=c(�Sulf�, �NAcCys�)).
5.	**MassAcc** � mass accuracy tolerance in parts per million for biomarker annotation dependent on mass spectrometer performance (default = 10 ppm).
6.	**wd** - name of study parent working directory in inverted commas (default="D:\\R_data_processing\\STUDY NAME\\").
7.	**unknowns.dir** - the directory address of your Auto.MV.Regress.results directory or MS/MS fragmentation data file folder (e.g. �MS_MS_mzXML�). 
