# Code for: The association between muscle quantity and overall survival depends on muscle radiodensity: a cohort study in non-small cell lung cancer

Newest version on [github](https://github.com/vanAmsterdam/bodycomp-survival-public/)

## How to reproduce

- most of the code assumes that there are data files in a subdirectory called 'data', these cannot be provided here due to privacy regulations
- all the code is in the 'code' directory

1. `01_merge_data.R` merge clinical data with measurements, selecting the right scan dates
2. `02_imputation.R` run imputation (can take some ~24 hours)
3. `03_analysis.R` this contains all the necessary analyses to reproduce the results in the paper
4. `power.R` this contains code for the power analysis. This is entirely self-contained and can be run without access to the data.

## Other files

`code/00_functions.R` and `code/pgwsource.R` contain function definitions that are used in the analysis.
The `pgwsource.R` file was adapted from the original publication:

Burke K, Jones MC, Noufaily A. A Flexible Parametric Modelling Framework for Survival Analysis. Journal of the Royal Statistical Society: Series C (Applied Statistics). 2020 Apr;69(2):429–57. 
https://doi.org/10.1111/rssc.12398p

The resource directory contains some helper files.
The package directory contains the sources of two packages that were amended for this analysis: `survival` (added support for mixed-stratification) and `smcfcs` (added support for stratified cox models and some steps to prevent numerical underflow).
