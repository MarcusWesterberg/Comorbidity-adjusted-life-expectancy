# Comorbidity-adjusted-life-expectancy
R code to compute the MDCI, DCI and health adjusted-life expectancy By Marcus Westerberg, Department of Surgical Sciences, Uppsala University, 2026

This document contains code to compute comorbidity indices and to compute health-adjusted life expectancy.
How to compute the comorbidity indices

See separate readme file "How to compute MDCI and DCI.txt" for how to compute the MDCI (DOI: 10.1371/journal.pone.0296804) and DCI (DOI: 10.1097/EDE.0000000000001358). Code for CCI used in this paper can be found here: https://github.com/bjoroeKI/Charlson-comorbidity-index-revisited.
How to compute life expectancy using the code underlying the paper?

See separate readme.txt in /Analysis
Fast computation of life expectancy using a compressed version of the method based on a GAM model (requires R package mgcv)

See example fastCALE_example.R script
