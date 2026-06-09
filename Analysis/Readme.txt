Code for 

Data processing: in folder "Data processing"
Plotting: in folder "Plotting"
Descriptive tables and figures: in folder "Descriptive analyses"

Analysis order: 

Requires initial data management to be completed, then

1. Run Data processing / Main.R 

2. Compute C-indices: compute_cindex.R

3. Life expectancy
	- life_expectancy_estimate_lifetables.R -- Life table based Estimates from Statistics Sweden
	- life_expectancy_model_fit.R -- FIT the parametric models
	- life_expectancy_model_calibration.R -- assess calibration

	- life_expectancy_estimate.R -- estimate life expectancy using MDCI and DCI
	- life_expectancy_estimate_ccidci.R -- estimate life expectancy using CCI and DCI

	- life_expectancy_mdci_dci_table.R -- for GAM estimate of life expectancy for fast computation
	- life_expectancy_cci_dci_table.R -- for GAM estimate of life expectancy for fast computation



