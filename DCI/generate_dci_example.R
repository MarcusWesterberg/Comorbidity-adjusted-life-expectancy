###
# Generate DCI
# By Marcus Westerberg
# Change INSERT PATH and INSERT NAME to the appropriate folder / file name


path_weights_men <- "INSERT PATH\\Men.RData"
path_weights_women <- "INSERT PATH\\Women.RData"

path_lm_reg <- "INSERT PATH"  # path to the files drug_partx.Rdata (see readme file)

to_path <- "INSERT PATH" # output folder where the DCI data should be saved

cohort_path <- "INSERT PATH" # path to the folder that contains the Rdata set containing the cohort (LopNr and date variable and sex (Kon))
cohort_file_name <- "INSERT NAME" # insert name of the file (exclude file ending and index)
cohort_colnames <- c("LopNr","indexdate","Kon") # update if required

source("INSERT PATH\\compute_dci_function.R") # path to script that performed the computations

compute_dci(path_weights_men=path_weights_men,
            path_weights_women=path_weights_women,
            path_lm_reg=path_lm_reg,
            to_path=to_path,
            cohort_path=cohort_path,
            cohort_colnames=cohort_colnames,
            cohort_file_name=cohort_file_name, 
            cluster_ids=1:32, # 1:X where X is the modulo number
            ncores=4) # number of available cores



###
# End
###