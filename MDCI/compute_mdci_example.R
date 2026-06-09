###
# Compute MDCI 
# By Marcus Westerberg
# Change INSERT PATH and INSERT NAME to the appropriate folder / file name

# load scripts
source("INSERT PATH\\compute_mdci_function.R")
source("INSERT PATH\\mdci_inner.R")
source("INSERT PATH\\process_data.R")

# run code
compute_mdci(from_path="INSERT PATH\\", # path to par_parts (see readme file)
             weights_path="INSERT PATH\\mdci_weights.Rdata", # path to weights
             to_path="INSERT PATH", # path to folder where output will be saved
             ncores=4, # number of cores
             cohort_path="INSERT PATH", # path to folder where the cohort is stored as Rdata file
             cohort_file_name="INSERT NAME", # name of the cohort file (excluding file extension)
             cohort_colnames=c("LopNr","indexdate"), # colnames to be used in cohort file (first is the unique identifier and second is a date)
             mdci_inner=mdci_inner,
             cluster_ids=1:32,  # 1:X where X is the modulo number
             process_data=process_data)




###
# end
###