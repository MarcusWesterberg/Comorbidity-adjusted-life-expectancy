
require(dplyr)


pdr <- "INSERT CODE TO LOAD DRUG REGISTER"

# Final file should have the following columns: LopNr, ATC and EDATUM (as date YYYY-MM-DD)

# For parallel computation, save separate files like this

nparts <- 8 # number of parts (larger than or equal to the number of cores to be used in the parallel computations later)

for(j in 1:nparts){
  print(j)
  
  start <- Sys.time()

  drug_part <- pdr %>% 
    filter(LopNr %% nparts == (j-1))
  
  print(dim(drug_part))
  save(drug_part,file=paste0("INSERT PATH",j,".Rdata"))
  print(Sys.time()-start)
  rm(drug_part)
  gc()
}
