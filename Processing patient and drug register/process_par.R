###
# Optional code to preprocess the patient register (par)

# load par
par <- "INSERT CODE TO LOAD"

require(dplyr)

pasteme <- function(x){
  x <- x[!x %in% ""]
  
  if(length(x)==0){
    x <- ""
  }
  return(paste0(x,collapse=" "))
}

par <- par %>%
  rowwise() %>%
  mutate(diagnos=pasteme(c(dia1,dia2,dia3,dia4,dia5,dia6,dia7,dia8,dia9,
                           dia10,dia11,dia12,dia13,dia14,dia15,dia16,dia17,dia18,dia19,dia20,
                           dia21,dia22,dia23,dia24,dia25,dia26,dia27,dia28,dia29,dia30))) %>%
  ungroup()


# Final file should have the following columns: lopnr,indatum,utdatum,hdia,diagnos,source
# where indatum and utdatum = Date (YYYY-MM-DD), and utdatum = NA for outpatient visits and source = "ov" / "sv" for out / inpatient visits

# For parallel computation, save separate files like this

nparts <- 8 # number of parts (larger than or equal to the number of cores to be used in the parallel computations later)
for(j in 1:nparts){
  par_temp <- par %>%
    filter(lopnr %% nparts == (j-1))
  
  save(par_temp,file=paste0("INSERT PATH\\par_part",j,".Rdata"))
}
