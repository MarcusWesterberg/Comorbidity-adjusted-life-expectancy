
mdci_inner <- function(PAR_inner,
                       weights_data){
  
  
  
  ###   ###   ###   ###   ###   ###  ###   ###   ###   ###
  # generate prognostic factors for each code length
  prognostic_factors <- list()
  library(tidyr)
  
  for(j in 2:5){
    print(paste0("Code depth ",j))
    prognostic_factors[[j]] <- process_data(PAR_inner=PAR_inner,
                                            weights_data=weights_data,
                                            n_positions=j) %>% 
      gather(predictor_names, value, -c(lopnr))
  }
  rm("PAR_inner")
  gc()

  prognostic_factors <- do.call(rbind,
                                prognostic_factors) %>% arrange(lopnr)
  

  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
  # add corresponding weights
  print("Merging coefficients")
  mdci_data <- merge(prognostic_factors,
                     weights_data %>% select(predictor_names,coefficients),
                     by="predictor_names",
                     all.x=TRUE,
                     all.y=FALSE)
  
  print("Summing")
  mdci_data <- mdci_data %>%
    group_by(lopnr) %>%
    mutate(MDCI=sum(value*coefficients)) %>%
    select(lopnr,MDCI) %>%
    #unique() %>%
    group_by(lopnr) %>%
    slice(1) %>%
    ungroup() %>% 
    arrange(lopnr)
  
  return(mdci_data)
}