###
# Process data

process_data <- function(PAR_inner,
                         weights_data,
                         n_positions){

  predictor_names <- unique(weights_data$predictor_names[weights_data$npos==n_positions])

  codes <- unique(weights_data$code[weights_data$npos==n_positions])
  
  PAR_inner$diagnos <- substr(PAR_inner$diagnos,1,n_positions)
  
  # fix code
  fix_code <- function(x,n_positions){
    xnew <- x
    if(!is.na(x)){
      nx <- nchar(x)
      if(nx<n_positions){
        xnew <- paste0( c(x,rep(9,n_positions-nx)),collapse="" )
      }
    }
    return(xnew)
  }

  diagnos_codes <- unique(PAR_inner$diagnos)
  
  diagnos_codes_updated <- sapply(FUN=fix_code,
                                  diagnos_codes,
                                  n_positions=n_positions)
  
  diagnos_codes <- data.frame("diagnos"=diagnos_codes,
                              "diagnosu"=diagnos_codes_updated) 
  
  PAR_inner <- merge(PAR_inner,
               diagnos_codes,
               by="diagnos",
               all.x=FALSE)
  
  PAR_inner <- PAR_inner %>% 
    mutate(diagnos=diagnosu) %>%
    select(-c(diagnosu)) %>%
    filter(!is.na(diagnos))
  
  # not needed here, only during development of mdci
  # # remove codes not unique at this code depth relative to the short depth
  # PAR_inner$diagnos_short <- substr($diagnos,1,n_positions-1)
  # 
  # PAR <- PAR %>%
  #   group_by(diagnos_short) %>%
  #   mutate(n_dx = length(unique(diagnos)) ) %>%
  #   filter(n_dx>1) %>%
  #   ungroup() %>%
  #   select(-c(n_dx,diagnos_short))

  ###
  # Extract only codes from PAR present in weights_data for the corresponding length
  # important that this is here and not before the processing above
  
  PAR_inner <- PAR_inner %>%
    filter(diagnos %in% codes)
  
  unique_lopnr <- PAR_inner %>% 
    select(lopnr) %>% 
    unique() %>%
    arrange(lopnr)
                                       
  data_temp <- PAR_inner  %>%
    
    group_by(lopnr,diagnos) %>%
    
    mutate(diagnos_hdia=any(is_hdia),
           diagnos_time3=any(date-indatum<=365 & is_hdia), 
           diagnos_time2=any(date-indatum<=180 & is_hdia), 
           diagnos_time1=any(date-indatum<=90 & is_hdia),
           diagnos_number=length(unique(indatum[is_hdia])), 
           duration=sum(0,duration[is_hdia],na.rm=TRUE)
    ) %>%
    
    ungroup()   %>%
    rename(diagnos_duration=duration)   
  
  data_temp$diagnos_number1 <- as.numeric(data_temp$diagnos_number >= 2)
  data_temp$diagnos_number2 <- as.numeric(data_temp$diagnos_number >= 3)
  data_temp$diagnos_number3 <- as.numeric(data_temp$diagnos_number >= 4)
  data_temp$diagnos_duration1 <- as.numeric(data_temp$diagnos_duration >= 7 )
  data_temp$diagnos_duration2 <- as.numeric(data_temp$diagnos_duration >= 14)                                     
                       
  # for each code, create data
  template <- data_temp %>%
    select(lopnr,
           diagnos,
           diagnos_hdia,
           diagnos_number1,
           diagnos_number2,
           diagnos_number3,
           diagnos_duration1,
           diagnos_duration2,
           diagnos_time1,
           diagnos_time2,
           diagnos_time3) %>%
    unique()                                      
                                       
  data_list <- lapply(X=codes,
                      FUN = function(dx){
                        # Select only this code
                        d <- template %>%
                          filter(diagnos %in% dx) %>%
                          select(-diagnos)
                        
                        stopifnot(length(d$lopnr)==length(unique(d$lopnr)))
                        
                        colnames(d)[-1] <- paste0(dx,"_",colnames(d)[-1])
                        d[[dx]]<-1
                        
                        # all individuals 
                        d <- merge(unique_lopnr %>% select(lopnr),
                                   d,
                                   by="lopnr",
                                   all.x=TRUE) %>%
                          arrange(lopnr) %>%
                          select(-lopnr)
                        
                        d[is.na(d)]<-0
                        
                        return(d)
                      })      
  
  data_temp <- do.call(cbind,data_list)
  rm("data_list")
  rm("template")
  
  data_temp <- cbind(unique_lopnr,
                     data_temp)
  
  rm("unique_lopnr")                                    
  

  data_temp <- data_temp[,colnames(data_temp) %in% c("lopnr",predictor_names)]
  
  return(data_temp)
}


###
# End
###