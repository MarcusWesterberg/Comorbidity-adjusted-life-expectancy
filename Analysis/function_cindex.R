###
# Compute C indices and pool

# with/without internal bootstrap per core


resample_individual <- function(data,
                                N=NULL){
  
  lpnr <- unique(data$LopNr)
  
  
  sampled_lpnr <- sample(lpnr,replace=TRUE,size=N) 
  
  count <- table(sampled_lpnr)
  count <- data.frame("LopNr"=as.numeric(names(count)),
                      "n_boot_repl"=as.numeric(count))
  
  return(count)
}

generate_bootstrapped_data <- function(data,count){
  
  data <- merge(data,
                count,by="LopNr")
  
  data <- data %>% tidyr::uncount(n_boot_repl,
                                  .id="boot_id")
  
  return(data)
}


# computes c-index and returns a data frame with the c-index, sample_size and potentially other things
helper <- function(data,
                   at_time,
                   weight=NULL,
                   boot_resample=FALSE){
  
  ret <- data.frame("timepoint"=at_time,
                    "N"=nrow(data),
                    "Nunique"=length(unique(data$LopNr)),
                    "cindex"=NA,
                    "cindex_var"=NA,
                    "concordant"=NA,
                    "discordant"=NA,
                    "tiedx"=NA,
                    "tiedy"=NA,
                    "tiedxy"=NA,
                    "computetime"=NA)
  
  
  require(survival)
  start <- Sys.time()
  
  if(!is.null(weight)){
    weight <- data[[weight]]
  }
  
  if(nrow(data)>0){
    
    if(boot_resample){
      stop("Bootstrap not supported")
      
    }
    
    # prepare for evaluation at at_time
    data$censor[data$timefu>at_time]<-0
    data$timefu[data$timefu>at_time]<-at_time
    
    cret <- tryCatch(
      {
        concordance(Surv(timefu, censor) ~ comorbidity_index,
                            data = data,
                            reverse = TRUE,
                            weights=weight,
                            cluster=data$LopNr)
      },
      error = function(cond) {
        message(paste("Error in concordance."))
        message("Here's the original error message:")
        message(conditionMessage(cond))
        # Choose a return value in case of error
        NA
      }
    )
    if(class(cret) %in% "concordance"){
      ret$cindex <- cret$concordance
      ret$cindex_var <- cret$var
      ret$concordant <- cret$count[1]
      ret$discordant <- cret$count[2]
      ret$tiedx <-cret$count[3]
      ret$tiedy <- cret$count[4]
      ret$tiedxy <- cret$count[5]
      ret$computetime <- as.numeric(Sys.time()-start,units="secs")
    } 
  }
  
  return(ret)
}
# End helper

###
# cindex package and function
compute_cindex_inner <- function(data,
                                 subgroup_function,
                                 subgroups,
                                 at_times,
                                 weight=NULL,
                                 boot_resample=FALSE){
  
  ret_inner <- list()
  total_iter <- length(at_times)*length(subgroups)
  index <- 1
  
  for(j in 1:length(subgroups)){
    
    data_temp <- subgroup_function(data,
                                   subset_g=subgroups[[j]]$subgroup,
                                   subset2_g=subgroups[[j]]$subgroup2,
                                   indexyear_g=subgroups[[j]]$indexyear,
                                   indexage_g=subgroups[[j]]$indexage,
                                   sex_g=subgroups[[j]]$sex,
                                   charlson_cat_g=subgroups[[j]]$charlson_cat,
                                   charlson_g=subgroups[[j]]$charlson,
                                   icd_n_chapters_g=subgroups[[j]]$icd_n_chapters,
                                   icd_chapters_g=subgroups[[j]]$icd_chapters,
                                   any_par_g=subgroups[[j]]$any_par)  
    
 
    
    
    ret_temp <- list()
    for(i in 1:length(at_times)){
    
      if(index %% 10 == 0 ){
        cat("\n")
        print(Sys.time())
        cat("\n")
      }
  
        print(paste0("Progress: ",index," out of ",total_iter,". Nrows=",nrow(data_temp)))
   
      ret_temp[[i]] <- helper(data=data_temp,
                              at_time=at_times[i],
                              weight=weight,
                              boot_resample=boot_resample)
      
      ret_temp[[i]]$subgroup_index <- j
      ret_temp[[i]]$subgroup <- paste0(subgroups[[j]]$subgroup,collapse=",")
      ret_temp[[i]]$subgroup2 <- paste0(subgroups[[j]]$subgroup2,collapse=",")
      ret_temp[[i]]$indexyear <- paste0(subgroups[[j]]$indexyear,collapse=",")
      ret_temp[[i]]$indexage <- paste0(subgroups[[j]]$indexage,collapse="-")
      ret_temp[[i]]$charlson <- paste0(subgroups[[j]]$charlson,collapse=",")
      ret_temp[[i]]$charlson_cat <- paste0(subgroups[[j]]$charlson_cat,collapse=",")
      ret_temp[[i]]$icd_n_chapters <- paste0(subgroups[[j]]$icd_n_chapters,collapse=",")
      ret_temp[[i]]$icd_chapters <- paste0(subgroups[[j]]$icd_chapters,collapse=",")
      
      index <- index + 1
    }
    ret_inner[[j]] <- do.call(rbind,ret_temp)
  }
  ret_inner <- do.call(rbind,ret_inner)
  
  return(ret_inner)
}








compute_cindex <- function(ncores,
                           cluster_indices,
                           log_path,
                           data_path,
                           subgroups,
                           subgroup_function,
                           at_times,
                           save_parts=TRUE,
                           save_path=NULL,
                           save_name=NULL,
                           risk_score_name,
                           neg_to_zero=rep(FALSE,length(risk_score_name))){
  
  
  
  if(save_parts){
    stopifnot(!is.null(save_path))
    stopifnot(!is.null(save_name))
  }

  
  stopifnot(!is.null(names(subgroups)))
  
  
  print("Starting...")
  
  require(dplyr)
  require(parallel)
  
  start_time <- Sys.time()
  
  print(cluster_indices)
  ncores <- min(ncores,length(cluster_indices))
  
  
  print(paste0("Ncores=",ncores))
  
  cl <- makeCluster(ncores)
  
  clusterExport(cl, c("log_path",
                      "data_path",
                      "subgroups",
                      "subgroup_function",
                      "at_times",
                      "save_parts",
                      "save_path",
                      "save_name",
                      "risk_score_name",
                      "neg_to_zero"),envir = environment())
  clusterExport(cl, c("compute_cindex_inner", "helper"),envir = .GlobalEnv)
  
  ###
  # Load data first and then iterate over subgroups
  ret <- parLapply(cl = cl,
                   cluster_indices,
                   fun = function(cluster_id){
                     
                     
                     cluster_id <<- cluster_id # make seed global
                     set.seed(cluster_id)
                     start_time <- Sys.time()
                     
                     require(tidyverse)
                     require(survival)
                     sink(file=paste0(log_path,"cindex_log_",save_name,"_",cluster_id,".txt"))
                     print(cluster_id)
                     print(Sys.time())
                     
                     cat("\n ... \n")
                     print(log_path)
                     print(at_times)
                     print(save_parts)
                     print(save_path)
                     print(save_name)
                     print(risk_score_name)
                     cat("\n ... \n")
                     
                     
                     print("Loading data...")
                     Sys.sleep(cluster_id/4)
                     
                     loads <- load(file=paste0(data_path,"cindex_cache_part",cluster_id,".Rdata"))
                     
                     ###
                     # keep only needed coluns
                     data <- data %>%
                       select(LopNr,rowid,indexdate,indexyear,subset,subset2,timefu,censor,sex,indexage,
                              starts_with("comorbidity_index"))
                
                     ###
                     # Drop stuff not needed
                     all_subgroups <- unique(unlist(lapply(FUN=function(x) x$subgroup ,subgroups)  ))
                     if(!is.null(all_subgroups)){
                       print("Filtering subgroup")
                       print(all_subgroups)
                       data <- data %>%
                         filter(subset %in% all_subgroups)
                     }
                     
                     all_subgroups2 <- unique(unlist(lapply(FUN=function(x) x$subgroup2 ,subgroups)  ))
                     if(!is.null(all_subgroups2)){
                       print("Filtering subgroup2")
                       print(all_subgroups2)
                       data <- data %>%
                         filter(subset2 %in% all_subgroups2)
                     }
                     
                     all_ages <- unique(unlist(lapply(FUN=function(x) x$indexage ,subgroups)  ))
                     if(!is.null(all_ages)){
                       print("Filtering ages")
                       print(all_ages)
                       data <- data %>%
                         filter(indexage >= (min(all_ages)-1) & indexage <= (max(all_ages)+1)  )
                     }
                     
                     all_years <- unique(unlist(lapply(FUN=function(x) x$indexyear ,subgroups)  ))
                     if(!is.null(all_years)){
                       print("Filtering years")
                       print(all_years)
                       data <- data %>%
                         filter(indexyear %in% all_years)
                     }
                     
                     all_cci <- unique(unlist(lapply(FUN=function(x) x$charlson ,subgroups)  ))
                     if(!is.null(all_cci)){
                       print("Filtering CCI")
                       print(all_cci)
                       data <- data %>%
                         filter(comorbidity_index_CCI10 %in% all_cci)
                     }
                     
                     
                     ###
                     all_par <- unique(unlist(lapply(FUN=function(x) x$any_par ,subgroups)  ))
                     if(!is.null(all_par)){
                       print("Filtering PAR")
                       print(all_par)
                       
                       loads <- load(file=paste0("X:\\ComorbidityBase\\Works\\Derived variables\\PAR\\","par_part",cluster_id,".Rdata"))
                       par_part <- par_part %>%
                         filter(LopNr %in% data$LopNr) %>%
                         select(LopNr,indatum) %>%
                         unique()
                       
                       in_par <- rep(FALSE,nrow(data))
                       
                       urid <- unique(data$rowid)
                       
                       for(j in 1:length(urid)){
                         print(j)
                         
                         data_temp <- data %>%
                           filter(rowid %in% urid[j]) %>%
                           select(LopNr,indexdate)
                         
                         par_part_temp <- par_part %>%
                           left_join(data_temp,
                                     by="LopNr") %>%
                           group_by(LopNr) %>%
                           filter(any(indatum <= indexdate & indatum >= (indexdate - 3652.4))) %>%
                           ungroup()
                         
                         in_par[data$rowid %in% urid[j] ][ data_temp$LopNr %in% par_part_temp$LopNr ] <- TRUE
                       }
                       
                       rm(par_part_temp)
                       rm(data_temp)
                       
                       data$in_par <- as.character(in_par)
                       
                       data <- data %>%
                         filter(in_par %in% all_par)
                       
                       rm(par_part)
                       gc()
                     }
                     
                     
                     
                     all_icd_n_chapt <- unique(unlist(lapply(FUN=function(x) x$icd_n_chapters ,subgroups)  ))
                     all_icd_chapt <- unique(unlist(lapply(FUN=function(x) x$icd_chapters ,subgroups)  ))
                     
                     if(!is.null(all_icd_n_chapt) |!is.null(all_icd_chapt)){
                       print("Loading ICD chapters")
                       icd_val <- get(load(file=paste0(data_path,"validation_ICD_lookback_10_part_",cluster_id,".Rdata")))
                       rm(cdata)
                       
                       icd_dev <- get(load(file=paste0(data_path,"development_ICD_lookback_10_part_",cluster_id,".Rdata")))
                       rm(cdata)
                       
                       icd <- rbind(icd_val,icd_dev) %>%
                         select(-rownr) %>%
                         rename(indexdate=date)
                       
                       rm(icd_val)
                       rm(icd_dev)
               
                       
                       icd[["AB"]] <- icd[["A"]] | icd[["B"]]
                       icd[["CD"]] <- icd[["C"]] | icd[["D"]]
                       icd[["ST"]] <- icd[["S"]] | icd[["T"]]
                       icd[["OPQUY"]] <- icd[["O"]] | icd[["P"]] | icd[["Q"]] | icd[["U"]] | icd[["Y"]]
                       
                       
                       if(!is.null(all_icd_n_chapt)){
                         # Compute n chapters
                         icd$icd_n_chapters <- 0
                         icd_chapters <- c("AB","CD","E","F","G","H","I","J","K","L","M", "N",  "ST")
                         
                         for(k in 1:length(icd_chapters)){
                           icd$icd_n_chapters <-   icd$icd_n_chapters + as.numeric(icd[[icd_chapters[k]]])
                         }
                         
                         data <- data %>%
                           left_join(icd %>% select(LopNr,indexdate,icd_n_chapters),by=c("LopNr","indexdate"))
                         data$icd_n_chapters[is.na(data$icd_n_chapters)] <- 0
                         
                         print("Filtering ICD n chapters")
                         
                         data <- data %>%
                           filter(icd_n_chapters %in% all_icd_n_chapt)
                         
                       }
                    
                       if(!is.null(all_icd_chapt)){
                         
                         data <- data %>%
                           left_join(icd %>% select(LopNr,indexdate,any_of(all_icd_chapt)),by=c("LopNr","indexdate"))
                         
                         print("Filtering ICD chapters")
                         
                         included <- FALSE
                         for(k in 1:length(all_icd_chapt)){
                           data[[all_icd_chapt[k]]][is.na(data[[all_icd_chapt[k]]])] <- FALSE
                           included <- included | data[[all_icd_chapt[k]]]
                         }
                         
                         print(paste0(sum(included)," included"))
                         data <- data[included,]
                         rm(included)
                       }
                       
                       rm(icd)
                     }
                     
                     
                     print(paste0("Nrow data=",nrow(data)))
                     print(colnames(data))
                     ret <- list()
                     
                     for(k in 1:length(risk_score_name)){
                       cat("\n")
                       print(paste0("Risk score ",risk_score_name[k]))
                      
                       ###
                       # Cindex
                       print( Sys.time() - start_time)
                       print("Computing C-index")
                       data$comorbidity_index <- data[[paste0("comorbidity_index_",risk_score_name[k])]]
                       
                       if(neg_to_zero[k]){
                         data$comorbidity_index[data$comorbidity_index<0]<-0
                       }
                       
                       ret_inner <- compute_cindex_inner(data=data,
                                                         subgroups=subgroups,
                                                         subgroup_function=subgroup_function,
                                                         boot_resample=FALSE,
                                                         at_times=at_times)
                       
                       ret_inner$risk_score_name <- risk_score_name[k]
                       
                       ret[[k]] <- ret_inner
                       
                     }
                     
                     rm(data)
                     gc()
                     
                     
                     names(ret) <- risk_score_name
                     
                     ret <- do.call(rbind,ret)
                     
                     ret$cluster_id <- cluster_id
                     
                     if(save_parts){
                       save(ret,file=paste0(save_path,"cindex_inner_",save_name,"_part",cluster_id,".Rdata"))
                     }
                     
                     return(ret)
                     
                   })
  print( Sys.time() - start_time)
  stopCluster(cl)
  gc(verbose = FALSE)
  
  point_ests <- do.call(rbind,ret)
  
  point_ests <- point_ests %>%
    group_by(risk_score_name,timepoint,subgroup_index) %>%
    mutate(nrow=n(),
           Ntot=sum(N),
           N_mean=mean(N),
           N_sd=sd(N),
           N_unique=sum(Nunique),
           N_unique_mean=mean(Nunique),
           N_unique_sd=sd(Nunique),
           N_na=sum(is.na(cindex)),
           Which_na=paste0(cluster_id[is.na(cindex)],collapse="|"),
           p = N/Ntot,
           p_sd = sd(p),
           cindex_mean=sum(p*cindex),
           cindex_sd=sd(cindex),
           cindex_varpooled=sum(cindex_var*p^2)) %>%
    ungroup()
  
  point_ests_aggr <- point_ests %>%
    group_by(risk_score_name,timepoint,subgroup_index) %>%
    slice(1) %>%
    select(risk_score_name,
           timepoint,
           subgroup_index,
           N_mean,
           N_sd,
           N_unique,
           N_unique_mean,
           N_unique_sd,
           cindex_mean,
           cindex_sd,
           cindex_varpooled,
           N_na,
           Which_na) %>%
    ungroup() %>%
    arrange(risk_score_name,subgroup_index, -timepoint)
  
  ret <- list("point_ests"=point_ests,
              "point_ests_aggr"=point_ests_aggr)
  
  cat("\n")
  print("Done!")
  print( Sys.time() - start_time)
  return(ret)
}




###
# End
###