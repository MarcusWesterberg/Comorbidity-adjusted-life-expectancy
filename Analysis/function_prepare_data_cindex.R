prepare_data_cindex <- function(ncores,
                                cluster_indices,
                                log_path,
                                data_path,
                                data_name,
                                comorb_data_path,
                                comorb_data_name,
                                comorb_colnames,
                                comorb_newcolnames,
                                background_data_path=NULL,
                                risk_score_name,
                                save_name){
  
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
                      "data_name",
                      "save_name",
                      "comorb_data_path",
                      "comorb_data_name",
                      "comorb_colnames",
                      "comorb_newcolnames",
                      "background_data_path",
                      "risk_score_name"),envir = environment())
  
  ###
  # Load data first and then iterate over subgroups
  ret <- parLapply(cl = cl,
                   cluster_indices,
                   fun = function(cluster_id){
                     
                     start_time <- Sys.time()
                     
                     require(tidyverse)
                     
                     
                   })
  
  print("Loading data...")
  nrow_data <- parLapply(cl,
                         cluster_indices,
                         fun=function(cluster_id){
                           sink(file=paste0(log_path,"cindex_log_",save_name,"_",cluster_id,".txt"))
                           start_time <- Sys.time()
                           data_list <- list()
                           
                           for(i in 1:length(data_name)){
                             print(paste0("i=",i))
                             ###
                             # load data
                             print( Sys.time() - start_time)
                             print("Loading data")
                             
                             loads <- load(paste0(data_path,data_name[i],"_part",cluster_id,".Rdata"))
                             data <- get(loads)
                             rm(list=loads)
                             gc()
                             
                             data$censor <- as.numeric(!data$censor %in% "0")
                             
                             ###
                             # load background data
                             if(!is.null(background_data_path)){
                               print("Loading background data")
                               
                               # load general
                               loads <- load(paste0(background_data_path,"general","_part",cluster_id,".Rdata"))
                               general <- get(loads)
                               rm(list=loads)
                               gc()
                               general <- general %>%
                                 select(LopNr,birthdate,birthyear,Kon,Fodelseland_EU28,UtlSvBakg )
                               
                               print("Joining data")
                               data <- data %>%
                                 left_join(general,by="LopNr")
                               
                               data <- data %>%
                                 mutate(indexage = (as.numeric(indexdate)-as.numeric(birthdate))/365.24) %>%
                                 rename(sex=Kon)
                               
                               rm("general")
                               gc()
                               # load others if needed?
                             }
                             
                             for(j in 1:length(comorb_data_name[[i]])){
                               print("Loading comorbidity")
                               print(j)
                               loads <- load(paste0(comorb_data_path,comorb_data_name[[i]][j],"_part_",cluster_id,".Rdata"))
                               CD <- get(loads)
                               if("lopnr" %in% colnames(CD)){
                                 CD <- CD %>%
                                   rename(LopNr=lopnr)
                               }
                               rm(list=loads)
                               gc()
                               
                               CD <- CD %>%
                                 select( all_of(comorb_colnames[[j]]))
                               
                               comorb_newcolnames_temp <- comorb_newcolnames
                               comorb_newcolnames_temp[3] <- paste0(comorb_newcolnames[3],"_",risk_score_name[j])
                               colnames(CD) <- comorb_newcolnames_temp
                               
                               print("Joining data")
                               data <- data %>%
                                 left_join(CD,
                                           by=c("LopNr","indexdate"))
                               rm(CD)
                               gc()
                               
                               # If no hosp visits or drugs then it is NA but should be zero
                               if(any(is.na(data[[ comorb_newcolnames_temp[3] ]]))){
                                 na <- is.na(data[[ comorb_newcolnames_temp[3] ]])
                                 data[[ comorb_newcolnames_temp[3] ]][na]<-0
                                 print(paste0("NA in ",comorb_newcolnames_temp[3],", N na=",sum(na)))
                                 #stop("Error: NA in comorbidity index not supported!")
                               }
                             }
                             
                             data_list[[i]] <- data
                             rm("data")
                             gc()
                           }
                           
                           data <- bind_rows(data_list)
                           rm(data_list) 
                           gc()
                           
                           save(data,file=paste0(data_path,"cindex_cache_part",cluster_id,".Rdata"))
                           
                           ret <- nrow(data)
                           rm("data")
                           gc()
                           sink()
                           return(ret)
                         })
  print(paste0("N total"))
  print(summary(unlist(nrow_data)))
  stopCluster(cl)
  ###
  print( Sys.time() - start_time)
}