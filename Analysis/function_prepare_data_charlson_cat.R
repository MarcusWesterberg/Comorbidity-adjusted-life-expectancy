prepare_data_charlson_cat <- function(ncores=ncores,
                                      cluster_indices,
                                      log_path,
                                      comorb_data_path,
                                      comorb_data_name,
                                      data_path,
                                      #data_name,
                                      save_name="charlson"){
  
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
                      #"data_name",
                      "save_name",
                      "comorb_data_path",
                      "comorb_data_name"),envir = environment())
  
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
                           
                           for(i in 1:length(comorb_data_name)){
                             print(paste0("i=",i))
                             ###
                             # load data
                             print( Sys.time() - start_time)
                             print("Loading data")
                       
                             if(grepl(comorb_data_name[[i]],pattern="CCI")){
                               print("Loading CCI")
                     
                               loads <- load(paste0(comorb_data_path,comorb_data_name[[i]],"_part_",cluster_id,".Rdata"))
                               CD <- get(loads)
                               rm(list=loads)
                               gc()
                               
                               CD[is.na(CD)] <- 0
                               CD <- CD %>%
                                 rename(indexdate=date)
                               
                               data_list[[i]] <- CD
                             }
                          
                             rm("CD")
                             gc()
                           }
                           
                           data <- bind_rows(data_list)
                           rm(data_list) 
                           gc()
                           
                           lb <- strsplit(comorb_data_name[[1]],split="_")[[1]]
                           lb <- lb[length(lb)]
                           save(data,file=paste0(data_path,save_name,lb,"_cache_part",cluster_id,".Rdata"))
                           
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
