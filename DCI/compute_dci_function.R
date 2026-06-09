
compute_dci <- function(path_weights_men,
                        path_weights_women,
                        path_lm_reg,
                        to_path,
                        cohort_path,
                        cohort_path_extra=NULL,
                        cohort_colnames,
                        cohort_file_name,
                        cohort_file_name_extra=NULL,
                        cohort_file_tail=NULL,
                        extraname="",
                        ncores=64,
                        cluster_ids=NULL,
                        collapse_rowid=FALSE){
  
  require(dplyr)
  require(parallel)
  
  cohort_newnames <- c("LopNr","date","sex") # forced to be this
  
  start_time <- Sys.time()
  
  cl <- makeCluster(ncores)
  if(is.null(cluster_ids)){
    cluster_ids <- 1:ncores
  }

  print(cluster_ids)
  clusterExport(cl, c("path_weights_men",
                      "path_weights_women",
                      "path_lm_reg",
                      "to_path",
                      "cohort_path",
                      "cohort_colnames",
                      "cohort_path_extra",
                      "cohort_file_name_extra",
                      "cohort_newnames",
                      "cohort_file_name",
                      "cohort_file_tail",
                      "extraname",
                      "extraname",
                      "collapse_rowid"),envir = environment())
  
  clusterEvalQ(cl, library(dplyr))
  
  a <- parLapply(cl = cl, cluster_ids , fun = function(cluster_id){
    cluster_id <<- cluster_id # make seed global
  
    Sys.sleep( runif(1,0.01,1) )
    
    # load weights
    loads <- load(path_weights_men) 
    loads <- load(path_weights_women) 
    
    
    # if cohort_file_name_extra # load Kon
    if(!is.null(cohort_path_extra)){
      cohort_extra <- get(load(paste0(cohort_path_extra,cohort_file_name_extra,cluster_id,".Rdata")))
      cohort_extra <- cohort_extra %>%
        select(LopNr,Kon)
      rm("D")
    }
    
    # Load cohort
    if(!is.null(cohort_file_tail)){
      cohort <- get(load(paste0(cohort_path,cohort_file_name,cluster_id,cohort_file_tail,".Rdata")))
    } else{
      cohort <- get(load(paste0(cohort_path,cohort_file_name,cluster_id,".Rdata")))
    }
    
    cohort <- cohort %>%
      select(all_of(cohort_colnames))
    
    
    if(!is.null(cohort_path_extra)){
      cohort <- cohort %>%
        left_join(cohort_extra,by="LopNr")
    }
    
    colnames(cohort) <- cohort_newnames
    
    cohort <- cohort %>%
      arrange(LopNr,date) %>%
      group_by(LopNr) %>%
      unique() %>%
      mutate(rownr=1:n())
    
    cohort <- cohort[!is.na(cohort$date),]
   
    
    n_rownr <- sort(unique(cohort$rownr))
    
    # load pdr
    drug_chunk <- get(load(paste0(path_lm_reg,"drug_part",cluster_id,".Rdata")))
    drug_chunk <- drug_chunk %>%
      select(LopNr,EDATUM,ATC) %>%
      rename(edatum=EDATUM,atc=ATC)
    
    ###
    # Helper function
    dci_inner <- function(lmreg,weights){
      
      ret <- merge(lmreg,
                   weights,
                   by.x="atc",
                   by.y="ATC5_365",
                   all.x=TRUE) %>%
        select(LopNr,atc,logHR,date,sex) %>% 
        filter(!is.na(logHR)) %>%
        group_by(LopNr,date,atc) %>%
        unique() %>%
        ungroup() 
      
      ret <- ret %>%
        group_by(LopNr,date) %>%
        unique() %>%
        mutate(dci=sum(logHR)) %>%
        ungroup() %>%
        select(LopNr,date,dci) %>%
        unique() %>%
        arrange(LopNr,date)
      
      return(ret)
    }
    
    dci_data_rowid <- list()
    
    for(j in n_rownr){
      cohort_temp <- cohort %>% 
        filter(!is.na(date)) %>%
        filter(rownr==j)
      
      ###
      # Compute for each row
      temp <- cohort_temp %>% 
        left_join(drug_chunk,by="LopNr") %>%
        filter(!is.na(edatum)) %>%
        filter(date <= 365.24 + edatum  & edatum < date) %>%
        mutate(atc=substr(atc,start=1,stop=5))
      
      print("Computing...")
      ###
      # Compute DCI for men
      dci_data <- dci_inner(lmreg=temp %>% 
                              filter(sex %in% "man"),
                            weights=Men)
      
      # Compute DCI for women
      dci_data_women <- dci_inner(lmreg=temp %>% 
                                    filter(sex %in% "woman"),
                                  weights=Women)
      rm(temp)
      
      dci_data <- rbind(dci_data,dci_data_women)
      rm(dci_data_women)
      gc()
      
      dci_data <- cohort_temp %>%
        left_join(dci_data,
                  by=c("LopNr","date"))
      
      dci_data$dci[is.na(dci_data$dci)]<-0
      
      
  
      if(!collapse_rowid){
        if(!extraname %in% ""){
          filename <- paste0(extraname,"_DCI_part_",cluster_id)
        } else{
          filename <- paste0("DCI_part_",cluster_id)
        }
        
        if(!is.null(cohort_file_tail)){
          filename <- paste0(filename,cohort_file_tail) 
        }
        
        save(dci_data,
             file=paste0(to_path,filename,"_rowid_",j,".Rdata"))
      } else{
        dci_data_rowid[[j]]<-dci_data
      }
    
      
    }

    if(collapse_rowid){
    
      dci_data <- do.call(rbind,dci_data_rowid)
     
      
      if(!extraname %in% ""){
        filename <- paste0(extraname,"_DCI_part_",cluster_id)
      } else{
        filename <- paste0("DCI_part_",cluster_id)
      }
      
      if(!is.null(cohort_file_tail)){
        filename <- paste0(filename,cohort_file_tail) 
      }
      
      save(dci_data,
           file=paste0(to_path,filename,".Rdata"))
    } 
    
    
    # extract data including only given atc codes
    rm("drug_chunk")
    rm("cohort")
    rm("cohort_temp")
    rm("temp")
    gc()

    
    NULL
  })
  
  
  
  print( Sys.time() - start_time)
  stopCluster(cl)
  gc()
}
