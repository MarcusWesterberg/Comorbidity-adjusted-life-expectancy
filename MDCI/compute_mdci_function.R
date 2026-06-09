
compute_mdci <- function(from_path,
                         weights_path,
                         to_path,
                         log_path="C:\\Marcus\\mdci_logs",
                         ncores=64,
                         cohort_path,
                         cohort_colnames,
                         cohort_file_name,
                         extraname="",
                         mdci_inner,
                         process_data,
                         max_years_back=10,
                         cluster_ids=NULL,
                         drop_icd=NULL,
                         dont_drop_icd=NULL,
                         collapse_rowid=FALSE,
                         rand_rowid=FALSE,
                         sources=c("ov","sv"),
                         drop_utdatum=FALSE,
                         min.per.inner=10000,
                         alt_instr_list=NULL){
  
  print("Starting...")
  stopifnot(length(cohort_colnames)==2)
  
  cohort_newnames <- c("lopnr","date") # must be this for the rest to work
  require(dplyr)
  require(parallel)
  
  start_time <- Sys.time()
  
  ###
  # Create cluster for parallel computations
  if(is.null(cluster_ids)){
    cluster_ids <- 1:64
  }
  print(cluster_ids)
  ncores <- min(ncores,length(cluster_ids))
  
  print(paste0("Ncores=",ncores))
  
  cl <- makeCluster(ncores)
  
  clusterExport(cl, c("from_path","cohort_path","cluster_ids","to_path","cohort_colnames","cohort_file_name",
                      "cohort_newnames","mdci_inner","process_data","weights_path","log_path","max_years_back",
                      "drop_icd","dont_drop_icd","rand_rowid","drop_utdatum","min.per.inner","alt_instr_list",
                      "sources"),envir = environment())
  
  if(!is.null(alt_instr_list)){
    cluster_ids <- alt_instr_list
  }
  
  cluster_id <- cluster_ids[[1]]
  a <- parLapply(cl = cl, cluster_ids , fun = function(cluster_id){
    start_time <- Sys.time()
    
    if(!is.null(alt_instr_list)){
      rowid_select <- cluster_id$rowids
      cluster_id <- cluster_id$cluster_id[1]
   
    } else{
      rowid_select <- NULL
    }
    
    require(tidyverse)
    cluster_id <<- cluster_id # make seed global
    
    set.seed(cluster_id)
    
    cluster_id_log <- cluster_id
    if(!is.null(rowid_select)){
      cluster_id_log <- paste0(cluster_id_log,"_",paste0(rowid_select,collapse=""))
    }
    
    sink(file=paste0(log_path,"mdci_log_",cluster_id_log,".txt"))
    print(Sys.time())
    
    if(!is.null(rowid_select)){print(rowid_select)}
    
    print("Loading...")
    # load weights
    weights_data <- get(load(weights_path)) # gives weights_data
    weights_data <- weights_data$summary
    rm("f10")
    
    
    # load PAR
    loads <- load(paste0(from_path,"par_part",cluster_id,".Rdata")) # par_part
    par_part <- par_part %>%
      rename(lopnr=LopNr)
    
    
    # read cohort
    loads <- load(paste0(cohort_path,cohort_file_name,cluster_id,".Rdata"))
    cohort <- get(loads)
    rm(list=loads)
    
    cohort <- cohort[,cohort_colnames]
    
    colnames(cohort) <- cohort_newnames
    
    cohort <- cohort[!is.na(cohort$date),]
    
    
    
    ### ### ### ### ### ### ### ### ### ### 
    ### Done reading 
    ### ### ### ### ### ### ### ### ### ### 
    
    print( Sys.time() - start_time)
    
    print("Merging...")
    
    cohort <- cohort %>%
      filter(lopnr %in% par_part$lopnr)
    
    par_part <- par_part %>%
      filter(lopnr %in% cohort$lopnr) %>%
      select(lopnr,indatum,utdatum,hdia,diagnos,source)  %>% 
      
      filter(!is.na(indatum)) %>%
      mutate(indatum=as.Date(indatum,origin="1970-01-01"))  %>%
      mutate(utdatum=as.Date(utdatum,origin="1970-01-01"))  
    
    cohort <- cohort %>%
      group_by(lopnr) %>%
      mutate(rowid=1:n()) %>%
      ungroup()
    
    nrowids <- sort(unique(cohort$rowid))
    maxrowid <- max(nrowids)
    
    if(rand_rowid){
      print("Random rowid...")
      ### 
      # Added uniform distribution over rowids 
      cohort <- cohort %>%
        group_by(lopnr) %>%
        mutate(rowid=sample(1:maxrowid,size=n(),replace = FALSE)) %>%
        ungroup()
    }
    
    
    # x <- cohort %>%
    #   group_by(lopnr) %>%
    #   mutate(m=length(unique(rowid))) %>%
    #   ungroup()
    ###
    
    if(!is.null(rowid_select)){
      nrowids <- rowid_select
    }
    
    
    datalist_rowid <- list()
    
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
    for(j in nrowids){
      print(paste0("rowid ",j))
      
      print( Sys.time() - start_time)
      
      ###
      # merge 
      PAR <- cohort %>%
        dplyr::rename_all(tolower) %>%
        
        filter(rowid %in% j) %>%
        
        left_join(par_part,by="lopnr") #%>% 
        
        #mutate(indatum=as.Date(indatum,origin="1970-01-01"))  
      
      
      ###   ###   ###   ###   ###   ###  ###   ###   ###   ###
      # look only 10 years back...
      PAR <- PAR %>% 
        filter(indatum < date & indatum >= (date-round(max_years_back*365.25) )) 
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
      
      
      if(drop_utdatum){
        print("Droppin utdatum...")
        PAR$utdatum2 <- PAR$utdatum
        PAR$utdatum2[is.na(PAR$utdatum2) & PAR$source %in% "sluten"] <- PAR$indatum[is.na(PAR$utdatum2) & PAR$source %in% "sluten"] + 30
        
        PAR$utdatum2[is.na(PAR$utdatum2) & PAR$source %in% "oppen"] <- PAR$indatum[is.na(PAR$utdatum2) & PAR$source %in% "oppen"]
        
        PAR <- PAR %>%
          mutate(utdatum2=as.Date(utdatum2,origin="1970-01-01"))  %>% 
          filter( utdatum2 < date  ) %>%
          select(- utdatum2)
      }
    
      
      #rm("par_part")
      
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
      print("Processing..")
      ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
      oppen <- PAR %>% 
        filter(source %in% "oppen") 
      
      sluten <- PAR %>% 
        filter(source %in% "sluten") 
      
      rm("PAR")
      gc()
      
      
      if("ov" %in% sources){
        oppen <- oppen %>%
          dplyr::rename_all(tolower) %>%
          select(lopnr,indatum,hdia,diagnos,date) %>%
          rowwise() %>%
          mutate(diagnos=paste0(c(hdia," ",diagnos),collapse="") ) %>%
          ungroup() %>%
          mutate(source="oppen") %>%
          #mutate(indatum=as.Date(indatum,origin="1970-01-01"))  %>% 
          mutate(duration=0.5) 
      }

      if("sv" %in% sources){
        sluten <- sluten %>%
          dplyr::rename_all(tolower) %>%
          select(lopnr,indatum,utdatum,hdia,diagnos,date) %>%
          rowwise() %>%
          mutate(diagnos=paste0(c(hdia," ",diagnos),collapse="") ) %>%
          ungroup() %>%
          mutate(source="sluten") %>%
          #mutate(indatum=as.Date(indatum,origin="1970-01-01"),
          #       utdatum=as.Date(utdatum,origin="1970-01-01")) %>%
          mutate(duration=as.numeric(utdatum-indatum)) %>% 
          select(-utdatum) %>%
          mutate(duration=ifelse(duration==0,0.5,duration)) %>%
          mutate(duration=ifelse(duration<0,0.5,duration)) 
      }
      
      if("sv" %in% sources & "ov" %in% sources){
        PAR <- rbind(oppen,sluten) 
      } else if("sv" %in% sources){
        PAR <- sluten
      } else if("ov" %in% sources){
        PAR <- oppen
      } else{
        stop("None of sv or ov are in sources!")
      }
      
      rm("oppen")
      rm("sluten")
      gc()
      
      PAR <- PAR %>% 
        mutate(diagnos = strsplit(as.character(diagnos), " ")) %>%
        tidyr::unnest(diagnos) %>%
        select(-c(source) ) %>%
        unique() %>%
        mutate(is_hdia=diagnos==hdia) %>%
        select(-hdia)
      
      # remove odd characters etc
      PAR <- PAR %>%
        filter(!diagnos %in% c(""," ","\t","\t\t","_atc","1","10","6500","AVVAKTA"))
      
      # remove special characters from codes
      PAR$diagnos <- stringr::str_replace_all(PAR$diagnos, "[^[:alnum:]]", "") 
      
      # remove digits, odd characters etc and prostate cancer diagnosis
      PAR <- PAR %>%
        filter(!diagnos %in% c("","\t","\t\t","_atc","1","10","6500","AVVAKTA")) %>%
        filter(!grepl("^[[:digit:]]+", diagnos) ) # %>%
      # filter(!diagnos %in% c("C61","C619"))
      
      if(!is.null(drop_icd)){
        print("Dropping ICD")
        print(nrow(PAR))
        print(drop_icd)
        print(dont_drop_icd)
        
        remove_these_icd <- rep(FALSE,nrow(PAR))
        dont_drop_these_icd <- rep(FALSE,nrow(PAR))
        
        if(!is.null(dont_drop_icd)){
          dont_drop_these_icd <- grepl(PAR$diagnos,pattern=paste0(dont_drop_icd,collapse="|"))
        } 
        
        remove_these_icd <- remove_these_icd | (  grepl(PAR$diagnos,pattern=paste0(drop_icd,collapse="|"))  & !dont_drop_these_icd  )
        
        PAR <- PAR[!remove_these_icd,]
        print(nrow(PAR))
        rm(remove_these_icd)
      }

      
      
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
      # MDCI
      print("Computing MDCI...")
      
      par_lopnr <- unique(PAR$lopnr)
      N <- ceiling(min(min.per.inner,length(par_lopnr)/min.per.inner))
      
      print(length(par_lopnr)/N)
      
      if(length(par_lopnr) > 0 & length(par_lopnr) < min.per.inner){
        N <- 1
      }
      
      selected <- rep(1:N,length.out=length(par_lopnr))
      
      mdci_data_list <- list()
      if(N>0){
        for(k in 1:N){
          print(k)
          print( Sys.time() - start_time)
          
          ###
          # Compute MDCI - using process_data function inside
          PAR_temp <- PAR %>% filter(lopnr %in% par_lopnr[selected==k])
          mdci_data <- mdci_inner(PAR_inner=PAR_temp,
                                  weights_data=weights_data)
          gc(verbose = FALSE)
          print( Sys.time() - start_time)
          
          ###
          # Join
          PAR_temp <- PAR_temp %>%
            select(lopnr,date) %>%
            group_by(lopnr,date) %>%
            slice(1) %>%
            ungroup()
          
          mdci_data_list[[k]] <- mdci_data %>%
            left_join(PAR_temp,by="lopnr") %>%
            mutate(rowid=j)
          rm("PAR_temp")
          rm("mdci_data")
          gc(verbose = FALSE)
        }
        rm("PAR")
        
        print( Sys.time() - start_time)
        print("Finalizing...")
        
        mdci_data <- do.call(rbind,mdci_data_list) %>% 
          arrange(lopnr)
        
        print(dim(mdci_data))
        rm("mdci_data_list")
        gc(verbose = FALSE)
        print( Sys.time() - start_time)
      } else{
        mdci_data <- NULL
      }
      print( Sys.time() - start_time)
      print(Sys.time())
      
      
      if(!collapse_rowid){
        print("Saving...")
        
        if(!extraname %in% ""){
          filename <- paste0(extraname,"_MDCI","_rowid_",j,"_part_",cluster_id)
        } else{
          filename <- paste0("MDCI","_rowid_",j,"_part_",cluster_id)
        }
        
        save(mdci_data,
             file=paste0(to_path,filename,".Rdata"))
        rm(mdci_data)
        gc(verbose = FALSE)
        print( Sys.time() - start_time)
      } else{
        datalist_rowid[[j]] <- mdci_data
      }
      
      
      
    } # End loop over j
    
    if(collapse_rowid){
      print("Collapse by rowid")
      print("Saving...")
      
      if(!extraname %in% ""){
        filename <- paste0(extraname,"_MDCI","_part_",cluster_id)
      } else{
        filename <- paste0("MDCI","_part_",cluster_id)
      }
      
      
      mdci_data <- do.call(rbind,
                           datalist_rowid)
      
      save(mdci_data,
           file=paste0(to_path,filename,".Rdata"))
      rm("mdci_data")
      gc(verbose = FALSE)
    } 
    
    rm("datalist_rowid")
    
    sink()
  })
  
  stopCluster(cl)
  
  
  
  
  ###
  # Done!
  print("Done!")
  print( Sys.time() - start_time)
  
}


