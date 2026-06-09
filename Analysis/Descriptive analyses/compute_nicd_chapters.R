# compute_
X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

ncores <- 16
log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")

setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\"))


cluster_indices<-1:64



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# All_ages 1 - focus on age


print("Starting...")

require(dplyr)
require(parallel)

start_time <- Sys.time()

print(cluster_indices)
ncores <- min(ncores,length(cluster_indices))


print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("log_path",
                    "data_path"),envir = environment())

###
# Load data first and then iterate over subgroups
ret <- parLapply(cl = cl,
                 cluster_indices,
                 fun = function(cluster_id){
                   
                   
                   cluster_id <<- cluster_id # make seed global
                   set.seed(cluster_id)
                   start_time <- Sys.time()
                   
                   require(tidyverse)
                   sink(file=paste0(log_path,"nicdchapters_log__",cluster_id,".txt"))
                   print(cluster_id)
                   print(Sys.time())
                   
                   cat("\n ... \n")
                   print(log_path)
                   cat("\n ... \n")
                   
                   
                   print("Loading data...")
                   Sys.sleep(cluster_id/4)
                   
                   loads <- load(file=paste0(data_path,"cindex_cache_part",cluster_id,".Rdata"))
                   
                   ###
                   # keep only needed columns
                   data <- data %>%
                     select(LopNr,indexdate,rowid,indexyear,sex,indexage)
                   gc()
                   
                   
                 
                   data <- data %>%
                     filter(indexyear>=2006 & indexyear<2023)
                   
                   data$indexage <- floor(data$indexage)
                   
                   icd_val <- get(load(file=paste0(data_path,"validation_ICD_lookback_10_part_",cluster_id,".Rdata")))
                   rm(cdata)
                   
                   icd_dev <- get(load(file=paste0(data_path,"development_ICD_lookback_10_part_",cluster_id,".Rdata")))
                   rm(cdata)
                   
                   icd <- rbind(icd_val,icd_dev) %>%
                     select(-rownr) %>%
                     rename(indexdate=date)
                   
                   rm(icd_val)
                   rm(icd_dev)
                   
                   data <- data %>%
                     left_join(icd,by=c("LopNr","indexdate"))
                   rm(icd)
                   
                   
                   gc()
                   ###
                   # Add ICD chapters indicators
                   icd_chapters <- c("A","B","C","D","CD_can","E","F","G","H","I","J","K","L","M", "N","R", "S","T", "Z","O","P","Q","U","Y")
                   
                   # Create indicator 2, 3, 4... or more
                   
                   for(k in 1:length(icd_chapters)){
                     data[[icd_chapters[k]]][is.na(data[[icd_chapters[k]]])] <- FALSE
                   }
                   
                   data[["AB"]] <- data[["A"]] | data[["B"]]
                   data[["CD"]] <- data[["C"]] | data[["D"]]
                   data[["ST"]] <- data[["S"]] | data[["T"]]
                   data[["OPQUY"]] <- data[["O"]] | data[["P"]] | data[["Q"]] | data[["U"]] | data[["Y"]]
                   
                   data$n_chapters <- 0
                
                   icd_chapters <- c("AB","CD","E","F","G","H","I","J","K","L","M", "N","ST")
                   
                   for(k in 1:length(icd_chapters)){
                     data$n_chapters <- data$n_chapters + as.numeric(data[[icd_chapters[k]]])
                   }
                   
                   d <- data %>%
                     group_by(indexage,sex) %>%
                     summarize(n=n(),
                               n0=sum(n_chapters %in% 0),
                               n1=sum(n_chapters %in% 1),
                               n2=sum(n_chapters %in% 2:3),
                               n3=sum(n_chapters %in% 4:5),
                               n4=sum(n_chapters >= 6),.groups="keep")
                   
                   rm(data) 
                   gc()
                   
                   return(d)
                   
                 })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)


n_icdchapters <- do.call(rbind,ret) %>%
  group_by(indexage,sex) %>%
  summarize(n=sum(n),
            n0=sum(n0),
            n1=sum(n1),
            n2=sum(n2),
            n3=sum(n3),
            n4=sum(n4),.groups="keep") %>%
  arrange(indexage,sex)
  

save(n_icdchapters,file=paste0(save_path,"n_icdchapters.Rdata"))





