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
                     select(LopNr,indexdate,rowid,indexyear,sex,indexage,comorbidity_index_MDCI)
                   gc()
                   
                   
                   d <- data  %>% 
                     filter(indexyear %in% 2014:2022) %>%
                     filter(indexage>59 & indexage < 61)
                   
                   # 0.88     1  1.11  1.64
                   data$MDCI <- exp(data$comorbidity_index_MDCI)
                   data$MDCI <- 1*(data$MDCI<0.88) + 2*(data$MDCI >= 0.88 & data$MDCI< 1) + 3*(data$MDCI >= 1 & data$MDCI<1.11) + 4*(data$MDCI>= 1.11 & data$MDCI<1.64)  +  5*(data$MDCI>=1.64) 
         
                   data <- data %>%
                     filter(indexyear>=2006 & indexyear<2023)
                   
                   data$indexage <- floor(data$indexage)
                   
                   d <- data %>%
                     group_by(indexage,sex) %>%
                     summarize(n=n(),
                               n0=sum(MDCI %in% 1),
                               n1=sum(MDCI %in% 2),
                               n2=sum(MDCI %in% 3),
                               n3=sum(MDCI %in% 4),
                               n4=sum(MDCI %in% 5),.groups="keep")
                   
                   rm(data) 
                   gc()
                   
                   return(d)
                   
                 })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)


n_mdci <- do.call(rbind,ret) %>%
  group_by(indexage,sex) %>%
  summarize(n=sum(n),
            n0=sum(n0),
            n1=sum(n1),
            n2=sum(n2),
            n3=sum(n3),
            n4=sum(n4),.groups="keep") %>%
  arrange(indexage,sex)

n_mdci
save(n_mdci,file=paste0(save_path,"n_mdci.Rdata"))





