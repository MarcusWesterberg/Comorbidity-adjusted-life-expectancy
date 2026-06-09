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
                     select(LopNr,indexdate,rowid,indexyear,sex,indexage,comorbidity_index_DCI)
                   gc()
                   
                   # 1  1.59  2.94  8.64
                   data$DCI <- exp(data$comorbidity_index_DCI)
                   data$DCI <- 1*(data$DCI<1) + 2*(data$DCI>=1 & data$DCI <1.59) + 3*(data$DCI>=1.59 & data$DCI< 2.94) + 4*(data$DCI>=2.94 & data$DCI <8.64)  + 5*(data$DCI >= 8.64) 
                   
                   data <- data %>%
                     filter(indexyear>=2006 & indexyear<2023)
                   
                   data$indexage <- floor(data$indexage)
                   
                   d <- data %>%
                     group_by(indexage,sex) %>%
                     summarize(n=n(),
                               n0=sum(DCI %in% 1),
                               n1=sum(DCI %in% 2),
                               n2=sum(DCI %in% 3),
                               n3=sum(DCI %in% 4),
                               n4=sum(DCI %in% 5),.groups="keep")
                   
                   rm(data) 
                   gc()
                   
                   return(d)
                   
                 })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)


n_dci <- do.call(rbind,ret) %>%
  group_by(indexage,sex) %>%
  summarize(n=sum(n),
            n0=sum(n0),
            n1=sum(n1),
            n2=sum(n2),
            n3=sum(n3),
            n4=sum(n4),.groups="keep") %>%
  arrange(indexage,sex)

n_dci
save(n_dci,file=paste0(save_path,"n_dci.Rdata"))





