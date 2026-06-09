###
# quartiles of MDCI and DCI in total and by age group

X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

ncores <- 32
data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\"))
save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")

###
# Load data...
require(parallel)

start_time <- Sys.time()
cluster_indices <- 1:64
print(cluster_indices)
ncores <- min(ncores,length(cluster_indices))


print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("X","data_path"),envir = environment())

print(nbatch)

cluster_id <- 1
ret <- parLapply(cl = cl,
                 cluster_indices,
                 fun = function(cluster_id){
                   require(dplyr)
                   setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\"))
                   
                   cluster_id <<- cluster_id # make seed global
                   
                   print("Loading data...")
                   Sys.sleep(cluster_id/16)
                   
                   
                   loads <- load(file=paste0(data_path,"cindex_cache_part",cluster_id,".Rdata"))
                   
                   data <- data %>%
                     filter(indexyear >= 2006 & indexyear <= 2022) %>%
                     filter(indexage >= 70 & indexage <= 79)
                   
                   data <- data %>%
                     select(LopNr,rowid,sex,indexage,comorbidity_index_MDCI,comorbidity_index_DCI) 
                   
                   
                   findq <- function(data,v){
                     quantile(data[[v]],probs=c(0,0.2,0.4,0.6,0.8,1))
                   }
                   
                   quantiles <- list("MDCI_all"=findq(data,v="comorbidity_index_MDCI"),
                                     "MDCI_men"=findq(data %>% filter(sex %in% "man"),v="comorbidity_index_MDCI"),
                                     "MDCI_women"=findq(data %>% filter(sex %in% "woman"),v="comorbidity_index_MDCI"),
                                     
                                     "DCI_all"=findq(data,v="comorbidity_index_DCI"),
                                     "DCI_men"=findq(data %>% filter(sex %in% "man"),v="comorbidity_index_DCI"),
                                     "DCI_women"=findq(data %>% filter(sex %in% "woman"),v="comorbidity_index_DCI"))
                   return(quantiles)
                   
                 })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)

ret_list <- ret

for(j in 1:64){
  ret_list[[j]] <- as.data.frame(do.call(rbind,ret_list[[j]]))
  colnames(ret_list[[j]]) <- c("Q0","Q2","Q4","Q6","Q8","Q10")
  ret_list[[j]]$type <- rownames( ret_list[[j]] )
  ret_list[[j]]$cluster_id <- j
}
ret_list <- do.call(rbind,ret_list)

data <- ret_list %>%
  group_by(type) %>%
  mutate(Q0=mean(Q0),
         Q2=mean(Q2),
         Q4=mean(Q4),
         Q6=mean(Q6),
         Q8=mean(Q8),
         Q10=mean(Q10)) %>%
  slice(1) %>%
  select(-cluster_id)

exp(data[1,-7]) # dci
exp(data[4,-7]) # mdci
save(data,file=paste0(save_path,"comorbidity_quintiles.Rdata"))

##
# End
##

###
# End
###





###
# End
###