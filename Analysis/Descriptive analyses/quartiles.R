###
# quartiles of MDCI and DCI in total and by age group

X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

ncores <- 16
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

ret_list <- list()

nbatch <- 64/ncores

print(nbatch)

cluster_id <- 1

for(j in 1:nbatch){
  print(j)
  ###
  # Load data first and then iterate over subgroups
  ret <- parLapply(cl = cl,
                   cluster_indices[cluster_indices %% nbatch == (j-1)],
                   fun = function(cluster_id){
                     require(dplyr)
                     setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\"))
                     
                     cluster_id <<- cluster_id # make seed global
                
                     print("Loading data...")
                     Sys.sleep(cluster_id/16)
             
                     
                     loads <- load(file=paste0(data_path,"cindex_cache_part",cluster_id,".Rdata"))
                     
                     data <- data %>%
                       filter(indexyear >= 2006 & indexyear <= 2022)
                     
                     data <- data %>%
                       select(LopNr,rowid,sex,indexage,comorbidity_index_MDCI,comorbidity_index_DCI) 
                    
                     
                     data$age_group <- cut(data$indexage,breaks=c(18,40,50,60,70,80,90,100),include.lowest = TRUE,right=FALSE)
             
                     
                     findq <- function(data,v){
                       quantile(data[[v]],probs=c(0,0.25,0.5,0.75,1))
                     }
                     
                     quantiles <- list("MDCI_all"=findq(data,v="comorbidity_index_MDCI"),
                                       "MDCI_men"=findq(data %>% filter(sex %in% "man"),v="comorbidity_index_MDCI"),
                                       "MDCI_women"=findq(data %>% filter(sex %in% "woman"),v="comorbidity_index_MDCI"),
                                       
                                       "DCI_all"=findq(data,v="comorbidity_index_DCI"),
                                       "DCI_men"=findq(data %>% filter(sex %in% "man"),v="comorbidity_index_DCI"),
                                       "DCI_women"=findq(data %>% filter(sex %in% "woman"),v="comorbidity_index_DCI"))
                     
                     ageg <- sort(unique(data$age_group))
                     for(j in 1:length(ageg)){
                       print(j)
                       
                       datat <- data %>%
                         filter(age_group %in% ageg[j])
                       
                       temp <- list("MDCI_all"=findq(datat,v="comorbidity_index_MDCI"),
                                    "MDCI_men"=findq(datat %>% filter(sex %in% "man"),v="comorbidity_index_MDCI"),
                                    "MDCI_women"=findq(datat %>% filter(sex %in% "woman"),v="comorbidity_index_MDCI"),
                                    
                                    "DCI_all"=findq(datat,v="comorbidity_index_DCI"),
                                    "DCI_men"=findq(datat %>% filter(sex %in% "man"),v="comorbidity_index_DCI"),
                                    "DCI_women"=findq(datat %>% filter(sex %in% "woman"),v="comorbidity_index_DCI"))
                       
                       names(temp)<-paste0(names(temp),"_",ageg[j])
                       
                       quantiles <- c(quantiles,
                                      temp)
                       
                     }
                     
                     quantiles_n <- names(quantiles)
                     
                     quantiles <- do.call(rbind,quantiles)
                     
                     colnames(quantiles)<-paste0("Q",c(0,25,50,75,100))
                     quantiles <- as.data.frame(quantiles)
                     quantiles$subgroup <- quantiles_n
                     
                     quantiles$sex <- c("all","men","women")
                     quantiles$age <- rep(c("all",as.character(ageg)),each=6)
                     
                     quantiles$comorbidity_index <- c("MDCI","MDCI","MDCI","DCI","DCI","DCI")
                     
                     return(quantiles)
                     
                   })
  ret_list[[j]]<-ret
  rm(ret)
}

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)

ret_list <- unlist(ret_list,recursive=FALSE)

for(j in 1:64){
  ret_list[[j]]$cluster_id <- j
}
data <- do.call(rbind,ret_list)
dim(data)
rm(ret_list)
gc()


data <- data %>%
  group_by(subgroup) %>%
  mutate(Q0=mean(Q0),
         Q25=mean(Q25),
         Q50=mean(Q50),
         Q75=mean(Q75),
         Q100=mean(Q100)) %>%
  slice(1) %>%
  select(-cluster_id)

data <- data %>%
  arrange(comorbidity_index,sex,age)

view(data %>% filter(sex %in% "all"))


quantiles <- data

save(quantiles,file=paste0(save_path,"comorbidity_quantiles.Rdata"))

##
# End
##

###
# End
###





###
# End
###