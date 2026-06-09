
X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
require(tidyverse)

log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"
background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")
figs_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")

cluster_indices <- 1:64

ncores <- 32

###
# Load data...
require(parallel)

start_time <- Sys.time()
print(cluster_indices)
ncores <- min(ncores,length(cluster_indices))


print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("X","data_path","save_path"),envir = environment())

ret <- parLapplyLB(cl = cl,
                   cluster_indices,
                   fun = function(cluster_id){
                     require(dplyr)
                     setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\"))
                     sink(paste0("C:\\Marcus\\Misc\\","le_data_",cluster_id,".txt"))
                     
                     cluster_id <<- cluster_id # make seed global
                     set.seed(cluster_id)
                     
                     print("Loading data...")
                     Sys.sleep(cluster_id/16)
                     
                     loads <- load(file=paste0(data_path,"general_part",cluster_id,".Rdata"))
                     
                     D <- D %>%
                       select(LopNr,birthdate,Kon,Fodelseland_EU28)
                     
                     loads <- load(file=paste0(data_path,"Validation_final_part",cluster_id,".Rdata"))
                     
                     Val <- Val %>%
                       filter(indexyear %in% 2006:2022) 
                     
                     loads <- load(file=paste0(data_path,"Development_final_part",cluster_id,".Rdata"))
                     
                     Dev <- Dev %>%
                       filter(indexyear %in% 2006:2022) 
                     
                     
                     data <- rbind(Dev,Val) %>%
                       left_join(D,by="LopNr")
                     
                     rm(Dev)
                     rm(Val)
                     rm(D)
                     
                     data <- data %>%
                       mutate(indexage= as.numeric((indexdate - birthdate))/365.24) %>%
                       filter(indexage >= 18)
               
                     data <- data %>%
                       select(-c(subset,subset2)) %>% 
                       rename(sex=Kon) %>%
                       mutate(sex=factor(sex,levels=c("man","woman"),labels=c("man","woman"))) %>%
                       rename(SwedOrigin=Fodelseland_EU28) %>%
                       mutate(SwedOrigin=SwedOrigin %in% "Sverige")
                     
                     ###
                     # Load DCI and MDCI and join
                     print("Loading MDCI")
                     
                     loads <- load(paste0(data_path,"Development_MDCI","_part_",cluster_id,".Rdata"))
                     CD <- get(loads)
                     rm(list=loads)
                     
                     loads <- load(paste0(data_path,"Validation_MDCI","_part_",cluster_id,".Rdata"))
                     CD2 <- get(loads)
                     rm(list=loads)
                     
                     mdci <- rbind(CD,CD2) %>%
                       select(LopNr,date,MDCI) %>%
                       rename(indexdate=date)
                     
                     data <- data %>%
                       left_join(mdci,
                                 by=c("LopNr","indexdate"))
                     rm(mdci)
                     gc()
                     
                     ### 
                     # DCI
                     print("Loading DCI")
                     
                     loads <- load(paste0(data_path,"Development_DCI","_part_",cluster_id,".Rdata"))
                     CD <- get(loads)
                     rm(list=loads)
                     
                     loads <- load(paste0(data_path,"Validation_DCI","_part_",cluster_id,".Rdata"))
                     CD2 <- get(loads)
                     rm(list=loads)
                     
                     dci <- rbind(CD,CD2) %>%
                       select(LopNr,date,dci) %>%
                       rename(indexdate=date)
                     
                     data <- data %>%
                       left_join(dci,
                                 by=c("LopNr","indexdate")) 
                     rm(dci)
                     gc()
                     
                     ### 
                     # DCI
                     print("Loading CCI")
                     
                     loads <- load(paste0(data_path,"Development_CCI_lookback_10","_part_",cluster_id,".Rdata"))
                     CD <- get(loads)
                     rm(list=loads)
                     
                     loads <- load(paste0(data_path,"Validation_CCI_lookback_10","_part_",cluster_id,".Rdata"))
                     CD2 <- get(loads)
                     rm(list=loads)
                     
                     cci <- rbind(CD,CD2) %>%
                       select(LopNr,date,CCIw) %>%
                       rename(indexdate=date,
                              CCI=CCIw)
                     
                     data <- data %>%
                       left_join(cci,
                                 by=c("LopNr","indexdate")) 
                     rm(cci)
                     gc()
                     
                     
                     data$caltime <- as.numeric(data$indexdate-as.Date("2006-01-01"))/365.24
                     data$censor <- as.numeric(!data$censor %in% "0")
                     
                     maxfu <- 15
                     data$censor[data$timefu>maxfu] <- 0
                     data$timefu[data$timefu>maxfu] <- maxfu
                     
                     data$age <- floor(data$indexage)
                     
                     ret <- c()
                     for(a in 18:100){
                       print(a)
                       print(Sys.time())
                       data_temp <- data %>% filter(age %in% a)
                       ret[a]<-nrow(data_temp)
                       save(data_temp,file=paste0("C:\\Marcus\\ComorbidityBase_cache\\cox_life_expectancy_",a,"_",cluster_id,".Rdata"))
                     }
                 
                     return(ret)
                     
                   })

ret <- do.call(rbind,ret)


print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)
