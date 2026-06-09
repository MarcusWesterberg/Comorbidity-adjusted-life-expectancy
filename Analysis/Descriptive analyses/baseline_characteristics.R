###
# Baseline characteristics 
X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

ncores <- 16
data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\"))

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
                     set.seed(cluster_id)
                     source(paste0(".\\Scripts\\Table_wrap.R"))
                     print("Loading data...")
                     Sys.sleep(cluster_id/16)
                     
                     
                     # 
                     # loads <- load(file=paste0(data_path,"charlson10_cache_part",cluster_id,".Rdata"))
                     # cci_data <- data
                     
                     loads <- load(file=paste0(data_path,"cindex_cache_part",cluster_id,".Rdata"))
                    
                     # data <- data %>%
                     #   left_join(cci_data %>% 
                     #               select(-any_of(c("rownr","CCIw","CCIunw"))),
                     #             by=c("LopNr","indexdate"))
                     # rm(cci_data)
                     # gc()
                     
                     # data$comorbidity_index_CCI10[is.na(data$comorbidity_index_CCI10)]<-0
                     
                     data$indexage <- floor(data$indexage)
                     
                     data <- data %>%
                       filter(indexyear>=2006 & indexyear<2023) %>%
                       filter(indexage >= 18 & indexage <= 120)
                       
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
                     
                     ###
                     # Load healthcare encounters
                     he_val <- get(load(file=paste0(data_path,"Validation_HE_part_",cluster_id,".Rdata")))
                     rm(he_data)
                     
                     he_dev <- get(load(file=paste0(data_path,"Development_HE_part_",cluster_id,".Rdata")))
                     rm(he_data)
                     
                     he <- rbind(he_val,he_dev) %>%
                       rename(indexdate=date) %>%
                       select(LopNr,indexdate,any1,any5,any10)
                     
                     rm(he_val)
                     rm(he_dev)
                     
                     data <- data %>%
                       left_join(he,by=c("LopNr","indexdate"))
                     rm(he)
                     
                     data$any1[is.na(data$any1)]<-0
                     data$any10[is.na(data$any10)]<-0
                     data$any5[is.na(data$any5)]<-0
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
                    
                     icd_chapters <- c("AB","CD","E","F","G","H","I","J","K","L","M", "N",  "ST")
                 
                     for(k in 1:length(icd_chapters)){
                       data$n_chapters <-   data$n_chapters + as.numeric(data[[icd_chapters[k]]])
                     }
                     
                     
                     
         
                     
                     # define columns
                     C = cbind(data$indexage >= 18 & data$indexage <= 120,
                               data$indexage <50,
                               data$indexage >=50 &  data$indexage <70,
                               data$indexage >=70,
                               data$sex %in% "man",
                               data$sex %in% "woman")
                     cn <- c("All","Age 18-49 years","Age 50-69 years","Age >=70","Men","Women")
                      
                    
                     ### 
                     # Create table
                     
                     dig_percent <- 0
                     dig <- 1
                     
                     table_1 <- list()
                     data$N <- TRUE
                     table_1$N   <- table_wrap(X=data,
                                               varname="N",
                                               variable_name="N",
                                               C=C,
                                               cn=cn,
                                               continuous=FALSE,
                                               cut_it=FALSE,
                                               breaks=NULL,
                                               labels=NULL,
                                               levels=NULL,
                                               add_info=FALSE,
                                               only_true=TRUE,
                                               missing=FALSE,
                                               percent=TRUE,
                                               dig=dig,
                                               dig_percent=dig_percent,
                                               add_lab = "")
                     
                     table_1$N <- table_1$N[2,]
                     table_1$N[1] <- "N"
                     
                     table_1$age <- table_wrap(X=data,
                                               varname="indexage",
                                               variable_name="Age",
                                               C=C,
                                               cn=cn,
                                               continuous=TRUE,
                                               cut_it=TRUE,
                                               breaks=c(0,50,70,80,Inf),
                                               labels=c("18-49","50-69","70-79","80+"),
                                               levels=NULL,
                                               add_info=TRUE,
                                               only_true=FALSE,
                                               missing=FALSE,
                                               percent=TRUE,
                                               dig=dig,
                                               dig_percent=dig_percent,
                                               add_lab = "")
                   
                     table_1$sex <- table_wrap(X=data,
                                                varname="sex",
                                                variable_name="Sex",
                                                C=C,
                                                cn=cn,
                                                continuous=FALSE,
                                                cut_it=TRUE,
                                                breaks=NULL,
                                                labels=c("Men","Women"),
                                                levels=c("man","woman"),
                                                add_info=FALSE,
                                                only_true=FALSE,
                                                missing=FALSE,
                                                percent=TRUE,
                                                dig=dig,
                                                dig_percent=dig_percent,
                                                add_lab = "")
                     
                     table_1$year <- table_wrap(X=data,
                                                varname="indexyear",
                                                variable_name="Index year",
                                                C=C,
                                                cn=cn,
                                                continuous=TRUE,
                                                cut_it=TRUE,
                                                breaks=c(0,2005,2012,2018,2022),
                                                labels=c("2003-2005","2006-2012","2013-2017","2018-2022"),
                                                levels=NULL,
                                                add_info=FALSE,
                                                only_true=FALSE,
                                                missing=FALSE,
                                                percent=TRUE,
                                                dig=dig,
                                                dig_percent=dig_percent,
                                                add_lab = "")
                     
                     table_1$bornse <- table_wrap(X=data %>% mutate(inse=Fodelseland_EU28 %in% "Sverige"),
                                                   varname="inse",
                                                   variable_name="Origin of birth",
                                                   C=C,
                                                   cn=cn,
                                                   continuous=FALSE,
                                                   cut_it=FALSE,
                                                   breaks=c(0,1,5,10,Inf),
                                                   labels=c("Sweden","Outside of Sweden"),
                                                   levels=c(TRUE,FALSE),
                                                   add_info=FALSE,
                                                   only_true=FALSE,
                                                   missing=FALSE,
                                                   percent=TRUE,
                                                   dig=dig,
                                                   dig_percent=dig_percent,
                                                   add_lab = "")
                     
                     table_1$yearsse <- table_wrap(X=data,
                                                varname="time_in_se",
                                                variable_name="Years in Sweden",
                                                C=C,
                                                cn=cn,
                                                continuous=TRUE,
                                                cut_it=TRUE,
                                                breaks=c(0,1,5,10,Inf),
                                                labels=c("0-1","1-5","5-10","10+"),
                                                levels=NULL,
                                                add_info=FALSE,
                                                only_true=FALSE,
                                                missing=FALSE,
                                                percent=TRUE,
                                                dig=dig,
                                                dig_percent=dig_percent,
                                                add_lab = "")
                     
                     table_1$HE10 <- table_wrap(X=data,
                                                varname="any10",
                                                variable_name="Number of healthcare encounters within 10 years",
                                                C=C,
                                                cn=cn,
                                                continuous=TRUE,
                                                cut_it=TRUE,
                                                breaks=c(0,1,4,Inf),
                                                labels=c("0","1-3",">=4"),
                                                levels=NULL,
                                                add_info=TRUE,
                                                only_true=FALSE,
                                                missing=FALSE,
                                                percent=TRUE,
                                                dig=dig,
                                                dig_percent=dig_percent,
                                                add_lab = "")[-2,]
                     
                     table_1$CCI <- table_wrap(X=data,
                                               varname="comorbidity_index_CCI10",
                                               variable_name="Charlson comorbidity index",
                                               C=C,
                                               cn=cn,
                                               continuous=TRUE,
                                               cut_it=TRUE,
                                               breaks=c(0,1,2,3,4,Inf),
                                               labels=c(0,1,2,3,"4+"),
                                               levels=NULL,
                                               add_info=FALSE,
                                               only_true=FALSE,
                                               missing=FALSE,
                                               percent=TRUE,
                                               dig=dig,
                                               dig_percent=dig_percent,
                                               add_lab = "")
                     
                     table_1$ECI <- table_wrap(X=data,
                                               varname="comorbidity_index_Elixhauser",
                                               variable_name="Elixhauser comorbidity index",
                                               C=C,
                                               cn=cn,
                                               continuous=TRUE,
                                               cut_it=TRUE,
                                               breaks=c(-Inf,0,1,6,14,Inf),
                                               labels=c("<0","0","1-5","6-13","14+"),
                                               levels=NULL,
                                               add_info=FALSE,
                                               only_true=FALSE,
                                               missing=FALSE,
                                               percent=TRUE,
                                               dig=dig,
                                               dig_percent=dig_percent,
                                               add_lab = "")
                     
                     # 0.88     1  1.11  1.64 for exp
                     data$MDCI <- exp(data$comorbidity_index_MDCI)
                     data$MDCI <- 1*(data$MDCI<0.88) + 2*(data$MDCI >= 0.88 & data$MDCI< 1) + 3*(data$MDCI >= 1 & data$MDCI<1.11) + 4*(data$MDCI>= 1.11 & data$MDCI<1.64)  +  5*(data$MDCI>=1.64)  
                     table_1$MDCI <- table_wrap(X=data,
                                               varname="MDCI",
                                               variable_name="Multidimensional diagnosis-based comorbidity index",
                                               C=C,
                                               cn=cn,
                                               continuous=FALSE,
                                               cut_it=FALSE,
                                               breaks=NULL,
                                               labels=NULL,
                                               levels=NULL,
                                               add_info=FALSE,
                                               only_true=FALSE,
                                               missing=FALSE,
                                               percent=TRUE,
                                               dig=dig,
                                               dig_percent=dig_percent,
                                               add_lab = "")
                     
                     # 0.881     1  1.11  1.64 for exp
                     data$MDCI <- exp(data$comorbidity_index_MDCI_sv)
                     data$MDCI <-  1*(data$MDCI<0.88) + 2*(data$MDCI >= 0.88 & data$MDCI< 1) + 3*(data$MDCI >= 1 & data$MDCI<1.11) + 4*(data$MDCI>= 1.11 & data$MDCI<1.64)  +  5*(data$MDCI>=1.64)  
                     table_1$MDCIsv <- table_wrap(X=data,
                                                  varname="MDCI",
                                                  variable_name="Multidimensional diagnosis-based comorbidity index (sv)",
                                                  C=C,
                                                  cn=cn,
                                                  continuous=FALSE,
                                                  cut_it=FALSE,
                                                  breaks=NULL,
                                                  labels=NULL,
                                                  levels=NULL,
                                                  add_info=FALSE,
                                                  only_true=FALSE,
                                                  missing=FALSE,
                                                  percent=TRUE,
                                                  dig=dig,
                                                  dig_percent=dig_percent,
                                                  add_lab = "")
                     
                     #   1  1.59  2.94  8.64 for exp
                     data$DCI <- exp(data$comorbidity_index_DCI)
                     data$DCI <- 1*(data$DCI<1) + 2*(data$DCI>=1 & data$DCI <1.59) + 3*(data$DCI>=1.59 & data$DCI< 2.94) + 4*(data$DCI>=2.94 & data$DCI < 8.64)  + 5*(data$DCI >= 8.64) 
             
                     table_1$DCI <- table_wrap(X=data,
                                               varname="DCI",
                                               variable_name="Drug comorbidity index",
                                               C=C,
                                               cn=cn,
                                               continuous=FALSE,
                                               cut_it=FALSE,
                                               breaks=NULL,
                                               labels=NULL,
                                               levels=NULL,
                                               add_info=FALSE,
                                               only_true=FALSE,
                                               missing=FALSE,
                                               percent=TRUE,
                                               dig=dig,
                                               dig_percent=dig_percent,
                                               add_lab = "")
                     
                     table_1[["n_chapters"]] <- table_wrap(X=data,
                                                           varname="n_chapters",
                                                           variable_name="Number of different ICD chapters",
                                                           C=C,
                                                           cn=cn,
                                                           continuous=TRUE,
                                                           cut_it=TRUE,
                                                           breaks=c(0,1,2,4,6,Inf),
                                                           labels=c(0,1,"2-3","4-5","6+"),
                                                           levels=NULL,
                                                           add_info=FALSE,
                                                           only_true=FALSE,
                                                           missing=FALSE,
                                                           percent=TRUE,
                                                           dig=dig,
                                                           dig_percent=dig_percent,
                                                           add_lab = "")
                     
                     n_unique <- length(unique(data$LopNr))
                     table_1$n_unique <- n_unique
                     
                     data2 <- data %>%
                       select(LopNr,sex) %>%
                       group_by(LopNr) %>%
                       slice(1) %>%
                       ungroup()
                     
                     table_1$n_unique_men <- sum(data2$sex %in% "man")
                     table_1$n_unique_women <- sum(data2$sex %in% "woman")
                     
                     rm(data) 
                     gc()
                     
                     return(table_1)
                     
                   })
  ret_list[[j]]<-ret
  rm(ret)
}

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)

###
# Combine and compute percents
tlist <- unlist(ret_list,recursive = FALSE)

n_unique <- unlist(lapply(FUN=function(x) x$n_unique,tlist))
sum(n_unique)

n_unique_men <- unlist(lapply(FUN=function(x) x$n_unique_men,tlist))
sum(n_unique_men)

n_unique_women <- unlist(lapply(FUN=function(x) x$n_unique_women,tlist))
sum(n_unique_women)


aged <- lapply(FUN=function(x) x$age[3,],tlist)
hes <- lapply(FUN=function(x) x$HE10[2,],tlist)

###
# Average the median IQR
aged <- do.call(rbind,aged)
hes <- do.call(rbind,hes)
openxlsx::write.xlsx(aged,file=paste0(".\\Results\\BaselineCharacteristics_aged.xlsx"))
openxlsx::write.xlsx(aged,file=paste0(".\\Results\\BaselineCharacteristics_hes.xlsx"))








tlist <- lapply(FUN=function(x){ x$n_unique <- NULL; x},tlist)
tlist <- lapply(FUN=function(x){ x$n_unique_men <- NULL; x},tlist)
tlist <- lapply(FUN=function(x){ x$n_unique_women <- NULL; x},tlist)

tlist <- lapply(FUN=function(x){ x$age <- x$age[-c(2:3),]; x},tlist)
tlist <- lapply(FUN=function(x){ x$HE10 <- x$HE10[-c(2),]; x},tlist)

table_1 <- do.call(rbind,tlist[[1]])
#table_1 <- table_1[-c(1:2),]
cns <- table_1[,c(1)]
cns[1]<-"N"
table_1 <- table_1[,c(2,4,6,8,10,12)]

conv_n <- function(x){
  isna <- x %in% ""
  y <- x[!isna]
  x <- rep(NA,length(x))
  y <- as.numeric(y)
  x[!isna] <- y
  return(x)
}

table_1 <- apply(FUN=conv_n,table_1,MARGIN=2)

for(j in 2:64){
  
  temp <- do.call(rbind,tlist[[j]])[,c(2,4,6,8,10,12)]
  #temp <- temp[-c(1:2),]
  temp <- apply(FUN=conv_n,temp,MARGIN=2)
  
  table_1 <- table_1 + temp
}

table_1p <- table_1

for(j in 1:ncol(table_1p)){
  table_1p[,j] <- round(100*table_1p[,j]/table_1p[1,j],dig=1)
}

table_1_f <- matrix(NA,ncol=ncol(table_1)*2 + 1,nrow=nrow(table_1))

table_1_f[,1] <- cns

fix2 <- function(x){
  x <- x / 100000
  for(j in 1:length(x)){
    
    if(!is.na(x[j])){
      if(x[j]<1){
        x[j] <- format(round(x[j],dig=1),nsmall=1)
      } else{
        x[j] <- round(x[j]) 
      }
    }
    
  }
  return(x)
}

table_1_f[,c(2,4,6,8,10,12)] <- apply(FUN=fix2,table_1,MARGIN=c(1,2))

table_1_f[,c(2,4,6,8,10,12)+1] <- apply(FUN=function(x){
  sapply(FUN=function(xx){
    if(is.na(xx)){
      ret <- ""
    } else{
      ret <- paste0("(",xx,")")
    }
    return(ret)
  },x)
  },table_1p,MARGIN=2)

table_1_f[is.na(table_1_f)]<-""



###
# Write to file
openxlsx::write.xlsx(table_1_f,file=paste0(".\\Results\\BaselineCharacteristics.xlsx"))



###
# End
###