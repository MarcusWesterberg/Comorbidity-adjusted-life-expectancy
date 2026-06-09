###
# Folow-up information
X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

ncores <- 16
data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\"))

###
# Load data...

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
                     setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment CCI\\"))
                     
                     cluster_id <<- cluster_id # make seed global
                     set.seed(cluster_id)
                     source(paste0(".\\Scripts\\Table_wrap.R"))
                     print("Loading data...")
                     Sys.sleep(cluster_id/16)
                     
                     
                     
                     #loads <- load(file=paste0(data_path,"charlson10_cache_part",cluster_id,".Rdata"))
                     #cci_data <- data
                     
                     loads <- load(file=paste0(data_path,"cindex_cache_part",cluster_id,".Rdata"))
                     
                     data <- data %>%
                       select(LopNr,indexdate,indexyear,indexage,timefu,censor)
                     gc()
                     
                     # data <- data %>%
                     #   left_join(cci_data %>% 
                     #               select(-any_of(c("rownr","CCIw","CCIunw"))),
                     #             by=c("LopNr","indexdate"))
                     # rm(cci_data)
                     # gc()
                     
                     data <- data %>%
                       filter(indexyear<2023 & indexyear >= 2006) %>%
                       filter(indexage >= 18 & indexage < 121)
                     
                     # data <- data %>%
                     #   group_by(LopNr) %>%
                     #   slice(sample(1:n(),size=1)) %>%
                     #   ungroup()
                     
                     ###
                     # keep only needed coluns
                     
                     data$indexage <- floor(data$indexage)
                     
                     # define columns
                     C = cbind(TRUE,
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
                     
                     require(survival)
                     revKM <- c()
                     revKMq1 <- c()
                     revKMq3 <- c()
                     for(j in 1:ncol(C)){
                       print(j)
                       temp <- data[C[,j],]
                       m <- survfit(Surv(timefu, censor) ~ 1,
                                    data=temp %>% mutate(censor=1-censor))
                       
                       q <- quantile(m, probs = c(0.5,0.25 , 0.75),conf.int = FALSE)
                       revKM[j]<-q[1]
                       revKMq1[j]<-q[2]
                       revKMq3[j]<-q[3]
                     }
                     
                     revKM <- rbind(revKM,revKMq1,revKMq3)
                     
                     attr(table_1,"revKM") <- revKM
                     
                     table_1$timefu <- table_wrap(X=data,
                                                  varname="timefu",
                                                  variable_name="Follow-up",
                                                  C=C,
                                                  cn=cn,
                                                  continuous=TRUE,
                                                  cut_it=TRUE,
                                                  breaks=c(0,1,5,10,Inf),
                                                  labels=c("<=1","1-5","5-10",">10"),
                                                  levels=NULL,
                                                  add_info=TRUE,
                                                  only_true=FALSE,
                                                  missing=FALSE,
                                                  percent=TRUE,
                                                  dig=dig,
                                                  dig_percent=dig_percent,
                                                  add_lab = "")
                     
                     table_1$censor <- table_wrap(X=data %>% mutate(censor = censor %in% 1),
                                                  varname="censor",
                                                  variable_name="Deaths",
                                                  C=C,
                                                  cn=cn,
                                                  continuous=FALSE,
                                                  cut_it=FALSE,
                                                  labels=c("Yes","Censored"),
                                                  levels=c(TRUE,FALSE),
                                                  add_info=FALSE,
                                                  only_true=FALSE,
                                                  missing=FALSE,
                                                  percent=TRUE,
                                                  dig=dig,
                                                  dig_percent=dig_percent,
                                                  add_lab = "")
                     
              
                     
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

revKMlist <- lapply(FUN=function(x) attr(x,"revKM"),tlist)
revKM <- revKMlist[[1]]
for(j in 2:length(tlist)){
  revKM <- revKM + revKMlist[[j]]
}

revKM <- revKM/length(tlist)
revKM <- round(revKM,dig=1)


table_1 <- do.call(rbind,tlist[[1]])
cns <- table_1[,c(1)]
table_1 <- table_1[,c(2,4,6,8)]


conv_n <- function(x){
  isna <- x %in% "" | grepl(x,pattern="[(]")
  y <- x[!isna ]
  
  # y2 <- y[grepl(y,pattern="[(]")]
  # 
  # y2a <- y2[!grepl(y2,pattern="[-]")]
  # y2b <- y2[grepl(y2,pattern="[-]")]
  # 
  # if(length(y2a)>0){
  #   y2a <- strsplit(y2a,split=" ")
  # }
  # if(length(y2b)>0){
  #   y2b <- strsplit(y2b,split=" ")
  # }
  # 
  x <- rep(NA,length(x))
  y <- as.numeric(y)
  x[!isna] <- y
  return(x)
}

table_1 <- apply(FUN=conv_n,table_1,MARGIN=2)
for(j in 2:64){
  
  temp <- do.call(rbind,tlist[[j]])[,c(2,4,6,8)]
  
  temp <- apply(FUN=conv_n,temp,MARGIN=2)
  
  table_1 <- table_1 + temp
}

table_1p <- table_1

for(j in 1:ncol(table_1p)){
  table_1p[1:11,j] <- round(100*table_1p[1:11,j]/table_1p[1,j],dig=1)
  #table_1p[12:22,j] <- round(100*table_1p[12:22,j]/table_1p[12,j],dig=1)
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

table_1_f[,c(2,4,6,8)] <- apply(FUN=fix2,table_1,MARGIN=c(1,2))

table_1_f[,c(2,4,6,8)+1] <- apply(FUN=function(x){
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
openxlsx::write.xlsx(table_1_f,file=paste0(".\\Results\\Followup.xlsx"))
openxlsx::write.xlsx(as.data.frame(revKM),file=paste0(".\\Results\\Followup_revKM.xlsx"),rowNames=TRUE)




###
# End
###