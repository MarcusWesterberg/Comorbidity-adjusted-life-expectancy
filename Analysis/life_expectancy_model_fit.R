###
# Survival models 

X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
require(tidyverse)

log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"
background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")
figs_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")

extract_obs <- function(timefu,censor){
  
  sf <- survfit(Surv(timefu, censor) ~ 1)
  
  sf$surv[sf$time==max(sf$time)]
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Crude model first

require(parallel)

start_time <- Sys.time()
cluster_indices <- 18:100
cluster_ids <- 1:64

print(cluster_ids)
ncores <- 32

print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("X","data_path","save_path","cluster_ids"),envir = environment())

ret <- parLapplyLB(cl = cl,
                   cluster_indices,
                   fun = function(a){
                     sink(paste0("C:\\Marcus\\Misc\\","LEcox_age_crude",a,".txt"))
                     require(dplyr)
                     X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
                     
                     source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\","est_model_function.R"))
                     
                     print(a)
                     print(cluster_ids)
                     
                     print("Loading data...")
                     
                     data_list <- list()
                     for(id in cluster_ids){
                       print(id)
                       print(Sys.time())
                       load(file=paste0("C:\\Marcus\\ComorbidityBase_cache\\cox_life_expectancy_",a,"_",id,".Rdata"))
                       data_list[[id]] <- data_temp
                     }
                     
                     data <- do.call(rbind,data_list)
                     rm(data_temp)
                     rm(data_list)
                     gc()
                  
                     data$indexage <- data$indexage - a - 0.5
             
                     require(survival)
                     require(pec)
                     require(flexsurv)
                     
               
                     coeff_list <- list()
                     current_inits <- NULL
                     maxfu <- 10
                     
                     for(sexstrata in c("man","woman")){
                       
                       cat("\n")
                       print(Sys.time())
                       print(sexstrata)
                       print(a)
                       print("Finding rate and shape")
                       
                       current_inits <- c(0.05, log(1/mean(data$timefu[data$sex %in% sexstrata])) )  
                       
                       print(paste0("Estimating model ",a,", sex ",sexstrata,", maxfu ",maxfu))
                       
                       m.init <- est_model(data=data,
                                           sexstrata=sexstrata,
                                           formula=as.formula("Surv(timefu, censor) ~ 1"),
                                           maxfu=maxfu,
                                           name="init",
                                           init_values=current_inits)
                       
                       init_values <- m.init$coefficients
                       init_values[2] <- exp( init_values[2] )
                       
                       model <- m.init
                       rm(m.init)
                       
                       modelname <- paste0(a,"_","CRUDE",
                                           "",
                                           maxfu,
                                           "",
                                           sexstrata)
                       
                       # model.frame(model)
                       
                       save(model,file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                       coeff_list[[modelname]] <- model$coefficients
                       rm(model)
                       
                       print("Initial values")
                       print(init_values)
                    
                       cat("\n")
                     }
                     
                     rmlist <- ls()
                     rmlist <- rmlist[!rmlist %in% "coeff_list"]
                     rm(list=rmlist)
                     gc()
                     
                     return(coeff_list)
                   })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)

save(ret,file=paste0(figs_path,"parameter.values.crude.Rdata"))

load(file=paste0(figs_path,"parameter.values.crude.Rdata"))

length(ret[[1]])
names(ret[[1]])

## 
estlist_crude_m <- do.call(rbind,lapply(FUN=function(x) x[[1]],ret))
estlist_crude_w <- do.call(rbind,lapply(FUN=function(x) x[[2]],ret))

a <- 18:100
w <- a
w[w>90]<-90
w[w<30]<-30
w <- w^2
w <- w/sum(w)
smooth_it <- function(y){
  mgcv::gam(formula=y ~ s(a,k=25,bs="cr"),family = "gaussian",method="REML",weights = w)$fitted.values
}
smooth_it_wrap <- function(y){
  apply(FUN=smooth_it,y,MARGIN=2)
}

estlist_crude_m_s <- smooth_it_wrap(estlist_crude_m)
estlist_crude_w_s <- smooth_it_wrap(estlist_crude_w)

# CRUDE
par(mfrow=c(2,2))
ages <- 18:100
for(j in 1:2){
  plot(ages,estlist_crude_m[,j],main=colnames(estlist_crude_m)[j])
  lines(ages,estlist_crude_m_s[,j],col="red")
  #m <- lm(estlist_mdci_m[,j]~ages + I(ages^2))
  #print(m)
  #lines(ages,m$fitted.values,col="red")
}
for(j in 1:2){
  plot(ages,estlist_crude_w[,j],main=colnames(estlist_crude_w)[j])
  #m <- lm(estlist_mdci_w[,j]~ages + I(ages^2))
  lines(ages,estlist_crude_w_s[,j],col="red")
  ##lines(ages,m$fitted.values,col="red")
}

rownames(estlist_crude_m_s)<- ages
rownames(estlist_crude_w_s)<- ages

parameter.values.smooth <- list("men"=estlist_crude_m_s,
                                "women"=estlist_crude_w_s)
save(parameter.values.smooth,file=paste0(figs_path,"parameter.values.crude.smooth.Rdata"))

















### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# MDCI DCI model
require(parallel)

start_time <- Sys.time()
cluster_indices <- 18:100
cluster_ids <- 1:64

print(cluster_ids)
ncores <- 32

print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("X","data_path","save_path","cluster_ids"),envir = environment())

ret <- parLapplyLB(cl = cl,
                   cluster_indices,
                   fun = function(a){
                     sink(paste0("C:\\Marcus\\Misc\\","LEcox_age_mdcidci",a,".txt"))
                     require(dplyr)
                     X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
                     
                     source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\","est_model_function.R"))
                     
                     print(a)
                     print(cluster_ids)
                     
                     print("Loading data...")
                     
                     data_list <- list()
                     for(id in cluster_ids){
                       print(id)
                       print(Sys.time())
                       load(file=paste0("C:\\Marcus\\ComorbidityBase_cache\\cox_life_expectancy_",a,"_",id,".Rdata"))
                       data_list[[id]] <- data_temp
                     }
                     
                     data <- do.call(rbind,data_list)
                     rm(data_temp)
                     rm(data_list)
                     gc()
                     
                     
                     trunc_it <- function(x,l,u){
                       x[x<l]<-l
                       x[x>u]<-u
                       return(x)
                     }
                     
                     MDCI_cutoffs <- c(-1,4)
                     DCI_cutoffs <- c(-0.5,8)
                     
                     
                     data$indexage <- data$indexage - a - 0.5
                     
                     data$MDCI <- trunc_it(data$MDCI,
                                           l=MDCI_cutoffs[1],
                                           u=MDCI_cutoffs[2])
                     
                     data$dci <- trunc_it(data$dci,
                                          l=DCI_cutoffs[1],
                                          u=DCI_cutoffs[2])
                     
                     ###
                     # Assess impact of follow-up
                     formulas <- list("MDCIDCI"=as.formula("Surv(timefu, censor) ~  MDCI + dci + I(MDCI*dci)"))
                     
                     require(survival)
                     require(pec)
                     require(flexsurv)
                     
                     
                     coeff_list <- list()
                     current_inits <- NULL
                     maxfu <- 10
                     
                     for(sexstrata in c("man","woman")){
                       
                       
                       cat("\n")
                       print(Sys.time())
                       print(sexstrata)
                       print("Finding rate and shape")
                       
                       modelname <- paste0(a,"_","CRUDE",
                                           "",
                                           maxfu,
                                           "",
                                           sexstrata)
                       
                       load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                       
                       init_values <- model$coefficients
                       init_values[2] <- exp( init_values[2] )
                       
                       print("Initial values")
                       print(init_values)
                       rm(model)
                       cat("\n")
                       
                       print(paste0("Estimating model ",a,", sex ",sexstrata,", maxfu ",maxfu))
                       for(k in 1:length(formulas)){
                         cat("\n")
                         modelname <- paste0(a,"_",names(formulas)[k],
                                             "",
                                             maxfu,
                                             "",
                                             sexstrata)
                         print(modelname)
                         
                         
                         model <- est_model(data=data,
                                            sexstrata=sexstrata,
                                            formula=formulas[[k]],
                                            maxfu=maxfu,
                                            name=names(formulas)[k],
                                            init_values=init_values)
                         coeff_list[[modelname]] <- model$coefficients
                         save(model,file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                         
                         rm(model)
                         gc()
                         cat("\n")
                       }
                     }
                     
                     rmlist <- ls()
                     rmlist <- rmlist[! rmlist %in% "coeff_list"]
                     rm(list=rmlist)
                     gc()
                     
                     return(coeff_list)
                   })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)

save(ret,file=paste0(figs_path,"parameter.values.mdcidci.Rdata"))

load(file=paste0(figs_path,"parameter.values.mdcidci.Rdata"))

length(ret[[1]])
names(ret[[1]])

## 
estlist_mdci_m <- do.call(rbind,lapply(FUN=function(x) x[[1]],ret))

estlist_mdci_w <- do.call(rbind,lapply(FUN=function(x) x[[2]],ret))


a<- 18:100
smooth_it <- function(y){
  mgcv::gam(formula=y ~ s(a,k=25,bs="cr"),family = "gaussian",method="REML")$fitted.values
}
smooth_it_wrap <- function(y){
  apply(FUN=smooth_it,y,MARGIN=2)
}

estlist_mdci_m_s <- smooth_it_wrap(estlist_mdci_m)
estlist_mdci_w_s <- smooth_it_wrap(estlist_mdci_w)

# MDCI DCI
par(mfrow=c(2,5))
ages <- 18:100
for(j in 1:5){
  plot(ages,estlist_mdci_m[,j],main=colnames(estlist_mdci_m)[j])
  lines(ages,estlist_mdci_m_s[,j],col="red")

}
for(j in 1:5){
  plot(ages,estlist_mdci_w[,j],main=colnames(estlist_mdci_w)[j])

  lines(ages,estlist_mdci_w_s[,j],col="red")

}

rownames(estlist_mdci_m_s)<- ages
rownames(estlist_mdci_w_s)<- ages

parameter.values.smooth <- list("men"=estlist_mdci_m_s,
                                "women"=estlist_mdci_w_s)
save(parameter.values.smooth,file=paste0(figs_path,"parameter.values.mdcidci.smooth.Rdata"))












### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# CCI model
require(parallel)

start_time <- Sys.time()
cluster_indices <- 18:100
cluster_ids <- 1:64

print(cluster_ids)
ncores <- 32

print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("X","data_path","save_path","cluster_ids"),envir = environment())

ret <- parLapplyLB(cl = cl,
                   cluster_indices,
                   fun = function(a){
                     sink(paste0("C:\\Marcus\\Misc\\","LEcox_age_cci",a,".txt"))
                     require(dplyr)
                     require(splines)
                     X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
                     
                     source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\","est_model_function.R"))
                     
                     print(a)
                     print(cluster_ids)
                     
                     print("Loading data...")
                     
                     data_list <- list()
                     for(id in cluster_ids){
                       print(id)
                       print(Sys.time())
                       load(file=paste0("C:\\Marcus\\ComorbidityBase_cache\\cox_life_expectancy_",a,"_",id,".Rdata"))
                       data_list[[id]] <- data_temp
                     }
                     
                     data <- do.call(rbind,data_list)
                     rm(data_temp)
                     rm(data_list)
                     gc()
                     
                     
                     trunc_it <- function(x,l,u){
                       x[x<l]<-l
                       x[x>u]<-u
                       return(x)
                     }
                     
                  
                     data$CCI[is.na(data$CCI)]<-0
                     data$CCI[data$CCI>6] <- 6
                     
                     
                     print(table(data$CCI,data$sex))
                     
                     data$indexage <- data$indexage - a - 0.5
                   
                     
                     ###
                     # Assess impact of follow-up
                     formulas <- list("CCI"=as.formula("Surv(timefu, censor) ~ ns(CCI, knots = c(1, 3), Boundary.knots = c(0,4), intercept = FALSE)"))
                     
                     require(survival)
                     require(pec)
                     require(flexsurv)
                     
                     
                     coeff_list <- list()
                     current_inits <- NULL
                     
                     maxfu <- 10
                     for(sexstrata in c("man","woman")){
                       
                       
                       
                       cat("\n")
                       print(Sys.time())
                       print(sexstrata)
                       print("Finding rate and shape")
                       
                       modelname <- paste0(a,"_","CRUDE",
                                           "",
                                           maxfu,
                                           "",
                                           sexstrata)
                       
                       load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                       
                       init_values <- model$coefficients
                       init_values[2] <- exp( init_values[2] )
                       
                       print("Initial values")
                       print(init_values)
                       rm(model)
                       cat("\n")
                       
                       print(paste0("Estimating model ",a,", sex ",sexstrata,", maxfu ",maxfu))
                       for(k in 1:length(formulas)){
                         cat("\n")
                         modelname <- paste0(a,"_",names(formulas)[k],
                                             "",
                                             maxfu,
                                             "",
                                             sexstrata)
                         print(modelname)
                         
                         
                         model <- est_model(data=data,
                                            sexstrata=sexstrata,
                                            formula=formulas[[k]],
                                            maxfu=maxfu,
                                            name=names(formulas)[k],
                                            init_values=init_values)
                         coeff_list[[modelname]] <- model$coefficients
                         save(model,file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                         
                         rm(model)
                         gc()
                         cat("\n")
                       }
                       
                     }
                     
                     rmlist <- ls()
                     rmlist <- rmlist[! rmlist %in% "coeff_list"]
                     rm(list=rmlist)
                     gc()
                     
                     return(coeff_list)
                   })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)

save(ret,file=paste0(figs_path,"parameter.values.cci.Rdata"))

loads <- load(file=paste0(figs_path,"parameter.values.cci.Rdata"))

length(ret[[1]])
names(ret[[1]])

## 
estlist_cci_m <- do.call(rbind,lapply(FUN=function(x) x[[1]],ret))

estlist_cci_w <- do.call(rbind,lapply(FUN=function(x) x[[2]],ret))


a <- 18:100
smooth_it <- function(y){
  mgcv::gam(formula=y ~ s(a,k=25,bs="cr"),family = "gaussian",method="REML")$fitted.values
}
smooth_it_wrap <- function(y){
  apply(FUN=smooth_it,y,MARGIN=2)
}

estlist_cci_m_s <- smooth_it_wrap(estlist_cci_m)
estlist_cci_w_s <- smooth_it_wrap(estlist_cci_w)


par(mfrow=c(2,5),mar=c(2,3,2,0.5))
ages <- 18:100
for(j in 1:5){
  plot(ages,estlist_cci_m[,j],main=colnames(estlist_cci_m)[j])
  lines(ages,estlist_cci_m_s[,j],col="red")

}
for(j in 1:5){
  plot(ages,estlist_cci_w[,j],main=colnames(estlist_cci_w)[j])
  lines(ages,estlist_cci_w_s[,j],col="red")
 
}
rownames(estlist_cci_m_s)<- ages
rownames(estlist_cci_w_s)<- ages

parameter.values.smooth <- list("men"=estlist_cci_m_s,
                                "women"=estlist_cci_w_s)
save(parameter.values.smooth,file=paste0(figs_path,"parameter.values.cci.smooth.Rdata"))





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# CCI + DCI model
require(parallel)

start_time <- Sys.time()
cluster_indices <- 18:100
cluster_ids <- 1:64

print(cluster_ids)
ncores <- 32

print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("X","data_path","save_path","cluster_ids"),envir = environment())

ret <- parLapplyLB(cl = cl,
                   cluster_indices,
                   fun = function(a){
                     sink(paste0("C:\\Marcus\\Misc\\","LEcox_age_ccidci",a,".txt"))
                     require(dplyr)
                     require(splines)
                     X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
                     
                     source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\","est_model_function.R"))
                     
                     print(a)
                     print(cluster_ids)
                     
                     print("Loading data...")
                     
                     data_list <- list()
                     for(id in cluster_ids){
                       print(id)
                       print(Sys.time())
                       load(file=paste0("C:\\Marcus\\ComorbidityBase_cache\\cox_life_expectancy_",a,"_",id,".Rdata"))
                       data_list[[id]] <- data_temp
                     }
                     
                     data <- do.call(rbind,data_list)
                     rm(data_temp)
                     rm(data_list)
                     gc()
                     
                     
                     trunc_it <- function(x,l,u){
                       x[x<l]<-l
                       x[x>u]<-u
                       return(x)
                     }
                     
                     
                     DCI_cutoffs <- c(-0.5,8)
                   
                     
                     data$dci <- trunc_it(data$dci,
                                          l=DCI_cutoffs[1],
                                          u=DCI_cutoffs[2])
                     
                     
                     data$CCI[is.na(data$CCI)]<-0
                     data$CCI[data$CCI>6] <- 6
                     
                     
                     print(table(data$CCI,data$sex))
                     
                     data$indexage <- data$indexage - a - 0.5
                     
                     
                     ###
                     # Assess impact of follow-up
                     formulas <- list("CCIDCI"=as.formula("Surv(timefu, censor) ~ ns(CCI, knots = c(1, 3), Boundary.knots = c(0,4), intercept = FALSE)*dci"))
                     
                     require(survival)
                     require(pec)
                     require(flexsurv)
                     
                     
                     coeff_list <- list()
                     current_inits <- NULL
                     
                     maxfu <- 10
                     for(sexstrata in c("man","woman")){
                       
                       
                       
                       cat("\n")
                       print(Sys.time())
                       print(sexstrata)
                       print("Finding rate and shape")
                       
                       modelname <- paste0(a,"_","CRUDE",
                                           "",
                                           maxfu,
                                           "",
                                           sexstrata)
                       
                       load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                       
                       init_values <- model$coefficients
                       init_values[2] <- exp( init_values[2] )
                       
                       print("Initial values")
                       print(init_values)
                       rm(model)
                       cat("\n")
                       
                       print(paste0("Estimating model ",a,", sex ",sexstrata,", maxfu ",maxfu))
                       for(k in 1:length(formulas)){
                         cat("\n")
                         modelname <- paste0(a,"_",names(formulas)[k],
                                             "",
                                             maxfu,
                                             "",
                                             sexstrata)
                         print(modelname)
                         
                         
                         model <- est_model(data=data,
                                            sexstrata=sexstrata,
                                            formula=formulas[[k]],
                                            maxfu=maxfu,
                                            name=names(formulas)[k],
                                            init_values=init_values)
                         coeff_list[[modelname]] <- model$coefficients
                         save(model,file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                         
                         rm(model)
                         gc()
                         cat("\n")
                       }
                       
                     }
                     
                     rmlist <- ls()
                     rmlist <- rmlist[! rmlist %in% "coeff_list"]
                     rm(list=rmlist)
                     gc()
                     
                     return(coeff_list)
                   })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)

save(ret,file=paste0(figs_path,"parameter.values.ccidci.Rdata"))

loads <- load(file=paste0(figs_path,"parameter.values.ccidci.Rdata"))

length(ret[[1]])
names(ret[[1]])

## 
estlist_ccidci_m <- do.call(rbind,lapply(FUN=function(x) x[[1]],ret))

estlist_ccidci_w <- do.call(rbind,lapply(FUN=function(x) x[[2]],ret))


a <- 18:100
smooth_it <- function(y){
  mgcv::gam(formula=y ~ s(a,k=25,bs="cr"),family = "gaussian",method="REML")$fitted.values
}
smooth_it_wrap <- function(y){
  apply(FUN=smooth_it,y,MARGIN=2)
}

estlist_ccidci_m_s <- smooth_it_wrap(estlist_ccidci_m)
estlist_ccidci_w_s <- smooth_it_wrap(estlist_ccidci_w)


par(mfrow=c(3,3),mar=c(2,3,2,0.5))
ages <- 18:100
for(j in 1:9){
  plot(ages,estlist_ccidci_m[,j],main=colnames(estlist_ccidci_m)[j])
  lines(ages,estlist_ccidci_m_s[,j],col="red")
  
}
for(j in 1:9){
  plot(ages,estlist_ccidci_w[,j],main=colnames(estlist_ccidci_w)[j])
  lines(ages,estlist_ccidci_w_s[,j],col="red")
  
}
rownames(estlist_ccidci_m_s)<- ages
rownames(estlist_ccidci_w_s)<- ages

parameter.values.smooth <- list("men"=estlist_ccidci_m_s,
                                "women"=estlist_ccidci_w_s)
save(parameter.values.smooth,file=paste0(figs_path,"parameter.values.ccidci.smooth.Rdata"))





