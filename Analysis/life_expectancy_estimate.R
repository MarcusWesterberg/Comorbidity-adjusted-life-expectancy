###
# Survival models 

X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
require(tidyverse)

log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"
background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")
figs_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")

###
# Load data...
require(parallel)

start_time <- Sys.time()
cluster_indices <- 18:100 # 18:100
cluster_ids <- 1:64

print(cluster_ids)
ncores <- 43

maxfu <- 10

print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("X","data_path","figs_path","save_path","cluster_ids","maxfu"),envir = environment())

ret <- parLapplyLB(cl = cl,
                   cluster_indices,
                   fun = function(a){
                     require(dplyr)
                     require(splines)
                     sink(paste0("C:\\Marcus\\Misc\\","LE_est_age",a,".txt"))
                     print(a)
                     print(cluster_ids)
                     print(Sys.time())
                     
                     print("Loading data...")
                     
                     
                     parameter.values.smooth.crude <- get(load(file=paste0(figs_path,"parameter.values.crude.smooth.Rdata")))
                     parameter.values.smooth.cci <- get(load(file=paste0(figs_path,"parameter.values.cci.smooth.Rdata")))
                     parameter.values.smooth.mdcidci <- get(load(file=paste0(figs_path,"parameter.values.mdcidci.smooth.Rdata")))
                     
                     names(parameter.values.smooth.crude) <- c("man","woman")
                     names(parameter.values.smooth.cci) <- c("man","woman")
                     names(parameter.values.smooth.mdcidci) <- c("man","woman")
                     parameter.values.smooth <- list("CRUDE"=parameter.values.smooth.crude,
                                                     "CCI"=parameter.values.smooth.cci,
                                                     "MDCIDCI"=parameter.values.smooth.mdcidci)
                     
                     ###
                     # Load individual-level data 
                     data_list_men <- list()
                     data_list_women <- list()
                     for(id in cluster_ids){
                       print(id)
                       print(Sys.time())
                       loads <- load(file=paste0("C:\\Marcus\\ComorbidityBase_cache\\cox_life_expectancy_",a,"_",id,".Rdata"))
                       
                       trunc_it <- function(x,l,u){
                         x[x<l]<-l
                         x[x>u]<-u
                         return(x)
                       }
                       
                       MDCI_cutoffs <- c(-1,4)
                       DCI_cutoffs <- c(-0.5,8)
                       
                       data_temp$CCI[is.na(data_temp$CCI)]<-0
                       data_temp$CCI[data_temp$CCI>6] <- 6
                       
                       data_temp$indexage <- data_temp$indexage - a - 0.5
                   
                       
                       data_temp$MDCI <- trunc_it(data_temp$MDCI,
                                                  l=MDCI_cutoffs[1],
                                                  u=MDCI_cutoffs[2])
                       
                       data_temp$dci <- trunc_it(data_temp$dci,
                                                 l=DCI_cutoffs[1],
                                                 u=DCI_cutoffs[2])
                       
                       
                       data_list_men[[id]] <- data_temp %>% 
                         filter(sex %in% "man") %>%
                         select(LopNr,rowid,timefu,censor,sex,indexdate,indexage,MDCI,dci,CCI,age)
                       
                       data_list_women[[id]] <- data_temp %>% 
                         filter(sex %in% "woman") %>%
                         select(LopNr,rowid,timefu,censor,sex,indexdate,indexage,MDCI,dci,CCI,age)
                     }
                     
                     rm(data_temp)
                     gc()
                     
                     
                     
                     extract_obs <- function(timefu,censor){
                       
                       sf <- survfit(Surv(timefu, censor) ~ 1)
                       
                       sf$surv[sf$time==max(sf$time)]
                     }
               
                     # extrapolates and provides life expectancy = integral of survival function
                     extrapolate_cumhaz <- function(cumhaz,
                                                    age,
                                                    scb_hazard,
                                                    plotit=FALSE){
                       
                       
                       ###
                       # Find average hazard 
                       scb_hazard_avr <- scb_hazard
                       inds <- scb_hazard_avr$x>= age
                       scb_hazard_avr$x <- scb_hazard_avr$x[inds]
                       scb_hazard_avr$y <- scb_hazard_avr$y[inds]
                       
                       ###
                       # Hazard
                       haz <- diff(cumhaz$.pred_cumhaz)
                       
                       # extend until end of scb_hazard_avr
                       haz <- c(haz,rep(0,length(scb_hazard_avr$y)-length(haz)))
                       
                       stopifnot(length(haz)==length(scb_hazard_avr$y))
                       
                       # interpolation between 10-30
                       min.length <- 30
                       interpol.window.start <- min(10,110-age)
                       interpol.window.end <- min(min.length, (110-age) ) 
                       
                       alpha <- rep(0,length(haz))
                       alpha[interpol.window.start:interpol.window.end] <- seq(0,1,length.out=interpol.window.end-interpol.window.start+1)
                       alpha[interpol.window.end:length(alpha)]<-1
                       
                       hazard_interpol <- haz*(1-alpha) + scb_hazard_avr$y*alpha
                       
                       survcurv <- c(1,exp(-cumsum(hazard_interpol))) # add time zero (surv=1)
                       
                       le <- (1 / 2) * (2*sum(survcurv) - survcurv[1] - survcurv[length(survcurv)]) 
                       
                       if(plotit){
                         
                         xx <- seq(0,length(survcurv)-1,by=1) + a
                         plot(xx,survcurv,col="darkgrey",type="l",lwd=2,xlab="Time",ylab="Survival")
                         
                         points( 1:interpol.window.start -1 + a,
                                 survcurv[1:interpol.window.start],
                                 col="black",pch=20,cex=2)
                         points( (interpol.window.start+1):interpol.window.end -1 + a,
                                 survcurv[ (interpol.window.start+1):interpol.window.end],
                                 col="red",pch=20,cex=2)
                         
                         mtext(side=3,text=paste0("LE: ",round(le,1) ),line=1)
                       }
                       
                       return(le)
                     }
                     
                     
                     
                     require(flexsurv)
                     require(pec)
                     require(dplyr)
                     
                     
                     ###
                     # Extract hazard functions from scb data
                     scb_hazard_men <- get(load(file=paste0(save_path,"crude_hazard_men.Rdata"))) # scb_hazard_men or crude age based model
                     scb_hazard_women <- get(load(file=paste0(save_path,"crude_hazard_women.Rdata"))) # scb_hazard_women or crude age based model
                     
                     # Hazard according to MDCI and DCI 
                     at.times <- seq(1,31,by=1)
                     
                     find_le <- function(model,data,sex,name,scb_hazard,plotit,simple=FALSE){
                       
                       if(simple){
                         model.h <- predict(model,
                                            newdata=data %>% slice(1),
                                            type="cumhaz",
                                            times=at.times)
                       } else{
                         model.h <- predict(model,
                                            newdata=data,
                                            type="cumhaz",
                                            times=at.times)
                       }
            
                       
                       if(plotit & !simple){
                         pdf(file=paste0("C:\\Marcus\\LEcox\\","Preds_",name,"_",sex,"_",a,".pdf"),width=8,height=8)
                         z <- data$MDCI
                         
                         qs <- quantile(z,probs=c(0.001,0.05,0.5,0.95,0.999))
                         
                         plot_these <- sapply(FUN=function(x) which(z %in% max(z[z<=x]))[1] ,qs)
                         plot_these <- c(plot_these,which.min(abs(mean(z)-z)))
                         par(mfrow=c(2,3))
                         
                         for(k in 1:6){
                           extrapolate_cumhaz(cumhaz=model.h[[1]][[plot_these[k]]],
                                              age=a,
                                              scb_hazard=scb_hazard,
                                              plotit=TRUE)
                         }
                         dev.off()
                       }
                       
                       
                       le <- lapply(FUN=extrapolate_cumhaz,
                                    X=model.h[[1]],
                                    age=a,
                                    scb_hazard=scb_hazard)
                       le <- unlist(le)
                       
                       print(summary(le))
                       if(simple){
                         le <- rep(le,nrow(data))
                       }
                  
                       # quantile(le.men,probs=seq(0,1,0.05))
                       rm(model.h)
                       gc()
                       return(le)
                     }
                     
                    
                     ###
                     # Load models
                     print("CRUDE men")
                     print(Sys.time())
                     modelname <- paste0(a,"_","CRUDE",
                                         "",
                                         maxfu,
                                         "",
                                         "man")
                     
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                     temp <-  parameter.values.smooth[["CRUDE"]][["man"]]
                     
                     model$coefficients <- temp[rownames(temp)==a,]
                     model$res[,1] <- temp[rownames(temp)==a,]
                     model$res[2,1]<-exp(  model$res[2,1])
                     
                     le.men.crude <- list()
                     for(i in 1:length(data_list_men)){
                       print(i)
                       le.men.crude[[i]] <- find_le(model=model,
                                                    data=data_list_men[[i]] %>% 
                                                      select(LopNr,sex,age),
                                                    sex="men",
                                                    name="CRUDE",
                                                    scb_hazard=scb_hazard_men,
                                                    plotit= FALSE,
                                                    simple=TRUE)
                     }
                     rm(model)
                     
                     
                  
                
                     #modelname CRUDE MDCIDCI CCI
                     print("MDCIDCI men")
                     print(Sys.time())
                     modelname <- paste0(a,"_","MDCIDCI",
                                         "",
                                         maxfu,
                                         "",
                                         "man")
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                     temp <-  parameter.values.smooth[[ "MDCIDCI" ]][["man"]]
                     
                     model$coefficients <- temp[rownames(temp)==a,]
                     model$res[,1] <- temp[rownames(temp)==a,]
                     model$res[2,1]<-exp(  model$res[2,1])
                     
                     le.men <- list()
                     for(i in 1:length(data_list_men)){
                       print(i)
                       le.men[[i]] <- find_le(model=model,
                                              data=data_list_men[[i]],
                                              sex="men",
                                              name="MDCIDCI",
                                              scb_hazard=scb_hazard_men,
                                              plotit= i==1)
                     }
                     rm(model)
                     
                     print("CCI men")
                     print(Sys.time())
                     modelname <- paste0(a,"_","CCI",
                                         "",
                                         maxfu,
                                         "",
                                         "man")
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                     temp <-  parameter.values.smooth[["CCI"]][["man"]]
                     
                     model$coefficients <- temp[rownames(temp)==a,]
                     model$res[,1] <- temp[rownames(temp)==a,]
                     model$res[2,1]<-exp(  model$res[2,1])
                     
                     le.men.cci <- list()
                     for(i in 1:length(data_list_men)){
                       print(i)
                       le.men.cci[[i]] <- find_le(model=model,
                                                  data=data_list_men[[i]],
                                                  sex="men",
                                                  name="CCI",
                                                  scb_hazard=scb_hazard_men,
                                                  plotit= i==1)
                     }
                     rm(model)
                     rm(scb_hazard_men)
                    
                     
                     
                     print("CRUDE women")
                     print(Sys.time())
                     modelname <- paste0(a,"_","CRUDE",
                                         "",
                                         maxfu,
                                         "",
                                         "woman")
                     
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                     temp <-  parameter.values.smooth[[ "CRUDE" ]][["woman"]]
                     
                     model$coefficients <- temp[rownames(temp)==a,]
                     model$res[,1] <- temp[rownames(temp)==a,]
                     model$res[2,1]<-exp(  model$res[2,1])
                     
                     le.women.crude <- list()
                     for(i in 1:length(data_list_men)){
                       print(i)
                       le.women.crude[[i]] <- find_le(model=model,
                                                      data=data_list_women[[i]],
                                                      sex="women",
                                                      name="CRUDE",
                                                      scb_hazard=scb_hazard_women,
                                                      plotit= FALSE,
                                                      simple=TRUE)
                     }
                     rm(model)
                     
                     
                     
                     print("MDCIDCI women")
                     print(Sys.time())
                     modelname <- paste0(a,"_","MDCIDCI",
                                         "",
                                         maxfu,
                                         "",
                                         "woman")
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                     temp <-  parameter.values.smooth[[ "MDCIDCI" ]][["woman"]]
                     
                     model$coefficients <- temp[rownames(temp)==a,]
                     model$res[,1] <- temp[rownames(temp)==a,]
                     model$res[2,1]<-exp(  model$res[2,1])
                     
                     le.women <- list()
                     for(i in 1:length(data_list_men)){
                       print(i)
                       le.women[[i]] <- find_le(model=model,
                                                data=data_list_women[[i]],
                                                sex="women",
                                                name="MDCIDCI",
                                                scb_hazard=scb_hazard_women,
                                                plotit= i==1)
                     }
                     rm(model)
                     
                     print("CCI women")
                     print(Sys.time())
                     modelname <- paste0(a,"_","CCI",
                                         "",
                                         maxfu,
                                         "",
                                         "woman")
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                     temp <-  parameter.values.smooth[[ "CCI" ]][["woman"]]
                     
                     model$coefficients <- temp[rownames(temp)==a,]
                     model$res[,1] <- temp[rownames(temp)==a,]
                     model$res[2,1]<-exp(  model$res[2,1])
                     
                     le.women.cci <- list()
                     for(i in 1:length(data_list_men)){
                       print(i)
                       le.women.cci[[i]] <- find_le(model=model,
                                                    data=data_list_women[[i]],
                                                    sex="women",
                                                    name="CCI",
                                                    scb_hazard=scb_hazard_women,
                                                    plotit= i==1)
                     }
                     rm(model)
                     rm(scb_hazard_women)
                     
                     
                     
                     data.m <- do.call(rbind,data_list_men)
                     rm(data_list_men)
                     
                     data.m$le <- unlist(le.men)
                     data.m$le.cci <- unlist(le.men.cci)
                     data.m$le.crude <- unlist(le.men.crude)
                     
                     rm(le.men)
                     rm(le.men.cci)
                     rm(le.men.crude)
                     
                     data.f <- do.call(rbind,data_list_women)
                     rm(data_list_women)
                     
                     data.f$le <- unlist(le.women)
                     data.f$le.cci <- unlist(le.women.cci)
                     data.f$le.crude <- unlist(le.women.crude)
                     
                     rm(le.women)
                     rm(le.women.cci)
                     rm(le.women.crude)
                     
                     LE <- rbind(data.m,data.f)
                     rm(data.m)
                     rm(data.f)
                     
                     print("Saving")
                     print(Sys.time())
                     save(LE,file=paste0("C:\\Marcus\\LEcox\\LE_age",a,".Rdata"))
                     
                     rm(list=ls())
                     gc()
                     
                     print("Done!")
                     print(Sys.time())
                     sink()
                     
                   })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)




















