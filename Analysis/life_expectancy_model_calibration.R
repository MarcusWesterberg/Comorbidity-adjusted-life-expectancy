###
# Survival models 

X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
require(tidyverse)

log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"
background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")
figs_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")


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

clusterExport(cl, c("X","data_path","figs_path","cluster_ids"),envir = environment())

calibration_list <- parLapplyLB(cl = cl,
                                cluster_indices,
                                fun = function(a){
                                  sink(paste0("C:\\Marcus\\Misc\\","LEcox_age_calibration",a,".txt"))
                                  require(dplyr)
                                  X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
                                  
                                  parameter.values.smooth.crude <- get(load(file=paste0(figs_path,"parameter.values.crude.smooth.Rdata")))
                                  parameter.values.smooth.cci <- get(load(file=paste0(figs_path,"parameter.values.cci.smooth.Rdata")))
                                  parameter.values.smooth.mdcidci <- get(load(file=paste0(figs_path,"parameter.values.mdcidci.smooth.Rdata")))
                                  parameter.values.smooth.ccidci <- get(load(file=paste0(figs_path,"parameter.values.ccidci.smooth.Rdata")))
                                  
                                  names(parameter.values.smooth.crude) <- c("man","woman")
                                  names(parameter.values.smooth.cci) <- c("man","woman")
                                  names(parameter.values.smooth.mdcidci) <- c("man","woman")
                                  names(parameter.values.smooth.ccidci) <- c("man","woman")
                                  parameter.values.smooth <- list("CRUDE"=parameter.values.smooth,
                                                                  "CCI"=parameter.values.smooth.cci,
                                                                  "MDCIDCI"=parameter.values.smooth.mdcidci,
                                                                  "CCIDCI"=parameter.values.smooth.ccidci)
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
                                  
                                  MDCI_cutoffs <- c(-1,4)
                                  DCI_cutoffs <- c(-0.5,8)
                                  
                                  
                                  data$indexage <- data$indexage - a - 0.5
                                  
                                  data$CCI[is.na(data$CCI)]<-0
                                  
                                  data$CCI[data$CCI>6] <- 6
                                  
                                  
                                  print(table(data$CCI,data$sex))
                                  
                                  trunc_it <- function(x,l,u){
                                    x[x<l]<-l
                                    x[x>u]<-u
                                    return(x)
                                  }
                                  
                                  data$MDCI <- trunc_it(data$MDCI,
                                                        l=MDCI_cutoffs[1],
                                                        u=MDCI_cutoffs[2])
                                  
                                  data$dci <- trunc_it(data$dci,
                                                       l=DCI_cutoffs[1],
                                                       u=DCI_cutoffs[2])
                                  
                                  require(survival)
                                  require(pec)
                                  require(flexsurv)
                                  require(splines)
                                  
                                  calibration_list <- list()
                                  current_inits <- NULL
                                  maxfu <- 10
                                  data$censor[data$timefu>maxfu] <- 0
                                  data$timefu[data$timefu>maxfu] <- maxfu
                                  
                                  
                                  models <- c("CCI",
                                              "MDCIDCI",
                                              "CCIDCI")
                                  
                                  
                                  for(sexstrata in c("man","woman")){
                                    print(sexstrata)
                                    
                                    data.s <- data %>%
                                      filter(sex %in% sexstrata)
                                    
                                    for(j in 1:length(models)){
                                      print(models[j] )
                                      
                                      modelname <- paste0(a,"_",models[j],
                                                          "",
                                                          maxfu,
                                                          "",
                                                          sexstrata)
                                      
                                      loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                                      temp <-  parameter.values.smooth[[ models[j] ]][[sexstrata]]
                                      
                                      model$coefficients <- temp[rownames(temp)==a,]
                                      model$res[,1] <- temp[rownames(temp)==a,]
                                      model$res[2,1]<-exp(  model$res[2,1])
                                      
                                      p <- predict(model,
                                                   newdata=data.s,
                                                   type="survival",
                                                   times=maxfu)
                                      data.s$p <- p$.pred_survival
                                      
                                      q <- unique(c(min(data.s$p),
                                                    quantile(data.s$p,probs=seq(0,1,length.out=10)) ))
                                      
                                      data.s$pcut <- cut(data.s$p,breaks=q,include.lowest = TRUE,right=TRUE)
                                      table(data.s$pcut)
                                      
                                      extract_obs <- function(timefu,censor){
                                        
                                        sf <- survfit(Surv(timefu, censor) ~ 1)
                                        
                                        sf$surv[sf$time== max(sf$time)]
                                      }
                                      
                                      cps <- data.s %>%
                                        group_by(pcut) %>%
                                        mutate(n=n(),
                                               p=mean(p),
                                               obs=extract_obs(timefu=timefu,censor=censor)) %>%
                                        select(sex,n,pcut,p,obs) %>%
                                        slice(1) %>%
                                        ungroup() %>%
                                        mutate(age=a,
                                               model_name=modelname)
                                      
                                      calibration_list[[modelname]] <- cps
                                      
                                    }
                                    
                                    cat("\n")
                                  }
                                  
                                  gc()
                                  
                                  return(do.call(rbind,calibration_list))
                                })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)


calibration_list <- do.call(rbind,calibration_list) %>%
  arrange(sex,age,model_name)
table(grepl(calibration_list$model_name,pattern="CCIDCI"))

save(calibration_list,file=paste0(figs_path,"calibration_list.Rdata"))

load(file=paste0(figs_path,"calibration_list.Rdata"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Calibration





###
# Summarize calibration plots 
ages <- c(20,30,40,50,60,70,80,90,100)
alpha <- 0.85
cols <- colorRampPalette(c(rgb(0,0.3,alpha),
                           rgb(0,0.8,alpha),
                           rgb(0.8,0.8,0,alpha),
                           rgb(0.8,0,0,alpha),
                           rgb(0.5,0,0,alpha),
                           rgb(0,0,0,alpha)) )(9)
pch <- c(16,17,16,17,
         16,17,16,17,
         16)

maxfu <- 10

cex <- 0.5

for(name in c("MDCIDCI","CCI","CCIDCI")){
  
  svg(paste0("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\calibrationplot",name,"_",maxfu,".svg"),width=10,height=6)
  par(mfrow=c(1,2),mar=c(3,3,4,1))
  
  for(j in 1:length(ages)){
    a <- ages[j]
    
    cps <- calibration_list %>%
      filter(sex=="man" & age %in% a & grepl(model_name,pattern=name) )
    
    if(j==1){
      plot(x=cps$p,
           y=cps$obs,
           xlab="",
           ylab="",
           main="Men",
           xlim=c(0,1),
           ylim=c(0,1),
           col=cols[j],
           pch=pch[j],
           cex=cex)
      lines(x=c(0,1),y=c(0,1),col="grey")
    } else{
      points(x=cps$p,
             y=cps$obs,
             col=cols[j],
             pch=pch[j],
             cex=cex)
    }
  }
  
  text(x=rep(0.1,9),
       y=seq(1,0.7,length.out=9),,labels=paste0("Age ",ages ))
  points(x=rep(0,9),
         y=seq(1,0.7,length.out=9),
         col=cols,
         pch=pch)
  mtext(side=1,text="Predicted probability",line=2)
  mtext(side=2,text="Observed probability",line=2)
  

  for(j in 1:length(ages)){
    a <- ages[j]
    
    cps <- calibration_list %>%
      filter(sex=="woman" & age %in% a & grepl(model_name,pattern=name) )
    
    if(j==1){
      plot(x=cps$p,
           y=cps$obs,
           xlab="",
           ylab="",
           main="Women",
           xlim=c(0,1),
           ylim=c(0,1),
           col=cols[j],
           pch=pch[j],
           cex=cex)
      lines(x=c(0,1),y=c(0,1),col="grey")
    } else{
      points(x=cps$p,
             y=cps$obs,
             col=cols[j],
             pch=pch[j],
             cex=cex)
    }
  }
  mtext(side=1,text="Prepicted probability",line=2)
  mtext(side=2,text="Observed probability",line=2)
  
  dev.off()
  
}

alpha <- 0.85
ages <- 18:100
cols <- colorRampPalette(c(rgb(0,0.3,alpha),
                           rgb(0,0.8,alpha),
                           rgb(0.8,0.8,0,alpha),
                           rgb(0.8,0,0,alpha),
                           rgb(0.5,0,0,alpha),
                           rgb(0,0,0,alpha)) )(length(ages))

pch <- rep(16,length(ages))
for(name in c("MDCIDCI","CCI","CCIDCI")){
  svg(paste0("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\calibrationplotALL",name,"_",maxfu,".svg"),width=10,height=6)
  par(mfrow=c(1,2),mar=c(3,3,4,1))
  
  for(j in 1:length(ages)){
    a <- ages[j]
    
    cps <- calibration_list %>%
      filter(sex=="man" & age %in% a & grepl(model_name,pattern=name) )
    if(j==1){
      plot(x=cps$p,
           y=cps$obs,
           xlab="",
           ylab="",
           main="Men",
           xlim=c(0,1),
           ylim=c(0,1),
           col=cols[j],
           pch=pch[j],
           cex=cex)
      lines(x=c(0,1),y=c(0,1),col="grey")
    } else{
      points(x=cps$p,
             y=cps$obs,
             col=cols[j],
             pch=pch[j],
             cex=cex)
    }
  }
  
  ages2 <- ages %in% c(18,30,40,50,60,70,80,90,100)
  text(x=rep(0.1,9),
       y=seq(1,0.7,length.out=9),
       labels=paste0("Age ",ages[ages2] ))
  points(x=rep(0,9),
         y=seq(1,0.7,length.out=9),
         col=cols[ages2],
         pch=pch[ages2])
  mtext(side=1,text="Predicted probability",line=2)
  mtext(side=2,text="Observed probability",line=2)
  
  
  for(j in 1:length(ages)){
    a <- ages[j]
    
    cps <- calibration_list %>%
      filter(sex=="woman" & age %in% a & grepl(model_name,pattern=name) )
    if(j==1){
      plot(x=cps$p,
           y=cps$obs,
           xlab="",
           ylab="",
           main="Women",
           xlim=c(0,1),
           ylim=c(0,1),
           col=cols[j],
           pch=pch[j],
           cex=cex)
      lines(x=c(0,1),y=c(0,1),col="grey")
    } else{
      points(x=cps$p,
             y=cps$obs,
             col=cols[j],
             pch=pch[j],
             cex=cex)
    }
  }
  mtext(side=1,text="Predicted probability",line=2)
  mtext(side=2,text="Observed probability",line=2)
  
  dev.off()
  
}


