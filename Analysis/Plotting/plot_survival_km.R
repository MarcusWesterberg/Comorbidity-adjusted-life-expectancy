###
# Plot survival curves (km)

###
# Settings

X <- "\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

save_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")

data_path <- "C:\\Marcus\\ComorbidityBase_cache\\"

setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\"))

results_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")

log_path<-"C:\\Marcus\\Misc\\"

source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\function_create_subgroups.R"))


###
# Subgroups
subgroup_function<-function(data,
                            subset_g,
                            subset2_g,
                            indexyear_g,
                            indexage_g,
                            sex_g,
                            charlson_cat_g,
                            charlson_g,
                            icd_n_chapters_g){
  
  if(!is.null(subset_g)){
    data<-data %>% filter(subset %in% subset_g)
  }
  
  if(!is.null(subset2_g)){
    data<-data %>% filter(subset2 %in% subset2_g)
  }
  
  if(!is.null(indexyear_g)){
    data<-data %>% filter(indexyear %in% indexyear_g)
  }
  
  if(!is.null(indexage_g)){
    if(length(indexage_g)==1){
      data <- data %>% filter( floor(indexage) %in% indexage_g)
    } else{
      data <- data %>% filter(indexage >= indexage_g[1] & indexage <= indexage_g[2] )
    }
    
  }
  
  if(!is.null(sex_g)){
    data<-data %>% filter(sex %in% sex_g)
  }
  
  if(!is.null(charlson_g)){
    data<- data %>% filter(comorbidity_index_CCI10 %in% charlson_g)
  }
  
  if(!is.null(charlson_cat_g)){

    for(k in 1:length(charlson_cat_g)){
      data <- data[data[[ charlson_cat_g[k] ]]>0,]
    }
  }
  
  if(!is.null(icd_n_chapters_g)){
    data<- data %>% filter(icd_n_chapters %in% icd_n_chapters_g)
  }
  
  
  return(data)
}



sex_groups<-list("man",
                 "woman",
                 c("man","woman"))

age_groups<-list(c(18,120),
                 c(18,49),
                 c(50,69),
                 c(70,120))

subgroup<-list(c("Development","Validation"))

years_groups <- c(list(c(2006:2022)))


subgroups<-create_subgroups(sex_groups=sex_groups,
                            age_groups=age_groups,
                            years=years_groups,
                            subgroup=subgroup,
                            subgroup2=NULL,
                            charlson=list(0:20,0,1:2,3:20))


length(subgroups$subgroups)
dim(subgroups$subgroups_matrix)
unique(subgroups$subgroups_matrix$icd_n_chapters )

subgroup<-list(c("Development","Validation"))

comorb_colnames<-list(c("LopNr","date","dci"),
                      c("LopNr","date","MDCI"))

risk_score_name<-c("DCI","MDCI")

save_parts<-TRUE



mutate_risk_scores <- list("DCI"=function(x){
 
  x <- exp(x)
  x <- 1*(x<1) + 2*(x>=1 & x <1.59) + 3*(x>=1.59 & x< 2.94) + 4*(x>=2.94 & x < 8.64)  + 5*(x >= 8.64) 
  return(x)
},
"MDCI"=function(x){
  
  x <- exp(x)
  x <- 1*(x<0.88) + 2*(x >= 0.88 & x< 1) + 3*(x >= 1 & x<1.11) + 4*(x>= 1.11 & x<1.64)  +  5*(x>=1.64) 
  return(x)
})

at_times <- c(10)
subgroups_matrix <-subgroups$subgroups_matrix
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
subgroups <-subgroups$subgroups
length(subgroups)


###
# parallel: load data and 
#   estimate survival function in each strata
#   according to comorbidity indices subgroups
#   also estimate HRs 

###
# Load data...
require(parallel)
compute_surv_inner <- function(data,
                               subgroup_function,
                               subgroups,
                               at_times,
                               weight=NULL,
                               boot_resample=FALSE){
  
  ret_inner <- list()
  total_iter <- length(at_times)*length(subgroups)
  index <- 1
  
  for(j in 1:length(subgroups)){
    
    data_temp <- subgroup_function(data,
                                   subset_g=subgroups[[j]]$subgroup,
                                   subset2_g=subgroups[[j]]$subgroup2,
                                   indexyear_g=subgroups[[j]]$indexyear,
                                   indexage_g=subgroups[[j]]$indexage,
                                   sex_g=subgroups[[j]]$sex,
                                   charlson_cat_g=subgroups[[j]]$charlson_cat,
                                   charlson_g=subgroups[[j]]$charlson,
                                   icd_n_chapters_g=subgroups[[j]]$icd_n_chapters)  
    
    ret_temp <- list()
    for(i in 1:length(at_times)){
      if(index %% 10 == 0 ){
        cat("\n")
        print(Sys.time())
        cat("\n")
      }
      print(paste0("Progress: ",index," out of ",total_iter,". Nrows=",nrow(data_temp)))
      
      helper <- function(data,
                         at_time,
                         weight=NULL,
                         boot_resample=FALSE){
        
        require(survival)
        
        ret <- list()
        
        
     
        
        start <- Sys.time()
        
        if(!is.null(weight)){
          weight <- data[[weight]]
        }
        
        if(nrow(data)>0 & length(unique(data$comorbidity_index))>1){
          
          if(boot_resample){
            stop("Bootstrap not supported")
  
          }
          
          # prepare for evaluation at at_time
          data$censor[data$timefu>at_time]<-0
          data$timefu[data$timefu>at_time]<-at_time
          
    
          
          sf <- tryCatch(
            {

              survfit(formula=Surv(timefu, censor) ~ comorbidity_index ,
                      data=data,
                      id=data$LopNr,
                      weight=weight)
              
            },
            error = function(cond) {
              message(paste("Error in survreg."))
              message("Here's the original error message:")
              message(conditionMessage(cond))
              # Choose a return value in case of error
              NA
            }
          )
  
          
          if(class(sf) %in% c("survfit") ){
            time_points <- seq(0,at_time,by=0.01)
            
            sfs <- summary(sf)
            
            sf_indices <- list()
            ind_temp <- 1
            sfss <- cumsum(table(sfs$strata))
            for(l in 1:length(sfss)){
              sf_indices[[l]] <- ind_temp:sfss[l]
              ind_temp <- sfss[l]+1
            }
            
            est <- matrix(NA,nrow=length(sf$strata),ncol=length(time_points))
            var <- matrix(NA,nrow=length(sf$strata),ncol=length(time_points))
            
            rownames(est)<-names(sf$strata)
            rownames(var)<-names(sf$strata)
            
            colnames(est)<-time_points
            colnames(var)<-time_points
            
            
            
            for(l in 1:length(sf$strata)){
              est[l,] <- approx(x=sfs$time[ sf_indices[[l]] ],
                                y=sfs$surv[ sf_indices[[l]] ],
                                xout=time_points,
                                method="constant",
                                yleft=1,
                                rule=2)$y
              
              var[l,] <- approx(x=sfs$time[ sf_indices[[l]] ],
                                y=sfs$std.err[ sf_indices[[l]] ],
                                xout=time_points,
                                method="constant",
                                yleft=0,
                                rule=2)$y^2
            }
            
            surv_curves <- list("est"=1-est,
                                "var"=var)
      
            
            surv_curves$time_points <- time_points
            surv_curves$comorbidity_index_levels <- levels(data$comorbidity_index)
            surv_curves$N <- table(data$comorbidity_index)
            surv_curves$N_at_risk_time_points <- c(0,sort(unique(round(time_points))))
            
            n_at_risk <- matrix(0,
                                nrow=length(surv_curves$N_at_risk_time_points),
                                ncol=length(surv_curves$comorbidity_index_levels))
            
            n_at_risk[1,] <- surv_curves$N
            
            for(k_time in 2:length(surv_curves$N_at_risk_time_points)  ){
              for(k_lvls in 1:length(surv_curves$comorbidity_index_levels)){
                n_at_risk[k_time,k_lvls] <- sum(data$timefu[data$comorbidity_index %in% surv_curves$comorbidity_index_levels[k_lvls]] >= surv_curves$N_at_risk_time_points[k_time]  )
   
              }
            }
            rownames(n_at_risk) <- surv_curves$N_at_risk_time_points
            colnames(n_at_risk) <- surv_curves$comorbidity_index_levels
            
            surv_curves$n_at_risk <- n_at_risk
            
            ret <- surv_curves
            
          } else{
            print("Not a cuminc object - probably failed.")
          }
        }
        return(ret)
      }
      
      
      ret_temp[[i]] <- helper(data=data_temp,
                              at_time=at_times[i],
                              weight=weight,
                              boot_resample=boot_resample)
      
      ret_temp[[i]]$subgroup_index <- j
      ret_temp[[i]]$subgroup <- paste0(subgroups[[j]]$subgroup,collapse=",")
      ret_temp[[i]]$subgroup2 <- paste0(subgroups[[j]]$subgroup2,collapse=",")
      ret_temp[[i]]$indexyear <- paste0(subgroups[[j]]$indexyear,collapse=", ")
      ret_temp[[i]]$indexage <- paste0(subgroups[[j]]$indexage,collapse="-")
      ret_temp[[i]]$charlson <- paste0(subgroups[[j]]$charlson,collapse=",")
      ret_temp[[i]]$charlson_cat <- paste0(subgroups[[j]]$charlson_cat,collapse=",")
      ret_temp[[i]]$icd_n_chapters <- paste0(subgroups[[j]]$icd_n_chapters,collapse=",")
      
      index <- index + 1
    }
    ret_inner[[j]] <- ret_temp
  }
  
  return(ret_inner)
}





ncores <-  32
nbatch <- 64/ncores

start_time <- Sys.time()
cluster_indices <- 1:64
print(cluster_indices)
ncores <- min(ncores,length(cluster_indices))



save_name <- "All"



print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("save_name","X","log_path","at_times","data_path","subgroups","subgroup_function","save_parts","save_path",
                    "compute_surv_inner","risk_score_name","mutate_risk_scores"),envir = environment())

ret_list <- list()



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
                     start_time <- Sys.time()
                     cluster_id <<- cluster_id # make seed global
                     
                     
                     sink(file=paste0(log_path,"surv_km_log_",save_name,"_",cluster_id,".txt"))
                     
                     print(cluster_id)
                     print(Sys.time())
                     
                     cat("\n ... \n")
                     print(log_path)
                     print(at_times)
                     print(save_parts)
                     print(save_path)
                     print(save_name)
                     print(risk_score_name)
                     cat("\n ... \n")
                     
                     
                     
                     print("Loading data...")
                     Sys.sleep(cluster_id/16)
                     
                 
                     
                     loads <- load(file=paste0(data_path,"cindex_cache_part",cluster_id,".Rdata"))
             
                     ###
                     # keep only needed coluns
                     data <- data %>%
                       select(LopNr,rowid,indexdate,indexyear,subset,subset2,timefu,censor,sex,indexage,starts_with("comorbidity_index"))
                     
                     
          
                     ###
                     # Drop stuff not needed
                     all_subgroups <- unique(unlist(lapply(FUN=function(x) x$subgroup ,subgroups)  ))
                     if(!is.null(all_subgroups)){
                       data <- data %>%
                         filter(subset %in% all_subgroups)
                     }
                     
                     all_subgroups2 <- unique(unlist(lapply(FUN=function(x) x$subgroup2 ,subgroups)  ))
                     if(!is.null(all_subgroups2)){
                       data <- data %>%
                         filter(subset2 %in% all_subgroups2)
                     }
                     
                     all_ages <- unique(unlist(lapply(FUN=function(x) x$indexage ,subgroups)  ))
                     if(!is.null(all_ages)){
                       data <- data %>%
                         filter(indexage >= min(all_ages) & indexage <= max(all_ages))
                     }
                     
                     all_years <- unique(unlist(lapply(FUN=function(x) x$indexyear ,subgroups)  ))
                     if(!is.null(all_years)){
                       data <- data %>%
                         filter(indexyear %in% all_years)
                     }
                     ###
                     
                     
                     all_icd_chapt <- unique(unlist(lapply(FUN=function(x) x$icd_n_chapters ,subgroups)  ))
                     if(!is.null(all_icd_chapt)){
                       print("Loading ICD chapters")
                       icd_val <- get(load(file=paste0(data_path,"validation_ICD_lookback_10_part_",cluster_id,".Rdata")))
                       rm(cdata)
                       
                       icd_dev <- get(load(file=paste0(data_path,"development_ICD_lookback_10_part_",cluster_id,".Rdata")))
                       rm(cdata)
                       
                       icd <- rbind(icd_val,icd_dev) %>%
                         select(-rownr) %>%
                         rename(indexdate=date)
                       
                       rm(icd_val)
                       rm(icd_dev)
                       
                       
                       icd[["AB"]] <- icd[["A"]] | icd[["B"]]
                       # icd[["CD"]] <- icd[["C"]] | icd[["D"]]
                       icd[["ST"]] <- icd[["S"]] | icd[["T"]]
                       
                       # Compute n chapters
                       icd$icd_n_chapters <- 0
                       icd_chapters <- c("AB","CD_can","E","F","G","H","I","J","K","L","M", "N",  "ST")
                       
                       for(k in 1:length(icd_chapters)){
                         icd$icd_n_chapters <-   icd$icd_n_chapters + as.numeric(icd[[icd_chapters[k]]])
                       }
                       
                       data <- data %>%
                         left_join(icd %>% select(LopNr,indexdate,icd_n_chapters),by=c("LopNr","indexdate"))
                       data$icd_n_chapters[is.na(data$icd_n_chapters)] <- 0
                       
                       rm(icd)
                       
                       print("Filtering ICD chapters")
                       
                       data <- data %>%
                         filter(icd_n_chapters %in% all_icd_chapt)
                       
                     }
                     
                     
                     ret <- list()
                     
                     for(k in 1:length(risk_score_name)){
                       cat("\n")
                       print(paste0("Risk score ",risk_score_name[k]))
                       
                       ###
                       # Cindex
                       print( Sys.time() - start_time)
                       print("Computing compute_surv_inner")
                       data$comorbidity_index <- data[[paste0("comorbidity_index_",risk_score_name[k])]]
                       
                       if(risk_score_name[k] %in% "CCI10"){
                         data$comorbidity_index <- factor(mutate_risk_scores$CCI10(data$comorbidity_index))
                       } else if(risk_score_name[k] %in% "DCI"){
                         data$comorbidity_index <- factor(mutate_risk_scores$DCI(data$comorbidity_index))
                       } else if(risk_score_name[k] %in% "MDCI"){
                         data$comorbidity_index <- factor(mutate_risk_scores$MDCI(data$comorbidity_index))
                       }
                 
                       table(data$comorbidity_index,data$comorbidity_index_CCI10)
                       
                       ret_inner <- compute_surv_inner(data=data,
                                                       subgroups=subgroups,
                                                       subgroup_function=subgroup_function,
                                                       boot_resample=FALSE,
                                                       at_times=at_times)
                       
                       ret_inner$risk_score_name <- risk_score_name[k]
                       ret_inner$cluster_id <- cluster_id
                       
                       ret[[k]] <- ret_inner
                       
                     }
                     
                     rm(data)
                     gc()
                     
                     print("compute_surv_inner done")
                     
                     names(ret) <- risk_score_name
                  
                   
                     if(save_parts){
                       print("Saving...")
                       save(ret,file=paste0(save_path,"surv_km_inner_",save_name,"_part",cluster_id,".Rdata"))
                     }
                     
                     print("Done.")
                     sink()
                     
                     return(ret)
                     
                   })
  ret_list[[j]]<-ret
  rm(ret)
}

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)


all_survs <- unlist(ret_list,recursive=FALSE)
save(all_survs,file=paste0(save_path,"all_survs_km.Rdata"))

load(file=paste0(save_path,"all_survs_km.Rdata"))

length(all_survs)
point_ests <- all_survs[[1]]



# For each risk score
for(k in 1:length(point_ests)){
  cat("\n")
  print(paste0("k=",k))
  print(names(point_ests)[k])
  
  point_ests[[k]]<-point_ests[[k]][names(point_ests[[k]]) %in% ""] # drop labels
  
  # for each subgroup
  for(l in 1:nrow(subgroups_matrix)  ){
    print(paste0("l=",l))
    
    print(subgroups_matrix %>%
            filter(subgroup_index==l))
    
    point_ests[[k]][[l]][[1]]$n_ok <- 0
    
    if(length(point_ests[[k]][[l]][[1]]$est)>0){
      point_ests[[k]][[l]][[1]]$est <-  point_ests[[k]][[l]][[1]]$est*0
      point_ests[[k]][[l]][[1]]$var <-  point_ests[[k]][[l]][[1]]$var*0
      point_ests[[k]][[l]][[1]]$n_at_risk <- point_ests[[k]][[l]][[1]]$n_at_risk*0
      
      N <- point_ests[[k]][[l]][[1]]$N*0
      
      # For each simulation
      for(m in 1:length(all_survs)){
        
        if(!is.null(all_survs[[m]][[k]][[l]][[1]]$n_at_risk)){
          point_ests[[k]][[l]][[1]]$n_ok <-  point_ests[[k]][[l]][[1]]$n_ok + 1
          
          point_ests[[k]][[l]][[1]]$n_at_risk <- point_ests[[k]][[l]][[1]]$n_at_risk + all_survs[[m]][[k]][[l]][[1]]$n_at_risk
          
          N <- N + all_survs[[m]][[k]][[l]][[1]]$N
          
          temp <- all_survs[[m]][[k]][[l]][[1]]$est
          
          for(i in 1:ncol(temp)){
            temp[,i] <- temp[,i]*all_survs[[m]][[k]][[l]][[1]]$N
          }
          
          point_ests[[k]][[l]][[1]]$est <-  point_ests[[k]][[l]][[1]]$est + temp
          
          temp <- all_survs[[m]][[k]][[l]][[1]]$var
          for(i in 1:ncol(temp)){
            temp[,i] <- temp[,i]*(all_survs[[m]][[k]][[l]][[1]]$N)^2
          }
          
          point_ests[[k]][[l]][[1]]$var <-  point_ests[[k]][[l]][[1]]$var + temp
        }
        
        # end loop m
      }
      
      # divide / average for each time point i 
      for(i in 1:ncol(point_ests[[k]][[l]][[1]]$est)){
        point_ests[[k]][[l]][[1]]$est[,i] <- point_ests[[k]][[l]][[1]]$est[,i]/N
        point_ests[[k]][[l]][[1]]$var[,i] <- point_ests[[k]][[l]][[1]]$var[,i]/(N^2)
      }
    }
    
    # end loop l
  }
  
  # end loop k
}

###
# reformat
for(k in 1:length(point_ests)){
  print(k)
  
  for(l in 1: (length(point_ests[[k]])) ){
    point_ests[[k]][[l]] <-  point_ests[[k]][[l]][[1]]
  }
}


### 
# Check N ok
n_ok <- lapply(FUN=function(x){
  do.call(rbind,lapply(FUN=function(y){
    data.frame(subgroup_index=y$subgroup_index,n_ok=y$n_ok)
  },x))
},point_ests)

print(n_ok)




plot_survq <- function(sp,
                       xlim=c(0,10),
                       ylim=c(0,1),
                       xlab="Years since index date",
                       ylab="Cumulative probability of death",
                       add=FALSE,
                       col=1:10,
                       n_at_risk_at=c(0,5,10),
                       n_at_risk_y=-0.25,
                       n_at_risk_x=-3,
                       n_at_risk_div=100000,
                       n_at_risk_yscale=0.1,
                       est_coords=NULL,
                       est_time=NULL,
                       xl = c("0","1","2","3","\u22654"),
                       xl2=xl,
                       topleft_corner_y=1,
                       replace.dot=NULL){
  
  
  if(!add){
    plot(1,xlim=xlim,ylim=ylim,type="n",axes=FALSE,xlab="",ylab="")
    
    xat <- n_at_risk_at
    axis(side=1,at=xat,padj=-0.95,line=-0.4)
    
    yat <- seq(ylim[1],ylim[2],length.out=5)
    axis(side=2,at=yat,labels=paste0(yat*100,"%"),padj=0.95,line=-0.4)
    
    mtext(side=1,line=1.75,xlab,cex=0.75)
    mtext(side=2,line=1.75,ylab,cex=0.75)
    
    for(j in 2:length(xat)){
      lines(c(xat[j],xat[j]),ylim,col="lightgrey",lty=2)
    }
    for(j in 2:length(yat)){
      lines(xlim,c(yat[j],yat[j]),col="lightgrey",lty=2)
    }
    
    
    
    
    
    ###
    # add n at risk
    n_at_risk <- sp$n_at_risk
    fixdig <- function(x){
      
      helper <- function(x){
        if(x<0.1){ 
          x <- format(round(x,dig=2),digits=2)
        } else if(x<1){
          x <- format(round(x,dig=1),digits=1)
        } else(
          x <- round(x)
        )
        return(x)
      }
      return(apply(FUN=helper,x,MARGIN=c(1,2)))
    }
    
    n_at_risk <- fixdig(n_at_risk/n_at_risk_div)
    
    n_at_risk <- n_at_risk[rownames(n_at_risk) %in% as.character(n_at_risk_at), ]
    n_at_risk <- n_at_risk[-1,]
    
    text(x=n_at_risk_x,
         y=n_at_risk_y,
         labels="N at risk",
         xpd=NA)
    
    for(j in 1:ncol(n_at_risk)){
      
      text(x=-1,
           y=rep(n_at_risk_y,length(n_at_risk_at))-(j-1)*n_at_risk_yscale,
           labels=paste0(xl[j],""),
           xpd=NA,
           col=col[j],
           pos=2)
      
 
      ntemp <- n_at_risk[,j]
      if(!is.null(replace.dot)){
        ntemp <- str_replace(ntemp,pattern="[.]",replacement=replace.dot)
      }
      text(x=n_at_risk_at,
           y=rep(n_at_risk_y,length(n_at_risk_at))-(j-1)*n_at_risk_yscale,
           labels=ntemp,
           xpd=NA,
           col=col[j])
    }
    
    
    
  }

  
  for(j in 1:nrow(sp$est)){
    
    est_temp <- sp$est[j,]
    var_temp <- sp$var[j,]
    cil <- est_temp - sqrt(var_temp)*qnorm(0.975)
    ciu <- est_temp + sqrt(var_temp)*qnorm(0.975)
    
    x <- sp$time_points
    
    lines(x=x,y=est_temp,col=col[j],lwd=2)
    
  }
  
  
  ###
  # Add estimates 
  
  if(!is.null(est_time)){
    if(est_coords[1]<5){
      polygon(x=c(0.1,9,9,0.1),
              y=topleft_corner_y - c(0,0,0.45,0.45),
              border=NA,
              col=rgb(1,1,1,alpha=0.75))
    } else{
      polygon(x=c(2,10,10,2),
              y=topleft_corner_y - c(0,0,0.45,0.45),
              border=NA,
              col=rgb(1,1,1,alpha=0.75))
    }
  }
  
  
  for(j in 1:nrow(sp$est)){
    
    est_temp <- sp$est[j,]
    var_temp <- sp$var[j,]
    cil <- est_temp - sqrt(var_temp)*qnorm(0.975)
    ciu <- est_temp + sqrt(var_temp)*qnorm(0.975)
    
    x <- sp$time_points
    
    if(!is.null(est_time)){
      wt <- which(x %in% est_time)
      
      rd <- function(x){
        if(x<1){
          x <- format(round(x,dig=2),digits =1,nsmall =1)
        } else if (x<10){
          x <- format(round(x,dig=2),digits =1,nsmall =1)
        }else{
          x <- format(round(x,dig=2),digits =1,nsmall =0)
        }
        return(x)
      } 

      
      est_info <- rd(est_temp[wt]*100)
      cil_info <- rd(cil[wt]*100)
      ciu_info <- rd(ciu[wt]*100)
      
      if(!is.null(replace.dot)){
        est_info <- str_replace(est_info,pattern="[.]",replacement=replace.dot)
        cil_info <- str_replace(cil_info,pattern="[.]",replacement=replace.dot)
        ciu_info <- str_replace(ciu_info,pattern="[.]",replacement=replace.dot)
      }
      
      est_info <- gsub(est_info,pattern=" ",replacement="")
      cil_info <- gsub(cil_info,pattern=" ",replacement="")
      ciu_info <- gsub(ciu_info,pattern=" ",replacement="")
      
      text(x=4,
           y=est_coords[2]-(j-1)*0.065,
           labels=paste0(xl2[j],""),
           col=col[j],
           xpd=NA,
           pos=2)
      
      text(x=est_coords[1],
           y=est_coords[2]-(j-1)*0.065,
           labels=paste0(est_info),
           #col=col[j],
           col="black",
           xpd=NA,
           pos=4)
      
      if(j==1){
        text(x=est_coords[1],
             y=est_coords[2]-(-1)*0.065,
             labels="Est",
             col="black",
             xpd=NA,
             pos=4)
        text(x=est_coords[1]+1.59,
             y=est_coords[2]-(-1)*0.065,
             labels=paste0("(95% CI)"),
             col="black",
             xpd=NA,
             pos=4)
      }
      
      text(x=est_coords[1]+1.59,
           y=est_coords[2]-(j-1)*0.065,
           labels=paste0("(", cil_info ,"-",  ciu_info , ")"),
           #col=col[j],
           col="black",
           xpd=NA,
           pos=4)
      
    }
  }
}



# plot(1:4,col=col)



ags <- unique(subgroups_matrix$indexage)

iys <- unique(subgroups_matrix$indexyear)


ylims_age <- list(c(0,1),
                  c(0,1),
                  c(0,1),
                  c(0,1))
x <- 2.85
y <- 0.85
est_coords <- list(c(x,y),
                   c(x,y),
                   c(x,y),
                   c(x,y))
est_time <- NULL

cs <- 1
a <- 1

cci_labs <- c("All",0,"1-2","\u22653")

ccigs <- unique(  subgroups_matrix$charlson)

cols <- c( colorRampPalette(c(rgb(0,0.5,0.2),"wheat"))(4)[2])

cols <- c(cols,
           colorRampPalette(c(cols,"tan2","firebrick1"))(4)[-1],
           rep("red4",1))

sexes <- c( "man, woman","man","woman")


mfrow <- c(4,4)
mar <- c(4.8,0.95,0.01,0.01)
oma <- c(1,1,1.55,0.5)

width <- 6
height <- 7
n_at_risk_at <- c(0,5,10)
n_at_risk_x <- -1000

n_at_risk_yscale <- 0.0875
n_at_risk_y <- -0.22





#####################################################


mains <- c("All ages","18-49 years","50-69 years","\u226570 years")


for(sx in 1:3){
  svg(file=paste0(results_path,"Surv","_",sexes[sx],"_MDCI.svg"),width=width,height=height)
  par(mfrow=mfrow,mar=mar,oma=oma)
  for(cci in 1:length(ccigs)){  
    for(a in 1:length(ags)){
      
      sg <- subgroups_matrix$indexage %in% ags[a] & 
        subgroups_matrix$sex %in% sexes[sx] &
        subgroups_matrix$charlson %in% ccigs[cci]
     
      yl <- "Cumulative probability of death"

      
      
      xl <- c("","","","","")
      
      topleft_corner_y <- 1

      plot_survq(sp=point_ests[["MDCI"]][[which(sg) ]]  ,
                 xlim=c(0,10),
                 ylim=ylims_age[[a]],
                 xlab="",
                 ylab="",
                 add=FALSE,
                 col=cols,
                 n_at_risk_at=n_at_risk_at,
                 n_at_risk_x=n_at_risk_x,
                 n_at_risk_yscale=n_at_risk_yscale*diff(ylims_age[[a]]),
                 n_at_risk_y=n_at_risk_y*diff(ylims_age[[a]]),
                 est_time=est_time,
                 est_coords=est_coords[[a]],
                 xl=xl,
                 xl2=c("","","","",""),
                 topleft_corner_y=topleft_corner_y,
                 replace.dot="\U00B7")
      if(cci==1){
        mtext(side=3,line=0,mains[a])
      }
     
    }
  }
  dev.off()
}







for(sx in 1:3){
  svg(file=paste0(results_path,"Surv","_",sexes[sx],"_DCI.svg"),width=width,height=height)
  par(mfrow=mfrow,mar=mar,oma=oma)
  for(cci in 1:length(ccigs)){  
    for(a in 1:length(ags)){
      
      sg <- subgroups_matrix$indexage %in% ags[a] & 
        subgroups_matrix$sex %in% sexes[sx] &
        subgroups_matrix$charlson %in% ccigs[cci]
      
      yl <- "Cumulative probability of death"
      
      
      
      xl <- c("","","","","")
      
      topleft_corner_y <- 1
      plot_survq(sp=point_ests[["DCI"]][[which(sg) ]]  ,
                 xlim=c(0,10),
                 ylim=ylims_age[[a]],
                 xlab="",
                 ylab="",
                 add=FALSE,
                 col=cols,
                 n_at_risk_at=n_at_risk_at,
                 n_at_risk_x=n_at_risk_x,
                 n_at_risk_yscale=n_at_risk_yscale*diff(ylims_age[[a]]),
                 n_at_risk_y=n_at_risk_y*diff(ylims_age[[a]]),
                 est_time=est_time,
                 est_coords=est_coords[[a]],
                 xl=xl,
                 xl2=c("","","","",""),
                 topleft_corner_y=topleft_corner_y,
                 replace.dot="\U00B7")
      if(cci==1){
        mtext(side=3,line=0,mains[a])
      }
      
    }
  }
  dev.off()
}













###
# Summarize survival into table with confidence intervals
# temp <- point_ests$MDCI

at_times <- c(1,5,10)

dig <- 3

surv_ci <- point_ests
for(j in 1:length(surv_ci)){
  
  temp_list <- list()
  for(k in 1:length(surv_ci[[j]])){
    
    temp <-  surv_ci[[j]][[k]]
    
    ests <- temp$est[,temp$time_points %in% at_times]
    vars <- temp$var[,temp$time_points %in% at_times]
    
    if(!is.null(ests)){
      sds <- sqrt(vars)
      
      cil <- ests - sds*qnorm(0.975)
      ciu <- ests + sds*qnorm(0.975)
      
      fix_dig <- function(x,dig){
        
        inner <- function(x){
          ret <- NA
          if(!is.na(x)){
            if(x<=1){
              ret <-    format(round(x,dig=dig),digits=dig,nsmall=dig)
            } else if(x<=10){
              ret <-    format(round(x,dig=dig-1),digits=dig-1,nsmall=dig-1)
            } else{
              ret <-    format(round(x,dig=dig-2),digits=dig-2,nsmall=dig-2)
            }
          }
          return(ret)
        }
        sapply(FUN=inner,x)
      }
      
      temp2 <- data.frame(comorbidity_index_levels=rep(temp$comorbidity_index_levels,
                                                       each=length(at_times)),
                          timepoint=at_times)
      
      for(l in 1:length(temp$comorbidity_index_levels)){
        temp2$est[temp2$comorbidity_index_levels %in% temp$comorbidity_index_levels[l]] <- fix_dig(100*ests[l,],dig=dig)
        
        temp2$ci[temp2$comorbidity_index_levels %in% temp$comorbidity_index_levels[l]] <- paste0("(",fix_dig(100*cil[l,],dig=dig),
                                                                                                 "-",
                                                                                                 fix_dig(100*ciu[l,],dig=dig),")")
      }
      
      temp2$comorbidity_index <- names(point_ests)[j]
      temp2$subgroup_index <-    temp$subgroup_index
      temp2$subgroup2 <-    temp$subgroup2
      temp2$indexyear <-   temp$indexyear
      temp2$indexage <-    temp$indexage
      temp2$charlson <-  temp$charlson
      temp2$charlson_cat <-    temp$charlson_cat
      temp_list[[k]] <- temp2
    }
  }
  surv_ci[[j]] <- do.call(rbind,temp_list)
}


temp <- lapply(FUN=function(x) x %>% filter(timepoint==1),surv_ci)
openxlsx::write.xlsx(temp,file=paste0(results_path,"surv_km_1y.xlsx"))

temp <- lapply(FUN=function(x) x %>% filter(timepoint==5),surv_ci)
openxlsx::write.xlsx(temp,file=paste0(results_path,"surv_km_5y.xlsx"))

temp <- lapply(FUN=function(x) x %>% filter(timepoint==10),surv_ci)
openxlsx::write.xlsx(temp,file=paste0(results_path,"surv_km_10y.xlsx"))


###
# Relative comparison
temp_rel <- temp
temp_rel <- lapply(FUN=function(d){
  d <- d %>%
    select(est,comorbidity_index,comorbidity_index_levels,indexage,charlson) %>%
    group_by(indexage,charlson) %>%
    filter(comorbidity_index_levels %in% c("1","5")) %>%
    mutate(est=as.numeric(est)) %>%
    mutate(relative_risk = round(est[2]/est[1]),
           comorbidity_index_levels="5 / 1") %>% 
    slice(1) %>%
    ungroup() %>%
    arrange(charlson,indexage) %>%
    select(-est)
},temp_rel)

openxlsx::write.xlsx(temp_rel,file=paste0(results_path,"relative_risk_10y.xlsx"))







###
# End
###