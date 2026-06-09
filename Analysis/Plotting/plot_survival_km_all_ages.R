###
# Plot survival curves (km)


###
# Settings

X <- "\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)
require(parallel)

save_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")

data_path <- "C:\\Marcus\\ComorbidityBase_cache\\"

setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\"))

results_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")

log_path<-"C:\\Marcus\\Misc\\"

source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\function_create_subgroups.R"))
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
    # ccats <- c(
    #   "Myocardial_infarction","Congestive_heart_failure","Peripheral_vascular_disease","Cerebrovascular_disease","Chronic_obstructive_pulmonary_disease",
    #   "Chronic_other_pulmonary_disease","Rheumatic_disease","Dementia","Hemiplegia","Diabetes_without_chronic_complication",
    #   "Diabetes_with_chronic_complication","Renal_disease","Mild_liver_disease" ,"Liver_special", "Severe_liver_disease",                 
    #   "Peptic_ulcer_disease","Malignancy","Metastatic_solid_tumor","Aids"
    # )
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
                 "woman")


age_groups <- as.list(c(30:90))
age_groups <- append(list(c(18,29)),age_groups)
age_groups <- append(age_groups,list(c(91,120)))

subgroup <- list(c("Development","Validation"))

years_groups <- list(c(2006:2022))



subgroups<-create_subgroups(sex_groups=sex_groups,
                            age_groups=age_groups,
                            years=years_groups,
                            subgroup=subgroup,
                            subgroup2=NULL,
                            charlson=NULL,
                            icd_n_chapters=NULL)
length(subgroups$subgroups)
dim(subgroups$subgroups_matrix)
unique(subgroups$subgroups_matrix$icd_n_chapters )

subgroups_matrix <-subgroups$subgroups_matrix
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
subgroups <-subgroups$subgroups
length(subgroups)


length(subgroups$subgroups)



subgroup<-list(c("Development","Validation"))

comorb_colnames<-list(c("LopNr","date","dci"),
                      c("LopNr","date","MDCI"))

risk_score_name<-c("DCI","MDCI")


save_parts<-FALSE


mutate_risk_scores <- list("DCI"=function(x){
  x <- exp(x)
  x <- 1*(x<1.07) + 2*(x>=1.07 & x <1.84) + 3*(x>=1.84 & x< 3.38) + 4*(x>=3.38 & x < 10.79)  + 5*(x >= 10.79) 
  return(x)
},
"MDCI"=function(x){
  x <- exp(x)
  x <- 1*(x<0.88) + 2*(x >= 0.88 & x< 1) + 3*(x >= 1 & x<1.18) + 4*(x>= 1.18 & x<1.76)  +  5*(x>=1.76) 
  return(x)
})


###
# parallel: load data and 
#   estimate survival function in each strata
#   according to comorbidity indices subgroups
#   also estimate HRs 

###
# Load data...

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
      #print(paste0("at_times i=",i))
      if(index %% 10 == 0 ){
        cat("\n")
        print(Sys.time())
        cat("\n")
      }
      #if(index %% 5 == 0 ){
      print(paste0("Progress: ",index," out of ",total_iter,". Nrows=",nrow(data_temp)))
      
      helper <- function(data,
                         at_time,
                         weight=NULL,
                         boot_resample=FALSE){
        
        require(survival)
        
        ret <- list()
        
        
        require(survival)
        #require(SurvRegCensCov)
        #require("cmprsk")
        
        
        # competing risk
        
        
        start <- Sys.time()
        
        if(!is.null(weight)){
          weight <- data[[weight]]
        }
        
        if(nrow(data)>0 & length(unique(data$comorbidity_index))>1){
          
          if(boot_resample){
            stop("Bootstrap not supported")
            # boot_w <- resample_individual(data,
            #                               N=N)
            # 
            # data <- data %>%
            #   left_join(boot_w,by="LopNr") %>%
            #   filter(n_boot_repl>0)  
            # 
            # weight <- data$n_boot_repl
          }
          
          # prepare for evaluation at at_time
          data$censor[data$timefu>at_time]<-0
          data$timefu[data$timefu>at_time]<-at_time
          
          # n_un <- unique(data$comorbidity_index)
          # if( all.equal(n_un, as.integer(n_un)) ){
          #   few <- which(table(data$comorbidity_index)<100)
          #   if(length(few)>0){
          #     first <-   sort(as.numeric(names(few)))[1]
          #     data$comorbidity_index[ data$comorbidity_index>first] <- first
          #   }
          # 
          #   data$comorbidity_index <- factor(data$comorbidity_index)
          # }
          
          # 
          
          # sf <- survfit(formula=Surv(timefu, censor) ~ comorbidity_index ,
          #               data=data,
          #               id=data$LopNr,
          #               weight=weight)
          
          
          sf <- tryCatch(
            {
              # survreg(Surv(timefu, censor) ~ comorbidity_index,
              #         data = data,
              #         weights=weight,
              #         cluster=data$LopNr,
              #         dist="weibull",
              #         control=survreg.control(maxiter=100),
              #         init=c(2.5,-1,1))
              
              # sf <- survfit(formula=Surv(timefu, censor) ~ comorbidity_index ,
              #               data=data,
              #               id=data$LopNr,
              #               weight=weight)
              # 
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
          # survreg scale parameter maps to 1/shape, linear predictor to log(scale)
          # The Weibull distribution with shape parameter alpha and scale parameter sigma  has densit
          
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
            
            #surv_curves <- timepoints(cret,times=time_points)
            
            
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
      
      
      #}
      #print(paste0("subgroups j=",j))
      # data_temp <- data_temp %>% arrange(-timefu)
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
  #ret_inner <- do.call(rbind,ret_inner)
  
  return(ret_inner)
}



at_times <- c(10)

ncores <-  32
nbatch <- 64/ncores

start_time <- Sys.time()
cluster_indices <- 1:64
print(cluster_indices)
ncores <- min(ncores,length(cluster_indices))



save_name <- "AllagesKM"



print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

# subgroups <- subgroups[1:4]

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
                     
                   
                     
                     sink(file=paste0(log_path,"surv_km2_log_",save_name,"_",cluster_id,".txt"))
                     mutate_risk_scores <<- mutate_risk_scores
                     print(mutate_risk_scores)
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
                     
                     
                     #loads <- load(file=paste0(data_path,"charlson10_cache_part",cluster_id,".Rdata"))
                     #cci_data <- data
                     
                     loads <- load(file=paste0(data_path,"cindex_cache_part",cluster_id,".Rdata"))
                     
                     # data <- data %>%
                     #   left_join(cci_data %>% 
                     #               select(-any_of(c("rownr","CCIw","CCIunw"))),
                     #             by=c("LopNr","indexdate"))
                     # rm(cci_data)
                     # gc()
                     
                     ###
                     # keep only needed coluns
                     data <- data %>%
                       select(LopNr,rowid,indexdate,indexyear,subset,subset2,timefu,censor,sex,indexage,starts_with("comorbidity_index"))
                     
                     data$indexage <- floor(data$indexage)
                     data$indexage[data$indexage<29]<-29
                     data$indexage[data$indexage>91]<-91

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
                     
                     #ret <- do.call(rbind,ret)
                     
                     
                     if(save_parts){
                       print("Saving...")
                       save(ret,file=paste0(save_path,"surv_km_all_ages_inner_",save_name,"_part",cluster_id,".Rdata"))
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
save(all_survs,file=paste0(save_path,"all_survs_km_all_ages.Rdata"))





load(file=paste0(save_path,"all_survs_km_all_ages.Rdata"))

length(all_survs)
point_ests <- all_survs[[1]]

atindex <- 1

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
    
    point_ests[[k]][[l]][[atindex]]$n_ok <- 0
    
    if(length(point_ests[[k]][[l]][[atindex]]$est)>0){
      point_ests[[k]][[l]][[atindex]]$est <-  point_ests[[k]][[l]][[atindex]]$est*0
      point_ests[[k]][[l]][[atindex]]$var <-  point_ests[[k]][[l]][[atindex]]$var*0
      point_ests[[k]][[l]][[atindex]]$n_at_risk <- point_ests[[k]][[l]][[atindex]]$n_at_risk*0
      
      N <- point_ests[[k]][[l]][[atindex]]$N*0
      
      # For each simulation
      for(m in 1:length(all_survs)){
        
        if(!is.null(all_survs[[m]][[k]][[l]][[atindex]]$n_at_risk)){
          point_ests[[k]][[l]][[atindex]]$n_ok <-  point_ests[[k]][[l]][[atindex]]$n_ok + 1
          
          point_ests[[k]][[l]][[atindex]]$n_at_risk <- point_ests[[k]][[l]][[atindex]]$n_at_risk + all_survs[[m]][[k]][[l]][[atindex]]$n_at_risk
          
          N <- N + all_survs[[m]][[k]][[l]][[atindex]]$N
          
          temp <- all_survs[[m]][[k]][[l]][[atindex]]$est
          
          for(i in 1:ncol(temp)){
            temp[,i] <- temp[,i]*all_survs[[m]][[k]][[l]][[atindex]]$N
          }
          
          point_ests[[k]][[l]][[atindex]]$est <-  point_ests[[k]][[l]][[atindex]]$est + temp
          
          temp <- all_survs[[m]][[k]][[l]][[atindex]]$var
          for(i in 1:ncol(temp)){
            temp[,i] <- temp[,i]*(all_survs[[m]][[k]][[l]][[atindex]]$N)^2
          }
          
          point_ests[[k]][[l]][[atindex]]$var <-  point_ests[[k]][[l]][[atindex]]$var + temp
        }
        
        # end loop m
      }
      
      # divide / average for each time point i 
      for(i in 1:ncol(point_ests[[k]][[l]][[atindex]]$est)){
        point_ests[[k]][[l]][[atindex]]$est[,i] <- point_ests[[k]][[l]][[atindex]]$est[,i]/N
        point_ests[[k]][[l]][[atindex]]$var[,i] <- point_ests[[k]][[l]][[atindex]]$var[,i]/(N^2)
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
    point_ests[[k]][[l]] <-  point_ests[[k]][[l]][[atindex]]
  }
  #point_ests[[k]] <- point_ests[[k]][1:(length(point_ests[[k]])-2) ]
}

### 
# Check N ok
n_ok <- lapply(FUN=function(x){
  do.call(rbind,lapply(FUN=function(y){
    data.frame(subgroup_index=y$subgroup_index,n_ok=y$n_ok)
  },x))
},point_ests)

n_ok



save(point_ests,file=paste0(save_path,"point_ests_km_all_ages.Rdata"))





load(file=paste0(save_path,"point_ests_km_all_ages.Rdata"))










ags <- unique(subgroups_matrix$indexage)


cci_labs <- 1:5

cols <- colorRampPalette(c(rgb(0,0.75,0,0.75),rgb(0.75,0,0,0.75)))(5)

sexes <- c("man","woman")

xlim <- c(0,1)
ylim <- c(29,91)

mar <- c(2,2,2,1)
oma <- c(0.1,0.1,0.1,0.1)

width <- 5
height <- 4

ps_men <- point_ests$DCI[1:63]
ps_women <- point_ests$DCI[63+1:63]

svg(file=paste0(results_path,"Surv_DCI_ages_men.svg"),width=width,height=height)
par(mar=mar,oma=oma)

plot(1,xlim=ylim,ylim=xlim,type="n",axes=FALSE,xlab="",ylab="")

xat <- seq(xlim[1],xlim[2],by=0.1)
axis(side=2,at=xat,labels=paste0(xat*100,"%"),padj=0.5,cex.axis=0.75)

yat <- c(29,seq(35,85,by=5),91)
yat_labs <-yat
yat_labs[1]<-"<30"
yat_labs[length(yat_labs)]<-">90"

#text(y=0,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)
axis(side=1,at=yat,labels=rep("",length(yat)),cex.axis=0.75,line=-0.5)
#text(y=yat2,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)
text(y=-0.03,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)


for(i in 2:length(xat)){
  lines(y=c(xat[i],xat[i]),
        x=ylim + c(-0.5,0.5),
        col="grey")
  
} 

for( a in 1:63){
  temp <- ps_men[[a]]
  est_temp <- temp$est
  var_temp <- temp$var
  
  cil <- est_temp - sqrt(var_temp)*qnorm(0.975)
  ciu <- est_temp + sqrt(var_temp)*qnorm(0.975)
  
  x <- temp$time_points
  x_sel <- which(x %in% 10)
  
  for( cci in 1:5){
    
    points(y=est_temp[cci,x_sel],
           x=a+28,
           col=cols[cci],
           pch=18,
           cex=0.75)
    
    lines(x=as.numeric(c(a+28,a+28)),
          y=c(  cil[cci,x_sel] ,ciu[cci,x_sel]),
          col=cols[cci])
    
  }
}

dev.off()
  




svg(file=paste0(results_path,"Surv_DCI_ages_women.svg"),width=width,height=height)
par(mar=mar,oma=oma)

plot(1,xlim=ylim,ylim=xlim,type="n",axes=FALSE,xlab="",ylab="")

xat <- seq(xlim[1],xlim[2],by=0.1)
axis(side=2,at=xat,labels=paste0(xat*100,"%"),padj=0.5,cex.axis=0.75)

yat <- c(29,seq(35,85,by=5),91)
yat_labs <-yat
yat_labs[1]<-"<30"
yat_labs[length(yat_labs)]<-">90"

#text(y=0,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)
axis(side=1,at=yat,labels=rep("",length(yat)),cex.axis=0.75,line=-0.5)
#text(y=yat2,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)
text(y=-0.03,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)


for(i in 2:length(xat)){
  lines(y=c(xat[i],xat[i]),
        x=ylim + c(-0.5,0.5),
        col="grey")
  
} 

for( a in 1:63){
  temp <- ps_women[[a]]
  est_temp <- temp$est
  var_temp <- temp$var
  
  cil <- est_temp - sqrt(var_temp)*qnorm(0.975)
  ciu <- est_temp + sqrt(var_temp)*qnorm(0.975)
  
  x <- temp$time_points
  x_sel <- which(x %in% 10)
  
  for( cci in 1:5){
    
    points(y=est_temp[cci,x_sel],
           x=a+28,
           col=cols[cci],
           pch=18,
           cex=0.75)
    
    lines(x=as.numeric(c(a+28,a+28)),
          y=c(  cil[cci,x_sel] ,ciu[cci,x_sel]),
          col=cols[cci])
    
  }
}

dev.off()





ps_men <- point_ests$MDCI[1:63]
ps_women <- point_ests$MDCI[63+1:63]

svg(file=paste0(results_path,"Surv_MDCI_ages_men.svg"),width=width,height=height)
par(mar=mar,oma=oma)

plot(1,xlim=ylim,ylim=xlim,type="n",axes=FALSE,xlab="",ylab="")

xat <- seq(xlim[1],xlim[2],by=0.1)
axis(side=2,at=xat,labels=paste0(xat*100,"%"),padj=0.5,cex.axis=0.75)

yat <- c(29,seq(35,85,by=5),91)
yat_labs <-yat
yat_labs[1]<-"<30"
yat_labs[length(yat_labs)]<-">90"

#text(y=0,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)
axis(side=1,at=yat,labels=rep("",length(yat)),cex.axis=0.75,line=-0.5)
#text(y=yat2,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)
text(y=-0.03,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)


for(i in 2:length(xat)){
  lines(y=c(xat[i],xat[i]),
        x=ylim + c(-0.5,0.5),
        col="grey")
  
} 

for( a in 1:63){
  temp <- ps_men[[a]]
  est_temp <- temp$est
  var_temp <- temp$var
  
  cil <- est_temp - sqrt(var_temp)*qnorm(0.975)
  ciu <- est_temp + sqrt(var_temp)*qnorm(0.975)
  
  x <- temp$time_points
  x_sel <- which(x %in% 10)
  
  for( cci in 1:5){
    
    points(y=est_temp[cci,x_sel],
           x=a+28,
           col=cols[cci],
           pch=18,
           cex=0.75)
    
    lines(x=as.numeric(c(a+28,a+28)),
          y=c(  cil[cci,x_sel] ,ciu[cci,x_sel]),
          col=cols[cci])
    
  }
}

dev.off()





svg(file=paste0(results_path,"Surv_MDCI_ages_women.svg"),width=width,height=height)
par(mar=mar,oma=oma)

plot(1,xlim=ylim,ylim=xlim,type="n",axes=FALSE,xlab="",ylab="")

xat <- seq(xlim[1],xlim[2],by=0.1)
axis(side=2,at=xat,labels=paste0(xat*100,"%"),padj=0.5,cex.axis=0.75)

yat <- c(29,seq(35,85,by=5),91)
yat_labs <-yat
yat_labs[1]<-"<30"
yat_labs[length(yat_labs)]<-">90"

#text(y=0,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)
axis(side=1,at=yat,labels=rep("",length(yat)),cex.axis=0.75,line=-0.5)
#text(y=yat2,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)
text(y=-0.03,x=yat,labels=yat_labs,pos=1,cex=0.75,xpd=NA)


for(i in 2:length(xat)){
  lines(y=c(xat[i],xat[i]),
        x=ylim + c(-0.5,0.5),
        col="grey")
  
} 

for( a in 1:63){
  temp <- ps_women[[a]]
  est_temp <- temp$est
  var_temp <- temp$var
  
  cil <- est_temp - sqrt(var_temp)*qnorm(0.975)
  ciu <- est_temp + sqrt(var_temp)*qnorm(0.975)
  
  x <- temp$time_points
  x_sel <- which(x %in% 10)
  
  for( cci in 1:5){
    
    points(y=est_temp[cci,x_sel],
           x=a+28,
           col=cols[cci],
           pch=18,
           cex=0.75)
    
    lines(x=as.numeric(c(a+28,a+28)),
          y=c(  cil[cci,x_sel] ,ciu[cci,x_sel]),
          col=cols[cci])
    
  }
}

dev.off()








###
# End
###