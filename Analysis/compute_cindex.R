###
# Script that computes c indices

X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

ncores<-64
log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")


results_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")



source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\function_cindex.R"))

source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\function_create_subgroups.R"))

cluster_indices<-1:64




dig <- 3


fix_dig <- function(x,dig){
  
  inner <- function(x){
    ret <- NA
    if(!is.na(x)){
      ret <-    format(round(x,dig=dig),digits=dig,nsmall=dig)
    }
    return(ret)
  }
  sapply(FUN=inner,x)
}

fix_table <- function(cindices,
                      subgroup_matrix,
                      dig){
  
  c_ci <- cindices$point_ests_aggr %>%
    select(risk_score_name,timepoint,subgroup_index,cindex_mean,cindex_varpooled)
  
  c_ci$est <- c_ci$cindex_mean
  
  c_ci$sd <- sqrt(c_ci$cindex_varpooled)
  
  c_ci$cil <- fix_dig(c_ci$est - c_ci$sd*qnorm(0.975),dig=dig)
  c_ci$ciu <- fix_dig(c_ci$est + c_ci$sd*qnorm(0.975),dig=dig)
  c_ci$est <- fix_dig(c_ci$est,dig=dig)
  
  c_ci$ci <- paste0("(",c_ci$cil,"-",c_ci$ciu,")")
  
  c_ci <- c_ci %>% 
    select(risk_score_name,timepoint,subgroup_index,est,ci) %>%
    left_join(subgroup_matrix,by="subgroup_index")
   
}








########################################
#Createsubgroups

subgroup_function<-function(data,
                            subset_g,
                            subset2_g,
                            indexyear_g,
                            indexage_g,
                            sex_g,
                            charlson_cat_g,
                            charlson_g,
                            icd_n_chapters_g,
                            icd_chapters_g,
                            any_par_g){
  
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
  
  if(!is.null(icd_chapters_g)){
    stopifnot(length(icd_chapters_g)==1)
    stopifnot(icd_chapters_g %in% colnames(data))
    data<- data[data[[icd_chapters_g]],]
  
  }
  
  if(!is.null(any_par_g)){
    data<- data %>% filter(in_par %in% any_par_g)
  }
  
  return(data)
}


comorb_colnames<-list(c("LopNr","date","CCIw"),
                      c("LopNr","date","dci"),
                      c("LopNr","date","MDCI"),
                      c("LopNr","date","MDCI"))

risk_score_name<-c("CCI10","DCI","MDCI","MDCI_sv")

save_parts<-TRUE




















### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# All 1  - focus on calendar time 
sex_groups<-list("man",
                 "woman",
                 c("man","woman"))

age_groups<-list(c(18,120))

subgroup<-list(c("Development","Validation"))

years_groups <- c(list(c(2006:2022)),
                  as.list(c(2003:2022)))


subgroups<-create_subgroups(sex_groups=sex_groups,
                            age_groups=age_groups,
                            years=years_groups,
                            subgroup=subgroup,
                            subgroup2=NULL,
                            charlson=NULL)
length(subgroups$subgroups)
dim(subgroups$subgroups_matrix)


cindices <- compute_cindex(ncores=24,
                           cluster_indices=cluster_indices,
                           log_path=log_path,
                           data_path=data_path,
                           subgroups=subgroups$subgroups,
                           subgroup_function=subgroup_function,
                           at_times=c(1),
                           save_parts=save_parts,
                           save_path=save_path,
                           save_name="All1",
                           risk_score_name=risk_score_name)
gc()
save(cindices,file=paste0(save_path,"cindices_all1.Rdata"))
subgroups_matrix <- subgroups$subgroups_matrix
save(subgroups_matrix,file=paste0(save_path,"subgroup_matrix_all1.Rdata"))
table(cindices$point_ests$N_na)
table(cindices$point_ests_aggr$N_na)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 











### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# All 5  - focus on calendar time 
sex_groups<-list("man",
                 "woman",
                 c("man","woman"))

age_groups<-list(c(18,120))

subgroup<-list(c("Development","Validation"))

years_groups <- c(list(c(2006:2022)),
                  as.list(c(2003:2022)))


subgroups<-create_subgroups(sex_groups=sex_groups,
                            age_groups=age_groups,
                            years=years_groups,
                            subgroup=subgroup,
                            subgroup2=NULL,
                            charlson=NULL)
length(subgroups$subgroups)
dim(subgroups$subgroups_matrix)


cindices <- compute_cindex(ncores=32,
                           cluster_indices=cluster_indices,
                           log_path=log_path,
                           data_path=data_path,
                           subgroups=subgroups$subgroups,
                           subgroup_function=subgroup_function,
                           at_times=c(5),
                           save_parts=save_parts,
                           save_path=save_path,
                           save_name="All5",
                           risk_score_name=risk_score_name)
gc()
save(cindices,file=paste0(save_path,"cindices_all5.Rdata"))
subgroups_matrix <- subgroups$subgroups_matrix
save(subgroups_matrix,file=paste0(save_path,"subgroup_matrix_all5.Rdata"))
table(cindices$point_ests$N_na)
table(cindices$point_ests_aggr$N_na)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 














### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# All_ages 1 - focus on age and comorbidity 
sex_groups<-list("man",
                 "woman",
                 c("man","woman"))

age_groups<- as.list(c(30:99))
age_groups <- append(list(c(18,29)),age_groups)
age_groups <- append(age_groups,list(c(100,120)))

subgroup<-list(c("Development","Validation"))

years_groups <- list(c(2006:2022))

subgroups<-create_subgroups(sex_groups=sex_groups,
                            age_groups=age_groups,
                            years=years_groups,
                            subgroup=subgroup,
                            subgroup2=NULL)
length(subgroups$subgroups)
dim(subgroups$subgroups_matrix)


cindices <- compute_cindex(ncores=32,
                           cluster_indices=cluster_indices,
                           log_path=log_path,
                           data_path=data_path,
                           subgroups=subgroups$subgroups,
                           subgroup_function=subgroup_function,
                           at_times=c(1,5),
                           save_parts=save_parts,
                           save_path=save_path,
                           save_name="All_ages",
                           risk_score_name=risk_score_name)
gc()
save(cindices,file=paste0(save_path,"cindices_all_ages1.Rdata"))
subgroups_matrix <- subgroups$subgroups_matrix
save(subgroups_matrix,file=paste0(save_path,"subgroup_matrix_all_ages1.Rdata"))
table(cindices$point_ests$N_na)
table(cindices$point_ests_aggr$N_na)





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Charson
sex_groups<-list("man",
                 "woman",
                 c("man","woman"))

age_groups<-list(c(18,120),
                 c(18,49),
                 c(50,69),
                 c(70,120))

subgroup<-list(c("Development","Validation"))

years_groups <- c(list(c(2006:2022)))

charlson <- list(0,1,2,3)
charlson[[5]]<-4:20

subgroups<-create_subgroups(sex_groups=sex_groups,
                            age_groups=age_groups,
                            years=years_groups,
                            subgroup=subgroup,
                            subgroup2=NULL,
                            charlson=charlson,
                            icd_n_chapters=NULL)
length(subgroups$subgroups)
dim(subgroups$subgroups_matrix)


cindices <- compute_cindex(ncores=32,
                           cluster_indices=cluster_indices,
                           log_path=log_path,
                           data_path=data_path,
                           subgroups=subgroups$subgroups,
                           subgroup_function=subgroup_function,
                           at_times=c(1,5),
                           save_parts=save_parts,
                           save_path=save_path,
                           save_name="CCI",
                           risk_score_name=risk_score_name)
gc()
save(cindices,file=paste0(save_path,"cindices_CCI.Rdata"))
subgroups_matrix <- subgroups$subgroups_matrix
save(subgroups_matrix,file=paste0(save_path,"subgroup_matrix_CCI.Rdata"))
table(cindices$point_ests$N_na)
table(cindices$point_ests_aggr$N_na)























### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# ICD chapters 1
sex_groups<-list("man",
                 "woman",
                 c("man","woman"))

age_groups<-list(c(18,120),
                 c(18,49),
                 c(50,69),
                 c(70,120))

subgroup<-list(c("Development","Validation"))

years_groups <- c(list(c(2006:2022)))

icd_n_chapters <- list(0,1,2,3,4,5,c(6:20))


subgroups<-create_subgroups(sex_groups=sex_groups,
                            age_groups=age_groups,
                            years=years_groups,
                            subgroup=subgroup,
                            subgroup2=NULL,
                            charlson=NULL,
                            icd_n_chapters=icd_n_chapters)
length(subgroups$subgroups)
dim(subgroups$subgroups_matrix)


cindices <- compute_cindex(ncores=24,
                           cluster_indices=cluster_indices,
                           log_path=log_path,
                           data_path=data_path,
                           subgroups=subgroups$subgroups,
                           subgroup_function=subgroup_function,
                           at_times=c(1),
                           save_parts=save_parts,
                           save_path=save_path,
                           save_name="ICD_chapters1",
                           risk_score_name=risk_score_name)
gc()
save(cindices,file=paste0(save_path,"cindices_ICD_chapters1.Rdata"))
subgroups_matrix <- subgroups$subgroups_matrix
save(subgroups_matrix,file=paste0(save_path,"subgroup_matrix_ICD_chapters1.Rdata"))
table(cindices$point_ests$N_na)
table(cindices$point_ests_aggr$N_na)



cindices <- compute_cindex(ncores=24,
                           cluster_indices=cluster_indices,
                           log_path=log_path,
                           data_path=data_path,
                           subgroups=subgroups$subgroups,
                           subgroup_function=subgroup_function,
                           at_times=c(5),
                           save_parts=save_parts,
                           save_path=save_path,
                           save_name="ICD_chapters5",
                           risk_score_name=risk_score_name)
gc()
save(cindices,file=paste0(save_path,"cindices_ICD_chapters5.Rdata"))
subgroups_matrix <- subgroups$subgroups_matrix
save(subgroups_matrix,file=paste0(save_path,"subgroup_matrix_ICD_chapters5.Rdata"))
table(cindices$point_ests$N_na)
table(cindices$point_ests_aggr$N_na)



###
# Not in par
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
                            any_par=FALSE)
length(subgroups$subgroups)
dim(subgroups$subgroups_matrix)


cindices <- compute_cindex(ncores=24,
                           cluster_indices=cluster_indices,
                           log_path=log_path,
                           data_path=data_path,
                           subgroups=subgroups$subgroups,
                           subgroup_function=subgroup_function,
                           at_times=c(1),
                           save_parts=save_parts,
                           save_path=save_path,
                           save_name="NoPAR",
                           risk_score_name=risk_score_name[2])
gc()
save(cindices,file=paste0(save_path,"cindices_NoPAR.Rdata"))
subgroups_matrix <- subgroups$subgroups_matrix
save(subgroups_matrix,file=paste0(save_path,"subgroup_matrix_NoPAR.Rdata"))
table(cindices$point_ests$N_na)
table(cindices$point_ests_aggr$N_na)


cindices <- compute_cindex(ncores=24,
                           cluster_indices=cluster_indices,
                           log_path=log_path,
                           data_path=data_path,
                           subgroups=subgroups$subgroups,
                           subgroup_function=subgroup_function,
                           at_times=c(5),
                           save_parts=save_parts,
                           save_path=save_path,
                           save_name="NoPAR5",
                           risk_score_name=risk_score_name[2])
gc()
save(cindices,file=paste0(save_path,"cindices_NoPAR5.Rdata"))
subgroups_matrix <- subgroups$subgroups_matrix
save(subgroups_matrix,file=paste0(save_path,"subgroup_matrix_NoPAR5.Rdata"))
table(cindices$point_ests$N_na)
table(cindices$point_ests_aggr$N_na)










### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Save as excel file

###
# All
load(file=paste0(save_path,"cindices_all1.Rdata")) # cindices
loads <- load(file=paste0(save_path,"subgroup_matrix_all1.Rdata")) # subgroup_matrix,
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)

cindices_table <- fix_table(cindices,subgroups_matrix,dig=dig)

temp <- split(cindices_table,f=cindices_table$risk_score_name)
openxlsx::write.xlsx(temp,file=paste0(results_path,"cindex_all1.xlsx"))

load(file=paste0(save_path,"cindices_all5.Rdata")) # cindices
loads <- load(file=paste0(save_path,"subgroup_matrix_all5.Rdata")) # subgroup_matrix,
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)

cindices_table <- fix_table(cindices,subgroups_matrix,dig=dig)

temp <- split(cindices_table,f=cindices_table$risk_score_name)
openxlsx::write.xlsx(temp,file=paste0(results_path,"cindex_all5.xlsx"))

###
# All ages
load(file=paste0(save_path,"cindices_all_ages1.Rdata")) # cindices
loads <- load(file=paste0(save_path,"subgroup_matrix_all_ages1.Rdata")) # subgroup_matrix,
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)

cindices_table <- fix_table(cindices,subgroups_matrix,dig=dig)

temp <- split(cindices_table,f=cindices_table$risk_score_name)
openxlsx::write.xlsx(temp,file=paste0(results_path,"cindex_AllAges.xlsx"))




### 
# CCI  cindices_CCI.Rdata
load(file=paste0(save_path,"cindices_CCI.Rdata")) # cindices
load(file=paste0(save_path,"subgroup_matrix_CCI.Rdata")) # subgroup_matrix,
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)

cindices_table <- fix_table(cindices,subgroups_matrix,dig=dig)

temp <- split(cindices_table,f=cindices_table$risk_score_name)
openxlsx::write.xlsx(temp,file=paste0(results_path,"cindex_CCI.xlsx"))













###
# ICD chapters
load(file=paste0(save_path,"cindices_ICD_chapters1.Rdata")) # cindices
load(file=paste0(save_path,"subgroup_matrix_ICD_chapters1.Rdata")) # subgroup_matrix,
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)

cindices_table <- fix_table(cindices,subgroups_matrix,dig=dig)

temp <- split(cindices_table,f=cindices_table$risk_score_name)
openxlsx::write.xlsx(temp,file=paste0(results_path,"cindex_ICD_chapters1.xlsx"))


load(file=paste0(save_path,"cindices_ICD_chapters5.Rdata")) # cindices
load(file=paste0(save_path,"subgroup_matrix_ICD_chapters5.Rdata")) # subgroup_matrix,
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)

cindices_table <- fix_table(cindices,subgroups_matrix,dig=dig)

temp <- split(cindices_table,f=cindices_table$risk_score_name)
openxlsx::write.xlsx(temp,file=paste0(results_path,"cindex_ICD_chapters5.xlsx"))



###
# Not in PAR
load(file=paste0(save_path,"cindices_NoPAR.Rdata")) # cindices
loads <- load(file=paste0(save_path,"subgroup_matrix_NoPAR.Rdata")) # subgroup_matrix,
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)

cindices_table <- fix_table(cindices,subgroups_matrix,dig=dig)

temp <- split(cindices_table,f=cindices_table$risk_score_name)
openxlsx::write.xlsx(temp,file=paste0(results_path,"cindex_NoPAR.xlsx"))

# Not in PAR
load(file=paste0(save_path,"cindices_NoPAR5.Rdata")) # cindices
loads <- load(file=paste0(save_path,"subgroup_matrix_NoPAR5.Rdata")) # subgroup_matrix,
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)

cindices_table <- fix_table(cindices,subgroups_matrix,dig=dig)

temp <- split(cindices_table,f=cindices_table$risk_score_name)
openxlsx::write.xlsx(temp,file=paste0(results_path,"cindex_NoPAR5.xlsx"))





