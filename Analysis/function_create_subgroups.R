create_subgroups <- function(sex_groups=list(c("man","woman")),
                             age_groups= list(c(60,69)),
                             years=list(c(2003:2012)),
                             subgroup=list(c("Development")),
                             subgroup2=NULL,
                             charlson=NULL,
                             charlson_cat=NULL,
                             icd_n_chapters=NULL,
                             icd_chapters=NULL,
                             any_par=NULL){
  
  
  subgroups <- list()
  
  for(sx in 1:length(sex_groups)){
    for(a in 1:length(age_groups) ){
      for(y in 1:length(years)){
        for(sg in 1:length(subgroup)){
          for(sg2 in 1:length(subgroup2)){
            for(cci in 1:length(charlson)){
              for(cc in 1:length(charlson_cat)){
                for(icd_n_c in 1:length(icd_n_chapters)){
                  for(icd_c in 1:length(icd_chapters)){
                    for(ap_c in 1:length(any_par)){
                      if(sx>0 & a>0 & y>0 & sg>0 & sg2>0 & cci>0 & cc>0 & icd_n_c>0 & icd_c>0 & ap_c>0){
                        tempname <- paste0("AG",a,"Y",y,"Sex",sx,"Sg",sg,"SgT",sg2,"CCI",cci,"CC",cc,"ICD",icd_n_c,"ICDc",icd_c,"AP"=ap_c)
                        subgroups[[tempname]] <- list("indexyear"=years[[y]],
                                                      "indexage"=age_groups[[a]],
                                                      "sex"=sex_groups[[sx]],
                                                      "subgroup"=subgroup[[sg]],
                                                      "subgroup2"=subgroup2[[sg2]],
                                                      "charlson"=charlson[[cci]],
                                                      "charlson_cat"=charlson_cat[[cc]],
                                                      "icd_n_chapters"=icd_n_chapters[[icd_n_c]],
                                                      "icd_chapters"=icd_chapters[[icd_c]],
                                                      "any_par"=any_par[[ap_c]])
                      }
                    }
               
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  subgroups_matrix <- subgroups
  subgroups_matrix <- lapply(FUN=function(x){
    x$indexyear <- paste0(x$indexyear,collapse=", ")
    x$indexage <- paste0(x$indexage,collapse="-")
    x$sex <- paste0(x$sex,collapse=", ")
    x$subgroup <- paste0(x$subgroup,collapse=", ")
    x$subgroup2 <- paste0(x$subgroup2,collapse=", ")
    x$charlson <- paste0(x$charlson,collapse=", ")
    x$charlson_cat <- paste0(x$charlson_cat,collapse=", ")
    x$icd_n_chapters <- paste0(x$icd_n_chapters,collapse=", ")
    x$icd_chapters <- paste0(x$icd_chapters,collapse=", ")
    x$any_par <- paste0(x$any_par,collapse=", ")
    return( do.call(cbind,x))
  }  , subgroups_matrix)
  subgroups_matrix <- do.call(rbind,subgroups_matrix)
  subgroups_matrix <- as.data.frame(subgroups_matrix)
  subgroups_matrix$subgroup_name <- names(subgroups)
  
  
  return(list("subgroups"=subgroups,
              "subgroups_matrix"=subgroups_matrix))
}
