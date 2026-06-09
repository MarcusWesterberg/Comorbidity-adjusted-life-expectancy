###
# Main
X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

ncores<-64
log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"
background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

comorb_data_path <- "C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")

source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\function_cindex.R"))
source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\function_prepare_data_cindex.R"))
source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\function_prepare_data_charlson_cat.R"))
source(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Scripts\\function_create_subgroups.R"))
cluster_indices <- 1:64

data_name<-c("Development_final","Validation_final")

comorb_data_name<-list(c("Development_CCI_lookback_10",
                         "Development_DCI",
                         "Development_MDCI",
                         "Development_sv_MDCI",
                         "Development_elixhauser"),
                       
                       c("Validation_CCI_lookback_10",
                         "Validation_DCI",
                         "Validation_MDCI",
                         "Validation_sv_MDCI",
                         "Validation_elixhauser"))

comorb_colnames<-list(c("LopNr","date","CCIw"),
                      c("LopNr","date","dci"),
                      c("LopNr","date","MDCI"),
                      c("LopNr","date","MDCI"),
                      c("LopNr","date","score"))

risk_score_name<-c("CCI10","DCI","MDCI","MDCI_sv","Elixhauser")


########################################
#Preapredata

# only run this ify ou want to prepared ata again...
notrun<-function(){
  comorb_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"
  prepare_data_cindex(ncores=16,
                      cluster_indices=cluster_indices,
                      log_path=log_path,
                      comorb_data_path=comorb_data_path,
                      comorb_data_name=comorb_data_name,
                      comorb_colnames=comorb_colnames,
                      comorb_newcolnames=c("LopNr","indexdate","comorbidity_index"),
                      data_path=data_path,
                      data_name=data_name,
                      background_data_path=background_data_path,
                      risk_score_name=risk_score_name,
                      save_name="SaveParts")#savenameonlyforlog
  
  prepare_data_charlson_cat(ncores=16,
                            cluster_indices=cluster_indices,
                            log_path=log_path,
                            comorb_data_path=comorb_data_path,
                            comorb_data_name=list(comorb_data_name[[1]][1]),
                            data_path=data_path,
                            #data_name=data_name,
                            save_name="charlson")#savenameonlyforlog
  
  prepare_data_charlson_cat(ncores=16,
                            cluster_indices=cluster_indices,
                            log_path=log_path,
                            comorb_data_path=comorb_data_path,
                            comorb_data_name=list("Development_sv_CCI_lookback_10",
                                                  "Validation_sv_CCI_lookback_10"),
                            data_path=data_path,
                            #data_name=data_name,
                            save_name="charlson_sv")#savenameonlyforlog
  
  
  gc()
}














###
#End
###