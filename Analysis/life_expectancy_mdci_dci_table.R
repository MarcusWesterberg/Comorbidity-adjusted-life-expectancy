

###
# Compress life expectancy estimates using a GAM

require(parallel)

start_time <- Sys.time()
cluster_indices <- 18:100
cluster_ids <- 1:64

print(cluster_ids)
ncores <- 32

print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("data_path","save_path","cluster_ids"),envir = environment())

ret <- parLapplyLB(cl = cl,
                   cluster_indices,
                   fun = function(a){
                     require(dplyr)
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_age",a,".Rdata"))
              
                     LE$sex[LE$sex %in% "men"]<-"man"
                     LE$sex[LE$sex %in% "women"]<-"woman"
                     LE$age <- a
                  
                     LE <- LE %>%
                       select(sex,age,MDCI,dci,le)  %>% 
                       mutate(MDCI=round(MDCI,dig=1),
                              dci=round(dci,dig=1)) %>%
                       group_by(sex,age,MDCI,dci) %>%
                       mutate(n=n(),
                              le=mean(le)) %>%
                       slice(1) %>%
                       ungroup()
                     
                     return(LE)
                     
                   })

stopCluster(cl)
gc(verbose = FALSE)
require(tidyverse)

le.table <- do.call(rbind,ret)

dim(le.table)

le.table <- le.table %>%
  arrange(sex,age,MDCI,dci)

###
# Save  as table
save(le.table,file="X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_mdci_dci_table.Rdata")

require(mgcv)
ERL.men <- le.table %>% filter(sex == "man") %>% mutate(le=log(le))
m.men <- gam(le ~ te(age,MDCI,dci),
             data=ERL.men
             ,weights=ERL.men$n,
             family="gaussian")
summary(m.men)

ERL.women <- le.table %>% filter(sex == "woman") %>% mutate(le=log(le))
m.women <- gam(le ~ te(age,MDCI,dci),
               data=ERL.women,
               weights=ERL.women$n,
               family="gaussian")
summary(m.women)

m.men$residuals <- c()
m.men$fitted.values <- c()
m.men$linear.predictors <- c()
m.men$model <- c()
m.men$y <- c()


m.women$residuals <- c()
m.women$fitted.values <- c()
m.women$linear.predictors <- c()
m.women$model <- c()
m.women$y <- c()

le.models <- list("men"=m.men,
                  "women"=m.women)
save(le.models,file="X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_mdci_dci_models.Rdata")





###
# End
###