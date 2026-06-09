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
ages <- 18:100 # 18:100

ncores <- 43

maxfu <- 10

print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("X","data_path","save_path","cluster_ids","maxfu"),envir = environment())

crude <- parLapplyLB(cl = cl,
                     ages,
                     fun = function(a){
                     require(dplyr)
                     setwd(paste0(X,"ComorbidityBase\\Works\\2024\\2024 Cause of death\\"))
                     
                     id <- 1
                     sink(paste0("C:\\Marcus\\Misc\\","LE_est_agemodel",a,".txt"))
                     print(a)
                     print(Sys.time())
                     loads <- load(file=paste0("C:\\Marcus\\ComorbidityBase_cache\\cox_life_expectancy_",a,"_",id,".Rdata"))
                     
                     data_list_men <- data_temp %>% 
                       select(LopNr,sex,age) %>%
                       filter(sex %in% "man") %>%
                       slice(1)
                     
                     modelname <- paste0(a,"_","CRUDE",
                                         "",
                                         10,
                                         "",
                                         "man")
                     
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                     
                     data_list_men$cumhaz <- predict(model,
                                                     newdata=data_list_men,
                                                     type="cumhaz",times=1)$.pred_cumhaz
                     
                     
                     data_list_women <- data_temp %>% 
                       select(LopNr,sex,age) %>%
                       filter(sex %in% "woman") %>%
                       slice(1)
                     
                     modelname <- paste0(a,"_","CRUDE",
                                         "",
                                         10,
                                         "",
                                         "woman")
                     
                     load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))
                     
                     data_list_women$cumhaz <- predict(model,
                                                       newdata=data_list_women,
                                                       type="cumhaz",times=1)$.pred_cumhaz
                     
                     return(list(men=data_list_men,
                                 women=data_list_women))
                   })

print( Sys.time() - start_time)
stopCluster(cl)
gc(verbose = FALSE)


crude_m <- do.call(rbind,lapply(FUN=function(x) x$men,crude)) %>%
  select(age,cumhaz)

crude_w <- do.call(rbind,lapply(FUN=function(x) x$women,crude))%>%
  select(age,cumhaz)

# Obtain 10-year estimate in 100-year olds

a <- 100
id <- 1
loads <- load(file=paste0("C:\\Marcus\\ComorbidityBase_cache\\cox_life_expectancy_",a,"_",id,".Rdata"))

data_list_men <- data_temp %>% 
  select(LopNr,sex,age) %>%
  filter(sex %in% "man") %>%
  slice(1)

modelname <- paste0(a,"_","CRUDE",
                    "",
                    10,
                    "",
                    "man")

loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))


p100_men <- predict(model,
             newdata=data_list_men,
             type="cumhaz",times=1:10)[[1]][[1]]
p100_men <- diff(c(0,p100_men$.pred_cumhaz))

data_list_women <- data_temp %>% 
  select(LopNr,sex,age) %>%
  filter(sex %in% "woman") %>%
  slice(1)

modelname <- paste0(a,"_","CRUDE",
                    "",
                    10,
                    "",
                    "woman")

load(file=paste0("C:\\Marcus\\LEcox\\LE_model_",modelname,".Rdata"))

p100_women<- predict(model,
                     newdata=data_list_women,
                     type="cumhaz",times=1:10)[[1]][[1]]
p100_women <- diff(c(0,p100_women$.pred_cumhaz))

load(file=paste0(save_path,"scb_hazard_women.Rdata"))
load(file=paste0(save_path,"scb_hazard_men.Rdata"))

plot(crude_m$cumhaz,
     scb_hazard_men$y[scb_hazard_men$x %in% 18:100])
lines(c(0,1),c(0,1))

summary(crude_m$cumhaz-scb_hazard_men$y[scb_hazard_men$x %in% 18:100])

plot(crude_w$cumhaz,
     scb_hazard_women$y[scb_hazard_women$x %in% 18:100])
lines(c(0,1),c(0,1))
summary(crude_w$cumhaz-scb_hazard_women$y[scb_hazard_women$x %in% 18:100])

###
# Extrapolate 
crude_m <- list(x=crude_m$age,
                y=crude_m$cumhaz)

crude_w <- list(x=crude_w$age,
                y=crude_w$cumhaz)


plot(crude_m$x,
     crude_m$y )
plot(crude_w$x,
     crude_w$y)

times_for_slope <- 0:1

###
# extrapolate after 101
p100_men
scb_hazard_men$y[scb_hazard_men$x %in% 101:110]

p100_women
scb_hazard_women$y[scb_hazard_women$x %in% 101:110]

crude_m$x <- c(crude_m$x,101:130)
crude_w$x <- c(crude_w$x,101:130)

crude_m$y <- c(crude_m$y,
               scb_hazard_men$y[scb_hazard_men$x %in% 101:130] )
crude_w$y <- c(crude_w$y,
               scb_hazard_women$y[scb_hazard_women$x %in% 101:130] )


plot(crude_m$x,
     crude_m$y,
     ylim=c(0,2) )
lines( scb_hazard_men$x, scb_hazard_men$y,col="red")

plot(crude_w$x,
     crude_w$y,
     ylim=c(0,2))
lines( scb_hazard_women$x, scb_hazard_women$y,col="red")

plot(crude_m$y,
     scb_hazard_men$y)
lines(c(0,10),c(0,10))
plot(crude_w$y,
     scb_hazard_women$y)
lines(c(0,10),c(0,10))

save(crude_w,file=paste0(save_path,"crude_hazard_women.Rdata"))
save(crude_m,file=paste0(save_path,"crude_hazard_men.Rdata"))

summary(diff(crude_m$y))
summary(diff(crude_w$y))
plot(diff(crude_m$y))
plot(diff(crude_w$y))
