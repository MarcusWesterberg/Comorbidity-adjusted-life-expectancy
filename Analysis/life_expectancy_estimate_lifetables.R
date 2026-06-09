###
#
X <- "\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"
require(tidyverse)

log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"
background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")
figs_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")


# LE from SCB - life tables
load("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\LE.Rdata")
load(file="\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\Popsize.Rdata")

D_LE <- D_LE %>%
  filter(year %in% 2006:2022) %>%
  left_join(P,by=c("age","sex","year"))
rm(P)

D_LE <- D_LE %>%
  group_by(age,sex) %>%
  summarise(LE=sum(N*LE)/sum(N),.groups="keep") %>%
  ungroup() %>%
  rename(le_scb=LE) 

D_LE$sex[D_LE$sex %in% "men"]<-"man"
D_LE$sex[D_LE$sex %in% "women"]<-"woman"

D_LE %>% filter(age==100)
###
# Extract hazard functions from scb data
loads <- load(file="\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\Popsize.Rdata")
loads <- load(file="\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\risk.Rdata")

D <- D %>% 
  filter(year %in% 2006:2022) %>%
  left_join(P,by=c("age","sex","year"))

###
# Extract hazard functions
scb_hazard <- D %>%
  group_by(sex,age) %>%
  summarise(risk=sum(risk*N)/sum(N),.groups="keep") %>%
  ungroup()

rm(D)

scb_hazard_men <- scb_hazard %>%
  filter(sex %in% "men")

scb_hazard_women <- scb_hazard %>%
  filter(sex %in% "women")

require(mgcv)


### Model only used for extrapolation
#men
scbm <- gam(log(risk) ~ s(age),
            data=scb_hazard_men)
scb_hazard_men$risk <- exp(scbm$fitted.values)

# women
scbw <- gam(log(risk) ~ s(age),
            data=scb_hazard_women)
scb_hazard_women$risk <- exp(scbw$fitted.values)


### 
# predictions
xout <- seq(18,110,by=1)

scb_hazard_men <- list("x"=xout,
                       "y"=exp(predict(scbm,newdata = data.frame(age=xout))))
scb_hazard_women <- list("x"=xout,
                         "y"=exp(predict(scbw,newdata = data.frame(age=xout))))

scb_hazard_men$y <- - log(1-scb_hazard_men$y) 
scb_hazard_women$y <- - log(1-scb_hazard_women$y) 

# slope before interpolation 110-130 years
times_for_slope <- 0:1

delta_men <- mean(diff(scb_hazard_men$y[ rev(length(scb_hazard_men$y)-times_for_slope) ]) )
delta_women <- mean(diff(scb_hazard_women$y[ rev(length(scb_hazard_women$y)-times_for_slope) ]) )

###
# extrapolate after 110
xout <- seq(110,130,by=1)
xout <- xout[-1]

delta_men <- delta_men*(1:length(xout))
delta_women <- delta_women*(1:length(xout))

scb_hazard_men$x <- c(scb_hazard_men$x,xout)
scb_hazard_women$x <- c(scb_hazard_women$x,xout)

scb_hazard_men$y <- c(scb_hazard_men$y,
                      scb_hazard_men$y[length( scb_hazard_men$y )] + delta_men )
scb_hazard_women$y <- c(scb_hazard_women$y,
                        scb_hazard_women$y[length( scb_hazard_women$y )] + delta_women)

save(scb_hazard_men,file=paste0(save_path,"scb_hazard_men.Rdata"))
save(scb_hazard_women,file=paste0(save_path,"scb_hazard_women.Rdata"))



###
rm(scbm)
rm(scbw)
scb_hazard <- scb_hazard_men
extrapolate_cumhaz <- function(age,
                               scb_hazard,
                               plotit=FALSE){
  
  ###
  # Find average hazard 
  scb_hazard_avr <- scb_hazard
  inds <- scb_hazard_avr$x>= age
  scb_hazard_avr$x <- scb_hazard_avr$x[inds]
  scb_hazard_avr$y <- scb_hazard_avr$y[inds]
  

  
  survcurv <- c(1,exp(-cumsum(scb_hazard_avr$y))) # add time zero (surv=1)
  
  le <- (1 / 2) * (2*sum(survcurv) - survcurv[1] - survcurv[length(survcurv)]) 
  
  if(plotit){
    
    xx <- seq(0,length(survcurv)-1,by=1) + a
    plot(xx,survcurv,col="darkgrey",type="l",lwd=2,xlab="Time",ylab="Survival")
    
    mtext(side=3,text=paste0("LE: ",round(le,1) ),line=1)
  }
  
  return(le)
}

le_men <- sapply(FUN=function(x) extrapolate_cumhaz(age=x,scb_hazard=scb_hazard_men), X=18:100)
le_women <- sapply(FUN=function(x) extrapolate_cumhaz(age=x,scb_hazard=scb_hazard_women), X=18:100)

le_men <- data.frame(age=18:100,
                     "est"=le_men,
                     "scb"=D_LE$le_scb[D_LE$sex == "man" & D_LE$age %in% 18:100])

le_women <- data.frame(age=18:100,
                     "est"=le_women,
                     "scb"=D_LE$le_scb[D_LE$sex == "woman" & D_LE$age %in% 18:100])

plot(le_men$est,type="l")
lines(le_men$scb,type="l",col="red")

plot(le_men$est-le_men$scb,type="l")


plot(le_women$est,type="l")
lines(le_women$scb,type="l",col="red")

plot(le_women$est-le_women$scb,type="l")
summary(le_women$est-le_women$scb)

