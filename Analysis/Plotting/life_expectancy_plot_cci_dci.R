

###
# Illustrate life expectancy
require(parallel)

start_time <- Sys.time()
cluster_indices <- 18:100
cluster_ids <- 1:64

print(cluster_ids)
ncores <- 32

print(paste0("Ncores=",ncores))

cl <- makeCluster(ncores)

clusterExport(cl, c("X","cluster_ids"),envir = environment())

ret <- parLapplyLB(cl = cl,
                   cluster_indices,
                   fun = function(a){
                     require(dplyr)
                     
                     # for crude 
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_age",a,".Rdata"))
                     LE_crude <- LE %>%
                       select(LopNr,rowid,le.crude)
                     
                     # for ccidci
                     loads <- load(file=paste0("C:\\Marcus\\LEcox\\LE_age_ccidci",a,".Rdata"))
                     
                     LE <- LE %>%
                       left_join(LE_crude,by=c("LopNr","rowid"))
                     gc()
                     
                     load("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\LE.Rdata")
                     load(file="\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\Popsize.Rdata")
                     
                     D_LE <- D_LE %>%
                       filter(year %in% 2006:2022) %>%
                       left_join(P,by=c("age","sex","year"))
                     
                     D_LE <- D_LE %>%
                       filter(age %in% a) %>%
                       group_by(sex) %>%
                       summarise(LE=sum(N*LE)/sum(N),.groups="keep") %>%
                       ungroup() %>%
                       rename(le_scb=LE) 
                     
                     D_LE$sex[D_LE$sex %in% "men"]<-"man"
                     D_LE$sex[D_LE$sex %in% "women"]<-"woman"
                     LE$age <- a
                     LE <- LE %>%
                       left_join(D_LE,by=c("sex"))
                     
                     ret_overall <- LE %>%
                       rename(le=le.ccidci) %>%
                       
                       summarise(n=n(),
                                 le_min=min(le),
                                 le_05=quantile(le,probs=0.05),
                                 le_25=quantile(le,probs=0.25),
                                 le_50=quantile(le,probs=0.5),
                                 le_75=quantile(le,probs=0.75),
                                 le_95=quantile(le,probs=0.95),
                                 le_max=max(le),
                                 le_mean=mean(le),
                                 
                                 le_scb_mean=mean(le_scb),
                                 le_scb_median=median(le_scb)) 
                     ret_overall$age <- a
                     
                     ret_all <- LE %>%
                       rename(le=le.ccidci) %>%
                       group_by(sex) %>%
                       summarise(n=n(),
                                 le_min=min(le),
                                 le_05=quantile(le,probs=0.05),
                                 le_25=quantile(le,probs=0.25),
                                 le_50=quantile(le,probs=0.5),
                                 le_75=quantile(le,probs=0.75),
                                 le_95=quantile(le,probs=0.95),
                                 le_max=max(le),
                                 le_mean=mean(le),
                              
                                 le_scb_mean=mean(le_scb),
                                 le_scb_median=median(le_scb),
                                 .groups="keep") %>%
                       ungroup()
                     ret_all$age <- a
                     
                     ret_cci <- LE %>%
                       rename(le=le.ccidci) %>%
                       group_by(sex,CCI) %>%
                       summarise(n=n(),
                                 le_min=min(le),
                                 le_05=quantile(le,probs=0.05),
                                 le_25=quantile(le,probs=0.25),
                                 le_50=quantile(le,probs=0.5),
                                 le_75=quantile(le,probs=0.75),
                                 le_95=quantile(le,probs=0.95),
                                 le_max=max(le),
                                 le_mean=mean(le),
                               
                                 le_scb_mean=mean(le_scb),
                                 le_scb_median=median(le_scb),
                                 .groups="keep") %>%
                       ungroup()
                     ret_cci$age <- a
                     
                     
                     ret_overall_d <- LE %>%
                       rename(le=le.ccidci) %>%
                       mutate(le=le-le_scb) %>%
                       summarise(n=n(),
                                 le_min=min(le),
                                 le_05=quantile(le,probs=0.05),
                                 le_25=quantile(le,probs=0.25),
                                 le_50=quantile(le,probs=0.5),
                                 le_75=quantile(le,probs=0.75),
                                 le_95=quantile(le,probs=0.95),
                                 le_max=max(le),
                                 le_mean=mean(le),
                                 
                                 le_scb_mean=mean(le_scb),
                                 le_scb_median=median(le_scb)) 
                     ret_overall_d$age <- a
                     
                     ret_all_d <- LE %>%
                       rename(le=le.ccidci) %>%
                       group_by(sex) %>%
                       
                       mutate(le=le-le_scb) %>%
                       
                       summarise(n=n(),
                                 le_min=min(le),
                                 le_05=quantile(le,probs=0.05),
                                 le_25=quantile(le,probs=0.25),
                                 le_50=quantile(le,probs=0.5),
                                 le_75=quantile(le,probs=0.75),
                                 le_95=quantile(le,probs=0.95),
                                 le_max=max(le),
                                 le_mean=mean(le),
                            
                                 le_scb_mean=mean(le_scb),
                                 le_scb_median=median(le_scb),
                                 .groups="keep") %>%
                       ungroup()
                     ret_all_d$age <- a
                     
                     ret_cci_d <- LE %>%
                       rename(le=le.ccidci) %>%
                       mutate(le=le-le_scb) %>%
                       group_by(sex,CCI) %>%
                       summarise(n=n(),
                                 le_min=min(le),
                                 le_05=quantile(le,probs=0.05),
                                 le_25=quantile(le,probs=0.25),
                                 le_50=quantile(le,probs=0.5),
                                 le_75=quantile(le,probs=0.75),
                                 le_95=quantile(le,probs=0.95),
                                 le_max=max(le),
                                 le_mean=mean(le),
                            
                                 
                                 le_scb_mean=mean(le_scb),
                                 le_scb_median=median(le_scb),
                                 .groups="keep") %>%
                       ungroup()
                     ret_cci_d$age <- a
                     
                     ### VS Crude
                     ret_all_dd <- LE %>%
                       rename(le=le.ccidci) %>%
                       mutate(le=le-le.crude) %>%
                       group_by(sex) %>%
                       summarise(n=n(),
                                 le_min=min(le),
                                 le_05=quantile(le,probs=0.05),
                                 le_25=quantile(le,probs=0.25),
                                 le_50=quantile(le,probs=0.5),
                                 le_75=quantile(le,probs=0.75),
                                 le_95=quantile(le,probs=0.95),
                                 le_max=max(le),
                                 le_mean=mean(le),
                                 
                                 le_scb_mean=mean(le_scb),
                                 le_scb_median=median(le_scb),
                                 .groups="keep") %>%
                       ungroup()
                     ret_all_dd$age <- a
                     
                     ret_cci_dd <- LE %>%
                       rename(le=le.ccidci) %>%
                       mutate(le=le-le.crude) %>%
                       group_by(sex,CCI) %>%
                       summarise(n=n(),
                                 le_min=min(le),
                                 le_05=quantile(le,probs=0.05),
                                 le_25=quantile(le,probs=0.25),
                                 le_50=quantile(le,probs=0.5),
                                 le_75=quantile(le,probs=0.75),
                                 le_95=quantile(le,probs=0.95),
                                 le_max=max(le),
                                 le_mean=mean(le),
                                 
                                 le_scb_mean=mean(le_scb),
                                 le_scb_median=median(le_scb),
                                 .groups="keep") %>%
                       ungroup()
                     ret_cci_dd$age <- a
                     
                     ret_overall_dd <- LE %>%
                       rename(le=le.ccidci) %>%
                       mutate(le=le-le.crude) %>%
                       summarise(n=n(),
                                 le_min=min(le),
                                 le_05=quantile(le,probs=0.05),
                                 le_25=quantile(le,probs=0.25),
                                 le_50=quantile(le,probs=0.5),
                                 le_75=quantile(le,probs=0.75),
                                 le_95=quantile(le,probs=0.95),
                                 le_max=max(le),
                                 le_mean=mean(le),
                                 
                                 le_scb_mean=mean(le_scb),
                                 le_scb_median=median(le_scb)) 
                     ret_overall_dd$age <- a
                     
                     
                     ret <- list("overall"=ret_overall,
                                 "all"=ret_all,
                                 "cci"=ret_cci,
                                 "overall_diff"=ret_overall_d,
                                 "all_diff"=ret_all_d,
                                 "cci_diff"=ret_cci_d,
                                 "overall_diff2"=ret_overall_dd,
                                 "all_diff2"=ret_all_dd,
                                 "cci_diff2"=ret_cci_dd)
                     
                     return(ret)
                     
                   })

stopCluster(cl)
gc(verbose = FALSE)
save(ret,file="X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_dist_ccidci.Rdata")


load(file="X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_dist_ccidci.Rdata")


ret_overall <- do.call(rbind,lapply(FUN=function(x) x$overall,ret))
ret_all <- do.call(rbind,lapply(FUN=function(x) x$all,ret))
ret_cci <- do.call(rbind,lapply(FUN=function(x) x$cci,ret))

ret_overall_d <- do.call(rbind,lapply(FUN=function(x) x$overall_diff,ret))
ret_all_d <- do.call(rbind,lapply(FUN=function(x) x$all_diff,ret))
ret_cci_d <- do.call(rbind,lapply(FUN=function(x) x$cci_diff,ret))


ret_overall_d2 <- do.call(rbind,lapply(FUN=function(x) x$overall_diff2,ret))
ret_all_d2 <- do.call(rbind,lapply(FUN=function(x) x$all_diff2,ret))
ret_cci_d2 <- do.call(rbind,lapply(FUN=function(x) x$cci_diff2,ret))


###
# CCI and DCI
svg("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_ccidci.svg",width=10,height=4)
par(mfrow=c(1,2),mar=c(3,3,4,1))

temp <- ret_all_d2 %>%
  filter(sex  %in% "man" )


plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="Men",
     type="l",
     lwd=2,
     xlim=c(18,100),
     ylim=c(-8,8),
     axes=FALSE)
axis(side=1,at=c(18,seq(20,100,by=10)))
axis(side=2,at=seq(-8,8,by=1))

polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))

lines(temp$age,temp$le_05,lty=3,col="darkgrey",lwd=3)
lines(temp$age,temp$le_95,lty=3,col="darkgrey",lwd=3)


legend(x=20,
       y=8,
       legend=c("Q 5%","Mean","Q 95%"),
       lty=c(3,1,3),
       col=c("darkgrey","black","darkgrey"),
       lwd=3,
       bty="n")


mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Difference (years)",line=2)

temp <- ret_all_d2 %>%
  filter(sex  %in% "woman" )
plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="Women",
     type="l",
     lwd=2,
     xlim=c(18,100),
     ylim=c(-8,8),
     axes=FALSE)
axis(side=1,at=c(18,seq(20,100,by=10)))
axis(side=2,at=seq(-8,8,by=1))

polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))

lines(temp$age,temp$le_05,lty=3,col="darkgrey",lwd=3)
lines(temp$age,temp$le_95,lty=3,col="darkgrey",lwd=3)

mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Difference (years)",line=2)

dev.off()




###
# cci and DCI
svg("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_ccidci_overall.svg",width=6,height=5)
par(mfrow=c(1,2),mar=c(3,3,4,1))

temp <- ret_overall_d2 

plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="",
     type="l",
     lwd=2,
     xlim=c(40,100),
     ylim=c(-8,8),
     axes=FALSE,
     col="black")
axis(side=1,at=seq(40,100,by=10))
axis(side=2,at=c(-8:8))
polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))


lines(temp$age,temp$le_05,lty=3,col="grey",lwd=3)
lines(temp$age,temp$le_95,lty=3,col="grey",lwd=3)


legend(x=40,
       y=8,
       legend=c("Q 5%","Mean","Q 95%"),
       lty=c(3,1,3),
       col=c("darkgrey","black","darkgrey"),
       lwd=c(3,3,3),
       bty="n",xpd=NA)


mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Difference (years)",line=2)


dev.off()


ret_all %>% filter(age == 70)


svg("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_ccidci_alt.svg",width=10,height=4)
par(mfrow=c(1,2),mar=c(3,3,4,1))

# CCI
temp <- ret_all %>%
  filter(sex  %in% "man" & age %in% 18:100)



plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="Men",
     type="l",
     lwd=2,
     xlim=c(18,100),
     ylim=c(0,70),
     axes=FALSE)
axis(side=1,at=c(18,seq(20,100,by=10)))
axis(side=2,at=seq(0,70,by=10))

polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))

lines(temp$age,temp$le_05,lty=3,lwd=3,col="darkgrey")
lines(temp$age,temp$le_95,lty=3,lwd=3,col="darkgrey")


mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Life expectancy (years)",line=2)


legend(x=70,
       y=68,
       legend=c("Q 5%","Mean","Q 95%"),
       lty=c(3,1,3),
       col=c("darkgrey","black","darkgrey"),
       lwd=3,
       bty="n")


temp <- ret_all %>%
  filter(sex  %in% "woman" & age %in% 18:100)
plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="Women",
     type="l",
     lwd=2,
     xlim=c(18,100),
     ylim=c(0,70),
     axes=FALSE)
axis(side=1,at=c(18,seq(20,100,by=10)))
axis(side=2,at=seq(0,70,by=10))

polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))

lines(temp$age,temp$le_05,lty=3,lwd=3,col="darkgrey")
lines(temp$age,temp$le_95,lty=3,lwd=3,col="darkgrey")

mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Life expectancy (years)",line=2)

dev.off()





###
# cci and DCI
svg("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_ccidci_alt.svg",width=10,height=4)
par(mfrow=c(1,2),mar=c(3,3,4,1))

temp <- ret_all %>%
  filter(sex  %in% "man" & age %in% 18:100)

plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="Men",
     type="l",
     lwd=2,
     xlim=c(18,100),
     ylim=c(0,70),
     axes=FALSE)
axis(side=1,at=c(18,seq(20,100,by=10)))
axis(side=2,at=seq(0,70,by=10))

polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))

lines(temp$age,temp$le_05,lty=3,lwd=3,col="darkgrey")
lines(temp$age,temp$le_95,lty=3,lwd=3,col="darkgrey")


mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Life expectancy (years)",line=2)

legend(x=70,
       y=68,
       legend=c("Q 5%","Mean","Q 95%"),
       lty=c(3,1,3),
       col=c("darkgrey","black","darkgrey"),
       lwd=3,
       bty="n")

temp <- ret_all %>%
  filter(sex  %in% "woman" & age %in% 18:100 )
plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="Women",
     type="l",
     lwd=2,
     xlim=c(18,100),
     ylim=c(0,70),
     axes=FALSE)
axis(side=1,at=c(18,seq(20,100,by=10)))
axis(side=2,at=seq(0,70,by=10))

polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))

lines(temp$age,temp$le_05,lty=3,lwd=3,col="darkgrey")
lines(temp$age,temp$le_95,lty=3,lwd=3,col="darkgrey")


mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Life expectancy (years)",line=2)

dev.off()




###
# cci and DCI
svg("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_ccidci_alt2.svg",width=6,height=5)
par(mfrow=c(1,2),mar=c(3,3,4,1))

temp <- ret_all %>%
  filter(sex  %in% "man" & age %in% 18:100)

plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="Men",
     type="l",
     lwd=2,
     xlim=c(40,100),
     ylim=c(0,40),
     axes=FALSE)
axis(side=1,at=seq(40,100,by=10))
axis(side=2,at=seq(0,40,by=10))

polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))

lines(temp$age,temp$le_05,lty=3,lwd=3,col="darkgrey")
lines(temp$age,temp$le_95,lty=3,lwd=3,col="darkgrey")


mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Life expectancy (years)",line=2)

legend(x=55,
       y=40,
       legend=c("Q 5%","Mean","Q 95%"),
       lty=c(3,1,3),
       col=c("darkgrey","black","darkgrey"),
       lwd=3,
       bty="n")

temp <- ret_all %>%
  filter(sex  %in% "woman" & age %in% 18:100 )
plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="Women",
     type="l",
     lwd=2,
     xlim=c(40,100),
     ylim=c(0,40),
     axes=FALSE)
axis(side=1,at=seq(40,100,by=10))
axis(side=2,at=seq(0,40,by=10))

polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))

lines(temp$age,temp$le_05,lty=3,lwd=3,col="darkgrey")
lines(temp$age,temp$le_95,lty=3,lwd=3,col="darkgrey")



mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Life expectancy (years)",line=2)

dev.off()





###
# cci and DCI
svg("X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_ccidci_overall_alt.svg",width=10,height=4)
par(mfrow=c(1,2),mar=c(3,3,4,1))

temp <- ret_overall %>%
  filter(age %in% 18:100) 

plot(x=temp$age,
     y=temp$le_mean,
     xlab="",
     ylab="",
     main="",
     type="l",
     lwd=2,
     xlim=c(40,100),
     ylim=c(0,40),
     axes=FALSE)
axis(side=1,at=seq(40,100,by=10))
axis(side=2,at=seq(0,40,by=10))

polygon(x=c(temp$age,rev(temp$age)),
        y=c(temp$le_05,rev(temp$le_95)),
        border=NA,
        col=rgb(0,0,0,0.1))

lines(temp$age,temp$le_05,lty=3,lwd=2,col="darkgrey")
lines(temp$age,temp$le_95,lty=3,lwd=2,col="darkgrey")



mtext(side=1,text="Age (years)",line=2)
mtext(side=2,text="Life expectancy (years)",line=2)

legend(x=60,
       y=60,
       legend=c("Q 5%","Mean","Q 95%"),
       lty=c(3,1,3),
       col=c("darkgrey","black","darkgrey"),
       lwd=3,
       bty="n")



dev.off()





###
# End
###