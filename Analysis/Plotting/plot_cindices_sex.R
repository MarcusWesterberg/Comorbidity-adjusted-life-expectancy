###
# Plot C-indices
require(tidyverse)
X <- "\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

log_path <- "C:\\Marcus\\Misc\\"

save_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")
results_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")



plot_cindex <- function(pd,
                        ages,
                        ages_lables,
                        xlabs_y=-0.05,
                        xlim=c(1,8),
                        ylim=c(0.5,1),
                        ylab="C-index",
                        cols=c(1,2),
                        yat=seq(0.5,1,by=0.05),
                        cex.labs=0.75){
  
  plot(1,type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=ylim)
  #axis(side=1,at=xat,labels=xat,cex.axis=cex.labs)
  axis(side=2,at=yat,labels=yat,cex.axis=cex.labs)
  #mtext(side=1,text=xlab,line=2,cex=cex.labs)
  mtext(side=2,text=ylab,line=2,cex=cex.labs)
  
  for(i in 1:length(yat)){
    lines(y=c(yat[i],yat[i]),x=xlim,col="lightgrey",lty=2)
  }
  
  x <- seq(xlim[1],xlim[2],length.out=length(ages))
  
  text(x=x,
       y=rep(0.5 + xlabs_y,lengt.out=length(ages)),
       labels=ages_lables,
       xpd=NA,
       cex=cex.labs)
  
  #mtext(side=3,text=paste0("Time=",time_point,", age=",ages,", sex=",sexes2),line=0,cex=cex.main)

  sx <- c("man","woman")

  for(k in 1:2){
    temp <- pd %>%
      filter(sex %in% sx[k]  ) %>%
      arrange(indexyear)
    
    y <- temp$cindex_mean
    xtemp <- temp$indexage
    stopifnot(all(xtemp==ages))
   
    points(x=x,
           y=y,
           col=cols[k],
           pch=20,
           cex=2)
    
    
    sp <- sqrt(temp$cindex_varpooled)
    cil <- temp$cindex_mean-sp*qnorm(0.975)
    ciu <- temp$cindex_mean+sp*qnorm(0.975)
    
    for(j in 1:length(x)){
      lines(y=as.numeric(c(x[j],x[j])),
            x=c(cil[j],ciu[j]),
            col=cols[k])
    }
    
  }
}




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
load(file=paste0(save_path,"cindices_all1.Rdata"))
load(file=paste0(save_path,"subgroup_matrix_All1.Rdata"))
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% 
  left_join(subgroups_matrix,by="subgroup_index")


cindices$point_ests_aggr <- cindices$point_ests_aggr %>% 
  filter(indexyear %in% "2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022" & 
           sex %in% c("man","woman"))

plot_data <- list("CCI10"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10" ),
                  "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI"),
                  "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI"))

ages <- c("18-120","18-39","40-49" ,"50-59","60-69","70-79" , "80-89",  "90-100")
ages_lables <- ages
ages_lables[1] <- "All"
ages_lables[-1] <- paste0(ages_lables[-1],"\nyears")


cols_men <- c(rgb(0,0,0.5,0.75),rgb(0,0,0.5,0.75),rgb(0,0,0.5,0.75))
cols_women <- c(rgb(0.25,0.5,0,0.75),rgb(0.25,0.5,0,0.75),rgb(0.25,0.5,0,0.75))
cnames <- c("CCI","MDCI","DCI")

svg(file=paste0(results_path,"cindex_sex_time1.svg"),width=7,height=4)
par(mfrow=c(1,3),mar=c(2,3,3,0.1),oma=c(1,1,1,1))
for(j in 1:3){
  plot_cindex(pd = plot_data[[j]],
              ages=ages,
              ages_lables=ages_lables,
              xlim=c(1,8),
              ylim=c(0.5,1),
              ylab="C-index",
              cols=c(cols_men[j],cols_women[j]),
              yat=seq(0.5,1,by=0.05),
              cex.labs=0.75)
  mtext(side=3,line=0,cnames[j])
  
}
dev.off()
  










###
###
# End
###