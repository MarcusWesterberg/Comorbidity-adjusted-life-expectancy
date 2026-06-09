###
# Plot C-indices
require(tidyverse)
X <- "\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

log_path <- "C:\\Marcus\\Misc\\"

save_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")
results_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")



plot_cindex <- function(plot_data,
                        xlim=c(2003,2012),
                        ylim=c(0.5,1),
                        xlab="Year",
                        ylab="C-index",
                        cols=1:length(plot_data),
                        ages="70-79",
                        sexes="man, woman",
                        time_point=10,
                        xat=c(2003:2012),
                        yat=seq(0.5,1,by=0.1),
                        cex.labs=0.75,
                        cex.main=1,
                        add_at=NULL,
                        add_cex=1){
  
  plot(1,type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=ylim)
  axis(side=1,at=xat,labels=xat,cex.axis=cex.labs)
  axis(side=2,at=yat,labels=yat,cex.axis=cex.labs)
  mtext(side=1,text=xlab,line=2,cex=cex.labs)
  mtext(side=2,text=ylab,line=2,cex=cex.labs)
  
  for(i in 1:length(xat)){
    lines(x=c(xat[i],xat[i]),y=ylim,col="lightgrey",lty=2)
  }
  for(i in 1:length(yat)){
    lines(y=c(yat[i],yat[i]),x=xlim,col="lightgrey",lty=2)
  }
  
  sexes2 <- sexes
  if(grepl(sexes2,pattern=", ")){
    sexes2<-"all"
  } else if(sexes == "man"){
    sexes2 <- "men"
  } else{
    sexes2 <- "women"
  }
  
  #mtext(side=3,text=paste0("Time=",time_point,", age=",ages,", sex=",sexes2),line=0,cex=cex.main)
  pnames <- names(plot_data)
  
  for(j in 1:length(plot_data)){
    temp <- plot_data[[j]] %>%
      filter(indexage %in% ages & sex %in% sexes & timepoint %in% time_point) %>%
      arrange(indexyear)
    
    lines(x=temp$indexyear,
          y=temp$cindex_mean,
          col=cols[j],lwd=1)
    
    points(x=temp$indexyear,
           y=temp$cindex_mean,
           col=cols[j],
           pch=20,
           cex=2)
    
    # sp <- sqrt(temp$cindex_varpooled)
    # cil <- temp$cindex_mean-sp*qnorm(0.975)
    # ciu <- temp$cindex_mean+sp*qnorm(0.975)
    # if(!is.null(add_at)){
    #   text(y=add_at[j],
    #        x=mean(xlim),
    #        labels=round(temp$cindex_mean[length(temp$cindex_mean)],dig=3),
    #        col=cols[j],
    #        cex=add_cex)
    # }
    
  }
}




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
load(file=paste0(save_path,"cindices_all.Rdata"))
load(file=paste0(save_path,"subgroup_matrix_All.Rdata"))
subgroup_matrix$subgroup_index <- 1:nrow(subgroup_matrix)
cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% left_join(subgroup_matrix,by="subgroup_index")

plot_data <- list("CCI10"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10"),
                  "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI"),
                  "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI"))


ages <- sort(unique(plot_data[[1]]$indexage))
sexes <- c("man","woman")
cols <- c(rgb(0,0,0.9,0.75),rgb(0.75,0,0,0.75),rgb(0,0.5,0,0.75))




index <- 1
for(a in ages){
  svg(file=paste0(results_path,"cindex_time10_",a,".svg"),width=7,height=4)
  par(mfrow=c(1,2),mar=c(2,3,3,0.1),oma=c(1,1,1,1))
  for(sx in sexes){
    plot_cindex(plot_data=lapply(FUN=function(x) x %>% filter(indexyear %in% 2003:2013),plot_data ),
                xlim=c(2003,2013),
                ylim=c(0.5,1),
                xlab="Year",
                ylab="C-index",
                cols=cols,
                ages=a,
                sexes=sx,
                time_point=10,
                xat=c(2003:2013),
                yat=seq(0.5,1,by=0.1),
                cex.labs=0.75,
                cex.main=0.75,
                add_at=c(0.95,0.9,0.85),
                add_cex=1)
    if(index==1){
      #legend(x=2003,y=1,legend=c("CCI","DCI","MDCI"),col=cols,lty=c(1,1,1),lwd=2)
    }
    index <- index + 1
   
  }
  dev.off()
}





index <- 1
for(a in ages){
  svg(file=paste0(results_path,"cindex_time5_",a,".svg"),width=7,height=4)
  par(mfrow=c(1,2),mar=c(2,3,3,0.1),oma=c(1,1,1,1))
  for(sx in sexes){
    plot_cindex(lapply(FUN=function(x) x %>% filter(indexyear %in% 2003:2018),plot_data ),
                xlim=c(2003,2018),
                ylim=c(0.5,1),
                xlab="Year",
                ylab="C-index",
                cols=cols,
                ages=a,
                sexes=sx,
                time_point=5,
                xat=c(2003:2018),
                yat=seq(0.5,1,by=0.1),
                cex.labs=0.75,
                cex.main=0.75,
                add_at=c(0.95,0.9,0.85),
                add_cex=1)
    if(index==1){
      #legend(x=2003,y=1,legend=c("CCI","DCI","MDCI"),col=cols,lty=c(1,1,1),lwd=2)
    }
    index <- index + 1
  }
  dev.off()
}






index <- 1
for(a in ages){
  svg(file=paste0(results_path,"cindex_time1_",a,".svg"),width=7,height=4)
  par(mfrow=c(1,2),mar=c(2,3,3,0.1),oma=c(1,1,1,1))
  for(sx in sexes){
    plot_cindex(lapply(FUN=function(x) x %>% filter(indexyear %in% 2003:2022),plot_data ),
                xlim=c(2003,2022),
                ylim=c(0.5,1),
                xlab="Year",
                ylab="C-index",
                cols=cols,
                ages=a,
                sexes=sx,
                time_point=1,
                xat=c(2003:2022),
                yat=seq(0.5,1,by=0.1),
                cex.labs=0.75,
                cex.main=0.75,
                add_at=c(0.95,0.9,0.85),
                add_cex=1)
    if(index==1){
      # legend(x=2003,y=1,legend=c("CCI","DCI","MDCI"),col=cols,lty=c(1,1,1),lwd=2)
    }
    index <- index + 1
  }
  dev.off()
}











### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
load(file=paste0(save_path,"cindices_CCI1p.Rdata"))
load(file=paste0(save_path,"subgroup_matrix_CCI1p.Rdata"))
subgroup_matrix$subgroup_index <- 1:nrow(subgroup_matrix)
cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% left_join(subgroup_matrix,by="subgroup_index")

plot_data <- list("CCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10"),
                  "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI"),
                  "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI"))



ages <- sort(unique(plot_data[[1]]$indexage))
sexes <- c("man","woman")
cols <- c(rgb(0,0,0.9,0.75),rgb(0.75,0,0,0.75),rgb(0,0.5,0,0.75))



index <- 1
for(a in ages){
  svg(file=paste0(results_path,"cindex_CCI1p_time10_",a,".svg"),width=7,height=4)
  par(mfrow=c(1,2),mar=c(2,3,3,0.1),oma=c(1,1,1,1))
  for(sx in sexes){
    plot_cindex(plot_data=lapply(FUN=function(x) x %>% filter(indexyear %in% 2003:2013),plot_data ),
                xlim=c(2003,2013),
                ylim=c(0.5,1),
                xlab="Year",
                ylab="C-index",
                cols=cols,
                ages=a,
                sexes=sx,
                time_point=10,
                xat=c(2003:2013),
                yat=seq(0.5,1,by=0.1),
                cex.labs=0.75,
                cex.main=0.75,
                add_at=c(0.95,0.9,0.85),
                add_cex=1)
    if(index==1){
      #legend(x=2003,y=1,legend=c("CCI","DCI","MDCI"),col=cols,lty=c(1,1,1),lwd=2)
    }
    index <- index + 1
    
  }
  dev.off()
}





index <- 1
for(a in ages){
  svg(file=paste0(results_path,"cindex_CCI1p_time5_",a,".svg"),width=7,height=4)
  par(mfrow=c(1,2),mar=c(2,3,3,0.1),oma=c(1,1,1,1))
  for(sx in sexes){
    plot_cindex(lapply(FUN=function(x) x %>% filter(indexyear %in% 2003:2018),plot_data ),
                xlim=c(2003,2018),
                ylim=c(0.5,1),
                xlab="Year",
                ylab="C-index",
                cols=cols,
                ages=a,
                sexes=sx,
                time_point=5,
                xat=c(2003:2018),
                yat=seq(0.5,1,by=0.1),
                cex.labs=0.75,
                cex.main=0.75,
                add_at=c(0.95,0.9,0.85),
                add_cex=1)
    if(index==1){
      #legend(x=2003,y=1,legend=c("CCI","DCI","MDCI"),col=cols,lty=c(1,1,1),lwd=2)
    }
    index <- index + 1
  }
  dev.off()
}






index <- 1
for(a in ages){
  svg(file=paste0(results_path,"cindex_CCI1p_time1_",a,".svg"),width=7,height=4)
  par(mfrow=c(1,2),mar=c(2,3,3,0.1),oma=c(1,1,1,1))
  for(sx in sexes){
    plot_cindex(lapply(FUN=function(x) x %>% filter(indexyear %in% 2003:2022),plot_data ),
                xlim=c(2003,2022),
                ylim=c(0.5,1),
                xlab="Year",
                ylab="C-index",
                cols=cols,
                ages=a,
                sexes=sx,
                time_point=1,
                xat=c(2003:2022),
                yat=seq(0.5,1,by=0.1),
                cex.labs=0.75,
                cex.main=0.75,
                add_at=c(0.95,0.9,0.85),
                add_cex=1)
    if(index==1){
      # legend(x=2003,y=1,legend=c("CCI","DCI","MDCI"),col=cols,lty=c(1,1,1),lwd=2)
    }
    index <- index + 1
  }
  dev.off()
}















### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
load(file=paste0(save_path,"cindices_all_ages.Rdata"))
loads <- load(file=paste0(save_path,"subgroup_matrix_all_ages.Rdata"))
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% left_join(subgroups_matrix,by="subgroup_index")

iys <- unique(cindices$point_ests_aggr$indexyear)

plot_data10 <- list("CCI10"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10" & timepoint == 10 & indexyear %in% iys[1]) ,
                    "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI" & timepoint == 10 & indexyear %in% iys[1]),
                    "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI" & timepoint == 10 & indexyear %in% iys[1]))

plot_data5 <- list("CCI10"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10" & timepoint == 5& indexyear %in% iys[1] ),
                    "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI" & timepoint == 5 & indexyear %in% iys[1]),
                    "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI" & timepoint == 5 & indexyear %in% iys[1]))

plot_data1 <- list("CCI10"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10" & timepoint == 5 & indexyear %in% iys[1]),
                    "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI" & timepoint == 5 & indexyear %in% iys[1] ),
                    "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI"& timepoint == 5 & indexyear %in% iys[1]))

ages <- sort(unique(plot_data[[1]]$indexage))
sexes <- c("man","woman")
cols <- c(rgb(0,0,0.9,0.75),rgb(0.75,0,0,0.75),rgb(0,0.5,0,0.75))







plot_data <- plot_data1

ylim <- as.numeric(c(min(ages),max(ages)))
xlim <- c(0.5,1)
yat <- seq(ylim[1],ylim[2],1)
xat <- seq(0.5,1,0.1)
index <- 1
cex.labs <- 0.75

xlab <- "C-index"
svg(file=paste0(results_path,"cindex_time1_all_ages.svg"),width=6,height=6)
par(mfrow=c(1,2),mar=c(2,3,3,0.1),oma=c(1,1,1,1))
for(sx in sexes){
  plot(1,type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=ylim)
  axis(side=1,at=xat,labels=xat,cex.axis=cex.labs)
  #axis(side=2,at=yat,labels=yat,cex.axis=cex.labs)
  text(x=0.5,y=yat,labels=yat,pos=2,cex=0.5,xpd=NA)
  mtext(side=1,text=xlab,line=2,cex=cex.labs)

  
  for(i in 1:length(yat)){
    if(i %% 2 == 0){
      polygon(y=c(yat[i]-0.5,yat[i]-0.5,
                  yat[i]+0.5,yat[i]+0.5),x=c(xlim,rev(xlim)),col="lightgrey",border =NA)
    }
  }
 
  
  for(i in 1:length(xat)){
    lines(x=c(xat[i],xat[i]),
          y=ylim,
          col="grey")
    
  }
  
  #mtext(side=3,text=paste0("Time=",time_point,", age=",ages,", sex=",sexes2),line=0,cex=cex.main)
  pnames <- names(plot_data)
  
  for(a in ages){
    for(j in 1:length(plot_data)){
      temp <- plot_data[[j]] %>%
        filter(indexage %in% a & sex %in% sx)
      
      sp <- sqrt(temp$cindex_varpooled)
      cil <- temp$cindex_mean-sp*qnorm(0.975)
      ciu <- temp$cindex_mean+sp*qnorm(0.975)
  
      lines(y=as.numeric(c(a,a)),
            x=c(cil,ciu),
            col=cols[j])
    
      points(x=temp$cindex_mean,
             y=as.numeric(a),
             col=cols[j],
             pch=20,
             cex=0.25)
    }
  }
}
dev.off()



###
###
# End
###