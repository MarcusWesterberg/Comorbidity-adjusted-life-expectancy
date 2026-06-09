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
                        xat2=c(2003:2012),
                        yat=seq(0.5,1,by=0.1),
                        yat2=ylim[1]-0.05,
                        cex.labs=0.75,
                        cex.main=1,
                        add_at=NULL,
                        add_cex=1,
                        pchs){
  
  plot(1,type="n",
       axes=FALSE,
       xlab="",
       ylab="",
       xlim=xlim,
       ylim=ylim)
  axis(side=1,at=xat,labels=rep("",length(xat)),cex.axis=cex.labs,line=-0.5,padj=-1)
  
  text(y=yat2,x=xat2-0.1,labels=xat2,pos=1,cex=0.75,xpd=NA) # ,srt=45
  
  print(yat)
  axis(side=2,at=yat,labels=yat,cex.axis=cex.labs,line=-0.4,padj=0.75)
  
  mtext(side=1,text=xlab,line=0.55,cex=cex.labs)
  mtext(side=2,text=ylab,line=1.25,cex=cex.labs)
  

  
  for(i in 2:length(yat)){
    lines(y=c(yat[i],yat[i]),x=xlim+c(-0.05,-0.05),col="lightgrey",lty=1)
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
    
    if(temp$risk_score_name[1] %in%"DCI"){
      temp <- temp %>% filter(indexyear>2005)
    }
    
    lines(x=temp$indexyear,
          y=temp$cindex_mean,
          col=cols[j],lwd=1)
    
    points(x=temp$indexyear,
           y=temp$cindex_mean,
           col=cols[j],
           pch=pchs[[j]],
           cex=1)
    
    sp <- sqrt(temp$cindex_varpooled)
    cil <- temp$cindex_mean-sp*qnorm(0.975)
    ciu <- temp$cindex_mean+sp*qnorm(0.975)
    
    for(i in 1:length(temp$indexyear)){
      lines(x=c(temp$indexyear[i],temp$indexyear[i]),
            y=c(cil[i],ciu[i]),
            col=cols[j])
    }

    
  }
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
load(file=paste0(save_path,"cindices_all1_Elixhauser.Rdata"))
load(file=paste0(save_path,"subgroup_matrix_All1_Elixhauser.Rdata"))
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% 
  left_join(subgroups_matrix,by="subgroup_index")

cindices_Elixhauser <- cindices 

load(file=paste0(save_path,"cindices_all1_Elixhauser2.Rdata"))
load(file=paste0(save_path,"subgroup_matrix_All1_Elixhauser2.Rdata"))
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% 
  left_join(subgroups_matrix,by="subgroup_index")

cindices_Elixhauser2 <- cindices 


load(file=paste0(save_path,"cindices_all1.Rdata"))
load(file=paste0(save_path,"subgroup_matrix_All1.Rdata"))
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% 
  left_join(subgroups_matrix,by="subgroup_index")

plot_data <- list("CCI10"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10"),
                  "ECI"=cindices_Elixhauser$point_ests_aggr ,
                  "ECI2"=cindices_Elixhauser2$point_ests_aggr ,
                  "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI"),
                  "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI"))


ages <- sort(unique(plot_data[[1]]$indexage))


cols <- c(rgb(0,0,0,0.75),rgb(0.55,0.25,0.25,0.75),rgb(0.35,0.35,0.35,0.75),rgb(0.65,0.45,0,0.75),rgb(0,0,0.75,0.75))
pchs <- rep(18,20)
pchs <- list(pchs)

pchs <- rep(pchs,5)
pchs[[1]][1:5]<-1
pchs[[2]][1:5]<-1
pchs[[3]][1:5]<-1
pchs[[5]][1:5]<-1

sexes <- c("man, woman")
index <- 1
a <- ages[1]
for(a in ages){
  svg(file=paste0(results_path,"cindex_years_time1_all_",a,".svg"),width=5,height=4)
  par(mfrow=c(1,1),mar=c(1,2,1,0.1),oma=0.5*c(1,1,1,1))
  for(sx in sexes){
    plot_cindex(plot_data=lapply(FUN=function(x) x %>% filter(indexyear %in% 2006:2022),plot_data ),
                xlim=c(2006,2022.05),
                ylim=c(0.5,0.9),
                xlab="",
                ylab="C-index",
                cols=cols,
                ages=a,
                sexes=sx,
                time_point=1,
                xat=c(2006:2022),
                xat2=c(2006,2008,2010,2012,2014,2016,2018,2020,2022),
                yat=seq(0.5,0.9,by=0.05),
                yat2=0.485,
                cex.labs=0.75,
                cex.main=0.75,
                add_at=c(0.95,0.9,0.85),
                add_cex=1,
                pchs=pchs)
    if(index==1){
      #legend(x=2003,y=1,legend=c("CCI","DCI","MDCI"),col=cols,lty=c(1,1,1),lwd=2)
    }
    index <- index + 1
    
  }
  dev.off()
}




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
loads <- load(file=paste0(save_path,"cindices_all5.Rdata"))
load(file=paste0(save_path,"subgroup_matrix_All5.Rdata"))
subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% 
  left_join(subgroups_matrix,by="subgroup_index")

plot_data <- list("CCI10"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10" & sex %in% "man, woman"),
                  "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI" & sex %in% "man, woman"),
                  "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI" & sex %in% "man, woman"))


ages <- sort(unique(plot_data[[1]]$indexage))

cols <- c(rgb(0,0,0,0.75),rgb(0.65,0.45,0,0.75),rgb(0,0,0.75,0.75))
pchs <- rep(18,20)
pchs <- list(pchs)

pchs <- rep(pchs,3)
pchs[[1]][1:5]<-1
pchs[[3]][1:5]<-1

sexes <- c("man, woman")
index <- 1
a <- ages[1]
for(a in ages){
  svg(file=paste0(results_path,"cindex_years_time5_all_",a,".svg"),width=5,height=4)
  par(mfrow=c(1,1),mar=c(1,2,1,0.1),oma=0.5*c(1,1,1,1))
  for(sx in sexes){
    plot_cindex(plot_data=lapply(FUN=function(x) x %>% filter(indexyear %in% 2003:2022),plot_data ),
                xlim=c(2003,2022.05),
                ylim=c(0.5,0.9),
                xlab="",
                ylab="C-index",
                cols=cols,
                ages=a,
                sexes=sx,
                time_point=5,
                xat=c(2003:2022),
                xat2=c(2003,2005,2007,2009,2011,2013,2015,2017,2019,2022),
                yat=seq(0.5,0.9,by=0.05),
                yat2=0.485,
                cex.labs=0.75,
                cex.main=0.75,
                add_at=c(0.95,0.9,0.85),
                add_cex=1,
                pchs=pchs)
    if(index==1){
      #legend(x=2003,y=1,legend=c("CCI","DCI","MDCI"),col=cols,lty=c(1,1,1),lwd=2)
    }
    index <- index + 1
    
  }
  dev.off()
}



















###
###
# End
###