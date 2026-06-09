X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"
background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")
results_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 





cols <- c(rgb(0,0,0.9,0.75),rgb(0.75,0,0,0.75))



fus <- c(1,5)
fu <- 1
for(fu in fus){
  
  loads <- load(file=paste0(save_path,"cindices_all_ages",1,".Rdata"))
  loads <- load(file=paste0(save_path,"subgroup_matrix_all_ages",1,".Rdata"))
  subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
  cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% 
    left_join(subgroups_matrix,by="subgroup_index") %>%
    filter(timepoint %in% fu)
  
  plot_data1 <-  list("CCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10") ,
                      "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI" ),
                      "MDCI_sv"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI_sv" ),
                      "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI"  ))
  
  
  for(k in 1:4){
    print(k)
    
    plot_data <- plot_data1[[k]]
    
    ages <- c("18-29",30:99,"100-120")
    
    cols <- c(rgb(0,0,0.9,0.75),rgb(0.75,0,0,0.75))
    
    ylim <- c(29,100)
    xlim <- c(0.5,0.9)
    #yat <- seq(ylim[1],ylim[2],5)
    yat <- c(29,seq(40,100,by=10))
    yat2 <- seq(ylim[1],ylim[2],0.85)
    yat_labs <-yat
    yat_labs[1]<-"<30"
    yat_labs[length(yat_labs)]<-"\u2265100"
    yat2 <- 0.485
    
    xat <- seq(0.5,0.9,0.05)
    index <- 1
    cex.labs <- 0.75
    xlab <- "C-index"
    
    
    svg(file=paste0(results_path,"cindex_time",fu,"_",names(plot_data1)[k],"_ages.svg"),width=3.5,height=4)
    par(mfrow=c(1,1),mar=c(1,0.5,1,0.1),oma=0.5*c(1,1,1,1))
    
    plot(1,type="n",axes=FALSE,xlab="",ylab="",xlim=ylim,ylim=xlim)
    axis(side=2,at=xat,labels=xat,cex.axis=cex.labs,line=-0.40,padj=0.75)
    axis(side=1,at=yat,labels=rep("",length(yat)),cex.axis=cex.labs,line=-0.5)
    text(y=yat2,x=yat,labels=yat_labs,pos=1,cex=cex.labs,xpd=NA)
    mtext(side=2,text=xlab,line=1.1,cex=cex.labs)
    #mtext(side=1,text="Age (years)",line=0.25,cex=cex.labs)
    
    # for(i in 1:length(yat2)){
    #   if(i %% 2 == 1){
    #     polygon(x=c(yat2[i]-0.5,yat2[i]-0.5,
    #                 yat2[i]+0.5,yat2[i]+0.5),
    #             y=c(xlim,rev(xlim)),
    #             col=rgb(0.1,0.1,0.1,alpha=0.1),border =NA)
    #   }
    # }
    # 
    for(i in 2:length(xat)){
      lines(y=c(xat[i],xat[i]),
            x=ylim + c(-0.5,0.5),
            col="grey")
      
    }  
    
    a <- "18-29"
    for(a in ages){
      sexes <- c("man","woman")
      for(sx in 1:2){
        
        temp <- plot_data %>%
          filter(indexage %in% a & sex %in% sexes[sx])
        
        sp <- sqrt(temp$cindex_varpooled)
        cil <- temp$cindex_mean-sp*qnorm(0.975)
        ciu <- temp$cindex_mean+sp*qnorm(0.975)
        
        if(a %in% "100-120"){
          ac <- 100
        } else if(a %in% "18-29"){
          ac <- 29
        } else{
          ac <- as.numeric(a)
        }
        # print(c(sx,ac,temp$cindex_mean))
        
        lines(x=as.numeric(c(ac,ac)),
              y=c(cil,ciu),
              col=cols[sx])
        
        points(y=temp$cindex_mean,
               x=as.numeric(ac),
               col=cols[sx],
               pch=18,
               cex=0.75)
        
        lines(x=as.numeric(c(ac,ac))+c(-0.25,+0.25),
              y=c(temp$cindex_mean,temp$cindex_mean),
              col=cols[sx])
        
      }
    }
    
    dev.off()
  }
}










































sexes <- c( "man, woman","man","woman")

cols <- c(rgb(0,0,0,0.75),
          rgb(0.55,0.25,0.25,0.75),
          rgb(0.45,0.45,0.45,0.75),
          rgb(0.65,0.45,0,0.75),
          rgb(0,0,0.75,0.75))
### 
# Combined no sex
fus <- c(1)
fu <- 1

for(sx in 1:3){
  for(fu in fus){
    
    loads <- load(file=paste0(save_path,"cindices_all_ages",1,"_Elixhauser.Rdata"))
    loads <- load(file=paste0(save_path,"subgroup_matrix_all_ages",1,"_Elixhauser.Rdata"))
    subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
    
    cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% 
      left_join(subgroups_matrix,by="subgroup_index") %>%
      filter(timepoint %in% fu)  %>% 
      filter(sex %in% sexes[sx] )
    
    cindices_elixhauser <- cindices
    
    loads <- load(file=paste0(save_path,"cindices_all_ages",1,"_Elixhauser2.Rdata"))
    loads <- load(file=paste0(save_path,"subgroup_matrix_all_ages",1,"_Elixhauser2.Rdata"))
    subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)
    
    cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% 
      left_join(subgroups_matrix,by="subgroup_index") %>%
      filter(timepoint %in% fu)  %>% 
      filter(sex %in% sexes[sx] )
    
    cindices_elixhauser2 <- cindices
    
    loads <- load(file=paste0(save_path,"cindices_all_ages",1,".Rdata"))
    loads <- load(file=paste0(save_path,"subgroup_matrix_all_ages",1,".Rdata"))
    subgroups_matrix$subgroup_index <- 1:nrow(subgroups_matrix)

    cindices$point_ests_aggr <-  cindices$point_ests_aggr %>% 
      left_join(subgroups_matrix,by="subgroup_index") %>%
      filter(timepoint %in% fu)  %>% 
      filter(sex %in% sexes[sx] )

    plot_data1 <-  list("CCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "CCI10") ,
                        "ECI"=cindices_elixhauser$point_ests_aggr ,
                        "ECI2"=cindices_elixhauser2$point_ests_aggr ,
                        "DCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "DCI"  ),
                        "MDCI"=cindices$point_ests_aggr %>% filter(risk_score_name %in% "MDCI" ))
    
    ages <- c("18-29",30:99,"100-120")
    
    ylim <- c(29,100)
    xlim <- c(0.5,0.9)
    #yat <- seq(ylim[1],ylim[2],5)
    yat <- c(29,seq(40,100,by=10))
    yat2 <- seq(ylim[1],ylim[2],0.85)
    yat_labs <-yat
    yat_labs[1]<-"<30"
    yat_labs[length(yat_labs)]<-"\u2265100"
    yat2 <- 0.485
    
    xat <- seq(0.5,0.9,0.05)
    index <- 1
    cex.labs <- 0.75
    xlab <- "C-index"
    
    
    svg(file=paste0(results_path,"cindex_time",fu,"_","comb",sexes[sx],"_ages.svg"),width=3.5,height=4)
    par(mfrow=c(1,1),mar=c(1,0.5,1,0.1),oma=0.5*c(1,1,1,1))
    
    plot(1,type="n",axes=FALSE,xlab="",ylab="",xlim=ylim,ylim=xlim)
    axis(side=2,at=xat,labels=xat,cex.axis=cex.labs,line=-0.40,padj=0.75)
    axis(side=1,at=yat,labels=rep("",length(yat)),cex.axis=cex.labs,line=-0.5)
    text(y=yat2,x=yat,labels=yat_labs,pos=1,cex=cex.labs,xpd=NA)
    mtext(side=2,text=xlab,line=1.1,cex=cex.labs)
    #mtext(side=1,text="Age (years)",line=0.25,cex=cex.labs)
    
    # for(i in 1:length(yat2)){
    #   if(i %% 2 == 1){
    #     polygon(x=c(yat2[i]-0.5,yat2[i]-0.5,
    #                 yat2[i]+0.5,yat2[i]+0.5),
    #             y=c(xlim,rev(xlim)),
    #             col=rgb(0.1,0.1,0.1,alpha=0.1),border =NA)
    #   }
    # }
    # 
    for(i in 2:length(xat)){
      lines(y=c(xat[i],xat[i]),
            x=ylim + c(-0.5,0.5),
            col="grey")
      
    }  
    
    
    for(k in 1:5){
      print(k)
      
      plot_data <- plot_data1[[k]]
      
      a <- "18-29"
      for(a in ages){
        
        temp <- plot_data %>%
          filter(indexage %in% a )
        
        sp <- sqrt(temp$cindex_varpooled)
        cil <- temp$cindex_mean-sp*qnorm(0.975)
        ciu <- temp$cindex_mean+sp*qnorm(0.975)
        
        if(a %in% "100-120"){
          ac <- 100
        } else if(a %in% "18-29"){
          ac <- 29
        } else{
          ac <- as.numeric(a)
        }
        # print(c(sx,ac,temp$cindex_mean))
        
        lines(x=as.numeric(c(ac,ac)),
              y=c(cil,ciu),
              col=cols[k])
        
        points(y=temp$cindex_mean,
               x=as.numeric(ac),
               col=cols[k],
               pch=18,
               cex=0.75)
        
        lines(x=as.numeric(c(ac,ac))+c(-0.25,+0.25),
              y=c(temp$cindex_mean,temp$cindex_mean),
              col=cols[k])
      }
    }
    dev.off()
  }
}

