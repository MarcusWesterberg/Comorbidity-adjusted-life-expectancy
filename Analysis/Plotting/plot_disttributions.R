###
# Plot distributions

X<-"\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\"

require(tidyverse)

log_path<-"C:\\Marcus\\Misc\\"

data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"
background_data_path<-"C:\\Marcus\\ComorbidityBase_cache\\"

save_path<-paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Cache\\")
results_path <- paste0(X,"ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


loads <- load(file=paste0(save_path,"n_icdchapters.Rdata")) #n_icdchapters
loads <- load(file=paste0(save_path,"n_cci.Rdata")) # n_cci
loads <- load(file=paste0(save_path,"n_mdci.Rdata"))  # n_mdci
loads <- load(file=paste0(save_path,"n_mdci_sv.Rdata"))  # n_mdci
loads <- load(file=paste0(save_path,"n_dci.Rdata")) #
loads <- load(file=paste0(save_path,"n_elixhauser.Rdata")) #


n_icdchapters <- n_icdchapters %>%
  filter(indexage %in% 18:120) %>% 
  mutate(indexage=ifelse(indexage>100,100,indexage)) %>%
  group_by(sex,indexage) %>%
  summarize(n=sum(n),
            n0=sum(n0),
            n1=sum(n1),
            n2=sum(n2),
            n3=sum(n3),
            n4=sum(n4),.groups="keep") %>%
  arrange(indexage,sex)

  
n_cci <- n_cci %>%
  filter(indexage %in% 18:120) %>% 
  mutate(indexage=ifelse(indexage>100,100,indexage)) %>%
  group_by(sex,indexage) %>%
  summarize(n=sum(n),
            n0=sum(n0),
            n1=sum(n1),
            n2=sum(n2),
            n3=sum(n3),
            n4=sum(n4),.groups="keep") %>%
  arrange(indexage,sex)


n_elixhauser <- n_elixhauser %>%
  filter(indexage %in% 18:120) %>% 
  mutate(indexage=ifelse(indexage>100,100,indexage)) %>%
  group_by(sex,indexage) %>%
  summarize(n=sum(n),
            n0=sum(n0),
            n1=sum(n1),
            n2=sum(n2),
            n3=sum(n3),
            n4=sum(n4),.groups="keep") %>%
  arrange(indexage,sex)

n_mdci <- n_mdci %>%
  filter(indexage %in% 18:120) %>% 
  mutate(indexage=ifelse(indexage>100,100,indexage)) %>%
  group_by(sex,indexage) %>%
  summarize(n=sum(n),
            n0=sum(n0),
            n1=sum(n1),
            n2=sum(n2),
            n3=sum(n3),
            n4=sum(n4),.groups="keep") %>%
  arrange(indexage,sex)

n_mdci_sv <- n_mdci_sv %>%
  filter(indexage %in% 18:120) %>% 
  mutate(indexage=ifelse(indexage>100,100,indexage)) %>%
  group_by(sex,indexage) %>%
  summarize(n=sum(n),
            n0=sum(n0),
            n1=sum(n1),
            n2=sum(n2),
            n3=sum(n3),
            n4=sum(n4),.groups="keep") %>%
  arrange(indexage,sex)

n_dci <- n_dci %>%
  filter(indexage %in% 18:120) %>% 
  mutate(indexage=ifelse(indexage>100,100,indexage)) %>%
  group_by(sex,indexage) %>%
  summarize(n=sum(n),
            n0=sum(n0),
            n1=sum(n1),
            n2=sum(n2),
            n3=sum(n3),
            n4=sum(n4),.groups="keep") %>%
  arrange(indexage,sex)

ages <- sort(unique(n_icdchapters$indexage))






plotf <- function(plotdata,cols){
  ylim <- c(18,100)
  xlim <- c(-1,1)
  
  yat <- c(18,seq(20,100,by=10))
  yat2 <- seq(ylim[1],ylim[2],0.85)
  yat_labs <-yat
  yat_labs[1]<-"18"
  yat_labs[length(yat_labs)]<-"\u2265100"
  
  xat <- seq(-1,1,0.2)
  index <- 1
  cex.labs <- 0.75
  xlab <- ""
  
  plot(1,type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  sexes <- c("man","woman")
  
  for(sx in 1:2){
    temp <- plotdata %>%
      filter(sex %in% sexes[sx])
    
    props <- temp %>%
      mutate(n0=n0/n,
             n1=n1/n,
             n2=n2/n,
             n3=n3/n,
             n4=n4/n)
    
    props <- as.matrix(props[,paste0("n",0:4)])
    props <- t(apply(FUN=cumsum,props,MARGIN=1))
    props <- cbind(0,props)
    
    for(a in 1:nrow(temp)){
      xdir <- ifelse(sx==1,-1,1)
      ydelta <- 0.5
      xtemp <- props[a,]
      
      for(j in 1:5){
        polygon(x=xdir*c(xtemp[c(j,j)],xtemp[c(j+1,j+1)] ),
                y=c(ages[a]-ydelta,ages[a]+ydelta,ages[a]+ydelta,ages[a]-ydelta),
                border=cols[j],
                col=cols[j])
      }
    }
    
  }
  
  lines(x=c(0,0),y=ylim + c(-1,0.5),lwd=2)
  axis(side=1,at=xat,labels=paste0(100*abs(xat),"%"),cex.axis=cex.labs,line=-0.25,padj=-1.5)
  axis(side=2,at=yat,labels=yat_labs,cex.axis=cex.labs,padj=0.75,line=-0.5)
  #text(x=0,y=yat,labels=yat_labs,pos=1,cex=cex.labs,xpd=NA)
  mtext(side=1,text=xlab,line=0.75,cex=cex.labs)
  mtext(side=2,text="Age (years)",line=1,cex=cex.labs)
}


# cols <- c(rgb(0,0,0.9,0.75),rgb(0.75,0,0,0.75),rgb(0,0.5,0,0.75))
cols <- c( colorRampPalette(c(rgb(0.7,0.95,0.75,alpha=1),"wheat"))(3)[2],
          colorRampPalette(c("wheat","tan2"))(3)[-3],
          colorRampPalette(c("tan2","firebrick1"))(2)[-1],
          rep("red4",1))

svg(file=paste0(results_path,"n_icd_chapters_ages.svg"),width=3.5,height=3.5)
par(mfrow=c(1,1),mar=c(1.5,1.25,0.01,0.01),oma=0.5*c(1,1,0.2,1))
plotf(plotdata=n_icdchapters,cols=cols)
dev.off()
  
svg(file=paste0(results_path,"n_cci_ages.svg"),width=3.5,height=3.5)
par(mfrow=c(1,1),mar=c(1.5,1.25,0.01,0.01),oma=0.5*c(1,1,0.2,1))
plotf(plotdata=n_cci,cols=cols)
dev.off()

svg(file=paste0(results_path,"n_elixhauser_ages.svg"),width=3.5,height=3.5)
par(mfrow=c(1,1),mar=c(1.5,1.25,0.01,0.01),oma=0.5*c(1,1,0.2,1))
plotf(plotdata=n_elixhauser,cols=cols)
dev.off()

svg(file=paste0(results_path,"n_mdci_ages.svg"),width=3.5,height=3.5)
par(mfrow=c(1,1),mar=c(1.5,1.25,0.01,0.01),oma=0.5*c(1,1,0.2,1))
plotf(plotdata=n_mdci,cols=cols)
dev.off()


svg(file=paste0(results_path,"n_mdci_sv_ages.svg"),width=3.5,height=3.75)
par(mfrow=c(1,1),mar=c(1.5,1.25,0.01,0.01),oma=0.5*c(1,1,0.2,1))
plotf(plotdata=n_mdci_sv,cols=cols)
dev.off()

svg(file=paste0(results_path,"n_dci_ages.svg"),width=3.5,height=3.5)
par(mfrow=c(1,1),mar=c(1.5,1.25,0.01,0.01),oma=0.5*c(1,1,0.2,1))
plotf(plotdata=n_dci,cols=cols)
dev.off()

###
# plot legend
svg(file=paste0(results_path,"n_distr_legend_ages.svg"),width=6,height=6)
par(mfrow=c(1,1),mar=c(3,3,1,0.1),oma=0.5*c(1,1,1,1))

plot(1,type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=ylim)

for(j in 1:5){
  ydelta <- 0.5
  polygon(x=2*c(c(j,j),c(j+1,j+1) )/6-1,
          y=50+c(-ydelta,ydelta,ydelta,-ydelta),
          border=NA,
          col=cols[j])
}   

dev.off()



###
# END
###