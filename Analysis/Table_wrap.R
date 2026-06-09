###
# table function



# helper
tablepart <- function(X,
                      varname="orig_diaage",
                      variable_name="",
                      continuous=TRUE,
                      cut_it=TRUE,
                      breaks=c(0,60,70,75,80,Inf),
                      labels=c("<60","60-69","70-74","75-80","80+"),
                      levels=NULL,
                      add_info=TRUE,
                      only_true=FALSE,
                      missing=FALSE,
                      percent=TRUE,
                      dig=2,
                      dig_percent=0,
                      keep_first_column=TRUE,
                      add_lab = ""){
  
  if(class(X)[1] %in% "data.frame"){
    y <- X[, varname ]
  } else{
    y <- X[, varname ][[1]]
  }

  n <- length(y)
  
  ret <- list()
  
  # variable name 
  temp <- cbind(variable_name,"","")
  colnames(temp)<-NULL
  ret[[1]] <- temp
  
  
  if(!percent){add_lab<-""}
  
  # if continuous variable
  if(continuous){
    
    # compute mean, max, etc..
    if(add_info){
      ai  <- c(paste0(round(mean(y,na.rm = TRUE),dig=dig)," (",round(sd(y,na.rm = TRUE),dig=dig),")"),"")
      q   <- round(quantile(y,p=c(0.5,0.25,0.75),na.rm=TRUE),dig=dig )
      ai  <- rbind(ai,c(paste0(q[1]," (",q[2],"-",q[3],")"),""))
      #ai  <- rbind(ai,c(round(min(y,na.rm=TRUE),dig=dig),round(max(y,na.rm=TRUE),dig=dig),""))
      ai  <- cbind(c("Mean (SD)","Median (Q1, Q3)"),ai)
      rownames(ai)<-NULL; colnames(ai)<-NULL
      ret[[length(ret)+1]] <- ai
    } 
    
    # compute the variables binned into the levels 
    if(cut_it){
      yc   <- cut(y,breaks=breaks,include.lowest = TRUE,right=FALSE)
      yct  <- table(yc)
      ycp  <- paste0("(",round(yct/n*100,dig=dig_percent),")")
      if(!percent){ycp<-""}
      temp <- cbind(labels=paste0(labels,add_lab),yct,ycp)
      #temp <- cbind(temp,"")
      rownames(temp)<-NULL; colnames(temp)<-NULL
      ret[[length(ret)+1]] <- temp
    }
    
    # if discrete variable
  } else{
    if(is.null(levels)){
      yt  <- table(y)
    } else{
      yt <- rep(0,length(levels))
      for(l in 1:length(levels)){
        yt[l] <- sum(y %in% levels[l])
      }
      names(yt) <- levels
    }
    
    if(is.null(labels)){labels <- names(yt)}
    
    if(only_true){yt <- yt["TRUE"]}
    
    yp  <- paste0("(",round(yt/n*100,dig=dig_percent),")")
    if(!percent){yp<-""}
    
    temp <- cbind(labels=paste0(labels,add_lab),yt,yp)
    #temp <- cbind(temp,"")
    rownames(temp)<-NULL; colnames(temp)<-NULL
    ret[[length(ret)+1]] <- temp
  }
  
  if(missing){
    mis <- sum(is.na(y))
    temp <- cbind("Missing",mis,paste0("(",round(100*mis/n,dig=dig_percent),")"))
    colnames(temp)<-NULL
    ret[[length(ret)+1]] <- temp
  }
  
  ret <- do.call(rbind,ret)
  if(!keep_first_column){
    ret <- ret[,-1]
  }
  
  return(ret)
}

# main function
table_wrap <- function(X=D,
                       varname="orig_diaage",
                       variable_name="Age at diagnosis, years",
                       C,
                       cn,
                       continuous=TRUE,
                       cut_it=TRUE,
                       breaks=c(0,60,70,75,80,Inf),
                       labels=c("<60","60-69","70-74","75-80","80+"),
                       levels=NULL,
                       add_info=TRUE,
                       only_true=FALSE,
                       missing=FALSE,
                       percent=TRUE,
                       dig=2,
                       dig_percent=0,
                       add_lab = ""){
  
  ret <- list()
  ret[[1]] <- tablepart(X=X[C[,1], ],
                        varname=varname,
                        variable_name=variable_name,
                        continuous=continuous,
                        cut_it=cut_it,
                        breaks=breaks,
                        labels=labels,
                        levels=levels,
                        add_info=add_info,
                        only_true=only_true,
                        missing=missing,
                        percent=percent,
                        dig=dig,
                        dig_percent=dig_percent,
                        keep_first_column=TRUE,
                        add_lab = add_lab)
  
  for(j in 2:ncol(C)){
    
    ret[[j]] <- tablepart(X=X[C[,j], ],
                          varname=varname,
                          variable_name=variable_name,
                          continuous=continuous,
                          cut_it=cut_it,
                          breaks=breaks,
                          labels=labels,
                          levels=levels,
                          add_info=add_info,
                          only_true=only_true,
                          missing=missing,
                          percent=percent,
                          dig=dig,
                          dig_percent=dig_percent,
                          keep_first_column=FALSE,
                          add_lab = add_lab)
  }
  
  ret <- do.call(cbind,ret)
  
  colnames(ret)[1:ncol(ret) %% 2==0] <- cn
  colnames(ret)[1:ncol(ret) %% 2==1] <- ""
  rownames(ret)<-NULL
  return(ret)
}

###
# end
###