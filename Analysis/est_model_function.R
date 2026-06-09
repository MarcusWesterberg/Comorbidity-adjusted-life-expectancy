est_model <- function(data,
                      sexstrata,
                      formula,
                      maxfu,
                      name,
                      init_values){
  
  
  if(init_values[2]<0){
    init_values[2] <- exp(init_values[2])
  }
  
  set.seed(1)
  
  data <- data %>% filter(sex %in% sexstrata)
  
  data$censor[data$timefu>maxfu] <- 0
  data$timefu[data$timefu>maxfu] <- maxfu
  
  require(flexsurv) 
  
  if(grepl(name,pattern="MDCI")){
    inits <- c(init_values,0,0,0)
  } else if(grepl(name,pattern="CCI") & grepl(name,pattern="DCI")  ){
    inits <- c(init_values,0,0,0,0,0,0,0)
  } else if(grepl(name,pattern="CCI")){
    inits <- c(init_values,0,0,0)
  } else{
    inits <- init_values
  }
  
  if(inits[2]<0){
    inits[2]<-exp(inits[2])
  }
  
  
  init_list <- list()
  
  n_init_guesses <- 11
  n <- min(20000,floor(nrow(data)/n_init_guesses))
  sample_inds <- sample(1:nrow(data),size=n*n_init_guesses,replace = FALSE)
  sample_inds <- split(sample_inds,f=1:length(sample_inds) %% n_init_guesses)
  
  for(k in 1:n_init_guesses){
    print(paste0("Guess ",k))
    #print(inits)
    
    data.s <- data[sample_inds[[k]],]
    
    model.s <- tryCatch(
      {
        flexsurvreg(formula,
                    data = data.s,
                    dist = "gompertz",
                    hessian=FALSE,
                    inits=inits,
                    sr.control = list(rel.tolerance=1e-06),
                    control=list(fnscale = nrow(data.s),
                                 trace=FALSE,
                                 maxit=500))
      },
      error = function(cond) {
        message(conditionMessage(cond))
        NA
      }
    )
    
    if(!class(model.s) %in% "logical"){
      if(length(coefficients(model.s))==length(inits)){
        init_list[[k]] <- coefficients(model.s)
      }
      
    } else{
      init_list[[k]] <- inits
    }
  }
  
  init_list <- do.call(rbind,init_list)
  
  inits2 <- apply(FUN=quantile,init_list,MARGIN=2,probs=0.5)
  
  if(grepl(name,pattern="CCI")){
    inits_temp <-  inits2[-c(1:2)]
    inits2[-c(1:2)] <- inits_temp
  }
  
  if(inits2[2]<0){
    inits2[2] <- exp(inits2[2])
  }
  
  inits <- 0.5 * (inits + inits2)
  
  if(inits[2]<0){
    inits[2] <- exp(inits[2])
  }
 
  
  print("Third guess")
  print(inits )
  
  rm(data.s)
  rm(model.s)
  
  model.m <- tryCatch(
    {
      flexsurvreg(formula,
                  data = data,
                  dist = "gompertz",
                  hessian=FALSE,
                  inits=inits,
                  sr.control = list(rel.tolerance=1e-08),
                  control=list(fnscale = nrow(data),
                               trace=TRUE,
                               maxit=1000) )
    },
    error = function(cond) {
      message(conditionMessage(cond))
      NA
    }
  )
  
 
  if(class(model.m) %in% "logical"){
    stop("------------------- Error in model -------------------")
  } 
  
  print("Final estimates")
  print(coefficients(model.m))

  
  model.m$logliki <- NULL
  model.m$data$Y <- NULL
  model.m$data$m <- NULL
  model.m$data$mml$rate <- NULL
  model.m$data$mml$shape <- NULL
  
  rm(data)
  gc()
  
  for(j in 1:length(model.m)){
    attr(model.m[[names(model.m)[j] ]],".Environment") <- c()
  }
  
  attr(model.m,".Environment") <- c()
  
  model.m$maxfu <- maxfu
  model.m$sex <- sexstrata
  model.m$name <- name
  
  gc()
  
  return(model.m)
}
