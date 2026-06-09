loads <- load(file="life_expectancy_mdci_dci_models.Rdata") 

# Your dataset is called DD and contains the columns age, MDCI and dci and sex (man/woman)

require(mgcv)

DD$age <- floor(DD$age)

MDCI_cutoffs <- c(-1,4)
DCI_cutoffs <- c(-0.5,8)

trunc_it <- function(x,l,u){
  x[x<l]<-l
  x[x>u]<-u
  return(x)
}

DD$MDCI <- trunc_it(DD$MDCI,
                    l=MDCI_cutoffs[1],
                    u=MDCI_cutoffs[2])

DD$dci <- trunc_it(DD$dci,
                   l=DCI_cutoffs[1],
                   u=DCI_cutoffs[2])



DD$MDCI <- round(DD$MDCI,dig=1)
DD$dci <- round(DD$dci,dig=1)

DD_men <- DD %>%
  filter(sex == "man") 
DD_men$le <- exp(predict(le.models$men, newdata = DD_men, type="response"))

DD_women <- DD %>%
  filter(sex == "woman") 
DD_women$le <- exp(predict(le.models$women, newdata = DD_women, type="response"))

LE <- rbind(DD_men,DD_women)






loads <- load(file="X:\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\Results\\life_expectancy_cci_dci_models.Rdata") 

# Your dataset is called DD and contains the columns age, MDCI and dci and sex (man/woman)

require(mgcv)




DD$age <- floor(DD$age)

DCI_cutoffs <- c(-0.5,8)

trunc_it <- function(x,l,u){
  x[x<l]<-l
  x[x>u]<-u
  return(x)
}

DD$dci <- trunc_it(DD$dci,
                   l=DCI_cutoffs[1],
                   u=DCI_cutoffs[2])


DD$CCI <- as.numeric(DD$CCI)
DD$CCI[DD$CCI>6]<-6
DD$dci <- round(DD$dci,dig=1)

DD_men <- DD %>%
  filter(sex == "man") 
DD_men$le <- exp(predict(le.models$men, newdata = DD_men, type="response"))

DD_women <- DD %>%
  filter(sex == "woman") 
DD_women$le <- exp(predict(le.models$women, newdata = DD_women, type="response"))

LE <- rbind(DD_men,DD_women)


