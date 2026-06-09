###
# Load SCB Data and estimate hazard function according to sex and age 
require(tidyverse)


# population size 
P <- openxlsx::read.xlsx("\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\Folkmängd 2006-2024.xlsx")
P[1:nrow(P) %% 2 == 0,1] <- P[1:nrow(P) %% 2 == 1,1]

P[,1]<-gsub(P[,1],pattern=" år",replacement="")
P[,1]<-gsub(P[,1],pattern="[+]",replacement="")

P[,1] <- as.numeric(P[,1])

P[1:nrow(P) %% 2 == 1,2]<-"men"
P[1:nrow(P) %% 2 == 0,2]<-"women"

P <- P %>% pivot_longer(cols=starts_with("20"),
                        names_to="year",
                        values_to="N") %>%
  mutate(N=as.numeric(N),
         year=as.numeric(year))

colnames(P)<-c("age","sex","year","N")

save(P,file="\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\Popsize.Rdata")

# death risks
D <- openxlsx::read.xlsx("\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\Dödsrisker 2006-2024.xlsx")

D_men <- D[1:93,]
D_women <- D[94:186,]

D_men$sex <- "men"
D_women$sex <- "women"

D <- rbind(D_men,D_women)
D$age <- gsub(D$age,pattern=" år",replacement="")
D$age <- as.numeric(D$age)

D_LE <- D[,c(1,2,22:40)]
D <- D[,c(1:21)]

D_LE <- D_LE %>% pivot_longer(cols=starts_with("20"),
                              names_to="year",
                              values_to="LE") %>%
  mutate(LE=as.numeric(LE),
         year=as.numeric(year))

D <- D %>% pivot_longer(cols=starts_with("20"),
                        names_to="year",
                        values_to="risk")%>%
  mutate(risk=as.numeric(risk)/1000,
         year=as.numeric(year)) # promille

save(D,file="\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\risk.Rdata")
save(D_LE,file="\\\\lul-net.ad.lul.se\\dfsLUL\\System\\UAS\\RCC-Forskning\\ComorbidityBase\\Works\\2024\\2024 Assessment of comorbidity indices\\SCB\\LE.Rdata")



###
# End
###