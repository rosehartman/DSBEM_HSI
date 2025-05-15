#How much do you need to increase Pseudodiaptomus in order to see a difference in smelt?

#In this case I think I"ll use average temp and turbidity, use different zoop scenarios instead of years


#Version with constant Temperature and turbididy
# 1. Read model parameters and data
# 3. Model bioenergetics and reference points
# 4. Get summaries
# 5. Graphs
library(tidyverse)
library(abind)
library(zoo)
library(lubridate)
library(lattice)
library(ggplot2)
library(gridExtra)
library(paletteer)
library(latticeExtra)
library(gsw)
library(chillR)

### 1. Read model parameters and data
n.strata <- 7 #actually number of temperature scenarios
ex.strata = sort(c( "Tempave", "Temp_05", "Temp_1", "Temp_15", "Temp05", "Temp1", "Temp15"))
n.prey <- 12 #may adjust this later
n.days<-153 # summer and fall
ln.a<-c(0.000005,0.00000183) # length-weight parameters
ln.b<-c(3,3.38)
#startY <- 2010
#endY <- 2022 #expand through 2024 when all data are in
#yr.seq<- seq(startY-2010,endY-2010,by=1)

n.yrs<-11 #actually the number of zooplankton levels
stage =4

#I'm not going to change any of these
VB.k<-2.72 #2.98 # VBGM parameters estimated in DSM TN 35
Linf<-78.39 #76.1
t0<--0.026 #-0.014


start.L = c(rep(23, 13)) #just use an average start lengths. 
terminal.L <- c((start.L[1]+30.025),(start.L[2]+30.025),49.6,53.5,52.1,57) # EDSM observed, (10/15-11/15),2017,2018,2019
start.Wt <- ln.a[2]*start.L^ln.b[2]
terminal.Wt <- ln.a[2]*terminal.L^ln.b[2]


terminal.Wtse <- c(0.299,0.356,0.367)
mn.terminal.Wtse <- mean(terminal.Wtse)
terminal.Wtse <- c(mn.terminal.Wtse ,mn.terminal.Wtse ,terminal.Wtse,mn.terminal.Wtse)

first.day <- c(1,(n.days+1)) # Jun 1 index for 2018 and 2019

# bioenergetics model parameters, given by Rose et al. 2013a
a.c <- c(0.18,0.18,0.18,0.1,0.1) # Rose
b.c <- c(-0.275,-0.275,-0.275,-0.54,-0.54)
CQ <- c(7,7,10,10,10)						# Temperature at CK1 of maximum (deg C)
T.L <- c(28,28,27,27,27)					# Temperature at CK4 of maximum (deg C)
CK.1 <- c(0.4,0.4,0.4,0.4,0.4)				# effect at temperature CQ
CK.4 <- c(0.01,0.01,0.01,0.001,0.001)			# effect at temperature T.L
T.0 <- c(17,17,20,20,20) # Temperature at 0.98 of maximum (deg C)
T.M <- c(20,20,23,23,23) # Temperature at 0.98 of maximum (deg C)
# Metabolism (R) parameters
a.r <- c(0.0027,0.0027,0.0027,0.0027,0.0027)	# weight multiplier
b.r <- c(-0.216,-0.216,-0.216,-0.216,-0.216)	# weight exponent
R.Q <- c(0.036,0.036,0.036,0.036,0.036)			# exponent for temperature effect
S.d <- c(0.175,0.175,0.175,0.175,0.175)			# Fraction of assimilated food lost to SDA (specific dynamic action)

# Egestion (F) and excretion (U) parameters
F.a <- c(0.16,0.16,0.16,0.16,0.16)				# Fraction of consumed food lost to egestion
U.a <- c(0.1,0.1,0.1,0.1,0.1)				# Fraction of assimilated food lost to excretion

e.s <- 4814									# J/g: convert g(prey)/g(delta smelt) to g(smelt)/g(smelt) - fixed
e.d <- c(2590,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590)		# energy density of prey items

V <- matrix(NA,n.prey,5)					# V.ij is vulnerability of prey type j to fish i. Set to 1 for all life stages eating all zooplankton types; except DS larvae values were 0 except for Limnoithona
V[1,] <- c(1,1,1,1,0) #limno ***
V[2,] <- c(1,1,1,1,1) #othcaljuv ***
V[3,] <- c(1,1,1,1,1) #pdiapjuv
V[4,] <- c(0,0,1,1,1) #othcalad ***
V[5,] <- c(0,0,1,1,1) #acartela
V[6,] <- c(0,0,1,1,1) #othclad
V[7,] <- c(1,1,1,0,0) #allcopnaup
V[8,] <- c(0,0,0,1,1) #daphnia
V[9,] <- c(0,1,1,1,1) #othcyc ***
V[10,] <- c(0,0,1,0,1) #other
V[11,] <- c(0,0,1,1,1) #eurytem ***
V[12,] <- c(0,0,1,1,1) #pdiapfor ***

K <- matrix(NA,n.prey,5)					# K.ik is half-saturation constant for fish i feeding on each prey type k - calculated outside model to obtain realistic diet and consumption rates
K[1,] <- c(NA,2.5,120,1.5,100) #limno				# From 2nd Rose model_apr2020
K[2,] <- c(NA,0.375,0.24,7.5,3) #othcaljuv
K[3,] <- c(NA,0.375,0.24,1.5,2) #pdiapjuv
K[4,] <- c(NA,250,6,0.75,0.6) #othcalad
K[5,] <- c(NA,250,36,0.75,0.25) #acartela
K[6,] <- c(NA,250,120,4.5,1) #othclad
K[7,] <- c(7.5,7.5,120,75,100) #allcopnaup
K[8,] <- c(NA,250,200,4.5,0.15) #daphnia
K[9,] <- c(NA,1.5,1.2,1.5,2) #othcyc
K[10,] <- c(NA,250,12,7.5,3) #other
K[11,] <- c(NA,250,6,0.375,0.25) #eurytem
K[12,] <- c(NA,250,2.4,0.375,0.25) #pdiapfor
for (i in 1:5) {
  K[,i] <- K[,i]
}

# parameters for new dbl logit temp fx 
logit.temp.b<-c(0.43,10.94)
logit.temp.c<-c(-2.78,23.62)

# turbidity effect on Cmax
# based on Hasenbein et al. (2016), Fig. 2
# Model turbidity effect on Cmax and natural mortality (M)
# based on Hasenbein et al. (2016), Fig. 2
max.turb <- 5 #248.51*5^-0.674 # min NTU measured=5, convert to Secchi
mid1.turb <- 35 #248.51*25^-0.674 # lower NTU at max feeding rate=25, convert to Secchi
min1 <- 0.65 #1.5*8/25 # feeding rate at minNTU=8, feeding rate at 35NTU=25
min2 <- 0.65 # lowest Cmax effect at low turbidity
turb.L.fx.mod <- glm(c(min1,min2)~c(20,45))
a.turb <- (max.turb+mid1.turb)/2 # glm logit regression parameters
b.turb <- 0.12


terminal.Wtse <- c(0.299,0.356,0.367)
mn.terminal.Wtse <- mean(terminal.Wtse)
terminal.Wtse <- c(mn.terminal.Wtse ,mn.terminal.Wtse ,terminal.Wtse,mn.terminal.Wtse)

first.day <- c(1,(n.days+1)) # Jun 1 index for 2018 and 2019



####Sensativity analysis ##########################
#See if I can change psdudos in each year, but hold other taxa the same. 

load("zoopsmwide.RData")

zoopsmwidef = filter(zoopsmAve, doy %in% c(153:305))%>%
  arrange(Region)

#average biomass of pseudodiaptomus by month in Suisun Bay
zoopsSuisun = filter(zoopsmwidef, Region %in% c("NE Suisun", "NW Suisun","SE Suisun")) %>%
  group_by(doy) %>%
  summarize(across(limno:pdiapfor, function(x) mean(x, na.rm = TRUE))) %>%
  mutate(level =1)
zoopsS = zoopsSuisun
#now replicate this for each biomass level
Levels = c(0.01, 0.5, 1.5, 2, 5, 10, 20, 40, 100, 200,500)
for(i in 1:10){
  Level = Levels[i]
  Pseudo = zoopsSuisun %>%
    mutate(level = Level, pdiapfor = level*pdiapfor, pdiapjuv = level*pdiapjuv)
  zoopsS = bind_rows(zoopsS, Pseudo)
}

pseudosave = ungroup(zoopsS) %>%
  group_by(level) %>%
  summarize(pdiapfor = mean(pdiapfor, na.rm =T))
yr.seq = pseudosave$level
#now the different temperature regimes

zoopsSS = bind_rows(mutate(zoopsS, Temp = "Temp_05"),
                    mutate(zoopsS, Temp = "Tempave"),
                    mutate(zoopsS, Temp = "Temp1"),
                    mutate(zoopsS, Temp = "Temp15"),
                    mutate(zoopsS, Temp = "Temp05"),
                    mutate(zoopsS, Temp = "Temp_1"),
                    mutate(zoopsS, Temp = "Temp_15"))  %>%
  select(Temp, level,doy, 
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) 

test_ps= zoopsSS%>%
  split(list(zoopsSS$Temp, zoopsSS$level))

zoop_ps= array(unlist(test_ps),dim=c(153,15,7, 11), 
               dimnames = list(c(153:305), names(zoopsSS), unique(zoopsSS$Temp),
                               unique(zoopsSS$level)))


#OK! Now i just need to get rid of the id columns
zoopx_ps = zoop_ps[(1:153), c(4:15), c(1:7), c(1:11)]
zoopx2_ps = apply(zoopx_ps, c(2,3,4), as.numeric)
PD.mn.array= zoopx2_ps #days by prey by strata by year!


### 1. Model temperature effect on Cmax
stage<-4
temp.fx <- function(new.temp.dat,stage,T.M1) {
  L.1<- exp((1/(T.0[stage]-CQ[stage]))*log(0.98*(1-CK.1[stage])/(CK.1[stage]*0.02))*(new.temp.dat-CQ[stage]))
  L.2 <- exp((1/(T.L[stage]-T.M1))*log(0.98*(1-CK.4[stage])/(CK.4[stage]*0.02))*(T.L[stage]-new.temp.dat))
  K.A <- CK.1[stage]*L.1/(1+CK.1[stage]*(L.1-1))
  K.B <- CK.4[stage]*L.2/(1+CK.4[stage]*(L.2-1))
  f.x <- K.A*K.B
  
  return(f.x)
}

### 2.

#Different temperature regimes
load("WaterQuality20102022.RData")

Tempave = mean(filter(AllWQmean2, Parameter == "watertemperature", DOY %in% c(153:305))$ValueImputed)

Tempwide_x = filter(AllWQmean2, Parameter == "watertemperature")%>%
  select(Parameter, DOY, Region) %>%
  distinct() %>%
  mutate(ValueImputed = case_when(Region == "Lower Sacramento River" ~ Tempave,
                                  Region == "Confluence" ~ Tempave-.5,
                                  Region == "Suisun Marsh" ~ Tempave-1,
                                  Region == "SW Suisun" ~ Tempave-1.5,
                                  Region == "SE Suisun" ~ Tempave+.5,
                                  Region == "NE Suisun" ~ Tempave+1,
                                  Region == "NW Suisun" ~ Tempave+1.5),
         Temp = case_when(Region == "Lower Sacramento River" ~ "Tempave",
                          Region == "Confluence" ~ "Temp_05",
                          Region == "Suisun Marsh" ~ "Temp_1",
                          Region == "SW Suisun" ~ "Temp_15",
                          Region == "SE Suisun" ~ "Temp05",
                          Region == "NE Suisun" ~ "Temp1",
                          Region == "NW Suisun" ~ "Temp15")) %>%
  arrange(Temp) %>%
   ungroup() %>%
  select(-Region) %>%
  pivot_wider(names_from = Temp, values_from = ValueImputed) 
 


Tempwidef_x = filter(Tempwide_x, DOY %in% c(153:305))
Tempwidef_x2 = bind_rows(mutate(Tempwidef_x, level = pseudosave$level[1]),
                         mutate(Tempwidef_x, level = pseudosave$level[2]),
                         mutate(Tempwidef_x, level = pseudosave$level[3]),
                         mutate(Tempwidef_x, level = pseudosave$level[4]),
                         mutate(Tempwidef_x, level = pseudosave$level[5]),
                         mutate(Tempwidef_x, level = pseudosave$level[6]),
                         mutate(Tempwidef_x, level = pseudosave$level[7]),
                         mutate(Tempwidef_x, level = pseudosave$level[8]),
                         mutate(Tempwidef_x, level = pseudosave$level[9]),
                         mutate(Tempwidef_x, level = pseudosave$level[10]),
                         mutate(Tempwidef_x, level = pseudosave$level[11]))

test_x = Tempwidef_x2  %>%
  split(list(Tempwidef_x2$level))

Temp_x = array(unlist(test_x),dim=c(153,10, 11), 
                      dimnames = list(c(153:305), names(Tempwidef_x2),
                                      unique(Tempwidef_x2$level)))

#OK! Now i just need to get rid of the id columns
Tempx_x = Temp_x[c(1:153), c(3:9), c(1:11)]
Tempx2_x = apply(Tempx_x, c(2,3), as.numeric)



#now turbidity (constant)


Turbave = mean(filter(AllWQmean2, Parameter == "turbidity", DOY %in% c(153:305))$ValueImputed)

turbwide_x = filter(AllWQmean2, Parameter == "turbidity")%>%
  select(Parameter, DOY, Region) %>%
  distinct() %>%
  mutate(ValueImputed = Turbave,
         Temp = case_when(Region == "Lower Sacramento River" ~ "Tempave",
                          Region == "Confluence" ~ "Temp_05",
                          Region == "Suisun Marsh" ~ "Temp_1",
                          Region == "SW Suisun" ~ "Temp_15",
                          Region == "SE Suisun" ~ "Temp05",
                          Region == "NE Suisun" ~ "Temp1",
                          Region == "NW Suisun" ~ "Temp15")) %>%
  arrange(Temp) %>%
  ungroup() %>%
  select(-Region, -Parameter) %>%
  pivot_wider(names_from = Temp, values_from = ValueImputed) 


Turbwidef_x = filter(turbwide_x, DOY %in% c(153:305))
Turbwidef_x2 = bind_rows(mutate(Turbwidef_x, level = pseudosave$level[1]),
                         mutate(Turbwidef_x, level = pseudosave$level[2]),
                         mutate(Turbwidef_x, level = pseudosave$level[3]),
                         mutate(Turbwidef_x, level = pseudosave$level[4]),
                         mutate(Turbwidef_x, level = pseudosave$level[5]),
                         mutate(Turbwidef_x, level = pseudosave$level[6]),
                         mutate(Turbwidef_x, level = pseudosave$level[7]),
                         mutate(Turbwidef_x, level = pseudosave$level[8]),
                         mutate(Turbwidef_x, level = pseudosave$level[9]),
                         mutate(Turbwidef_x, level = pseudosave$level[10]),
                         mutate(Turbwidef_x, level = pseudosave$level[11]),)


test_xt = Turbwidef_x2  %>%
  split(list(Turbwidef_x2$level))

Turb_x = array(unlist(test_xt),dim=c(153,9, 11), 
               dimnames = list(c(153:305), names(Turbwidef_x2),
                               unique(Turbwidef_x2$level)))

#OK! Now i just need to get rid of the id columns
turbx_x = Turb_x[c(1:153), c(2:8), c(1:11)]
turbx2_x = apply(turbx_x, c(2,3), as.numeric)


#next step is cdec turbitiy and temperature

X2 <- read.table(file="data/X2_daily.txt",header=T)
obs.temp.dat <- Tempx2_x 
obs.turb.dat <- turbx2_x
daylight = daylength(38.15, c(130:310))$Daylength*60
LT.fx <- daylight/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace

#day <- CDEC[,2]
#LT.fx <- daylength/max(daylength)
wtr.yr = c("BN", "W", "BN", "D", "C", "C", "BN", "W", "BN", "W", "D", "C", "C") #water year types 2010-2022
#wtr.yr<-c(1,1,1,1,1,1,4,4,2,3,3,1,4,5,4,3,1,3,4,5,5) # Sacto WY type wet = 1,above normal = 2, below normal = 3, dry = 4, critical = 5
#rep.yrs <- seq(21,26,by=1) # years 2015-2020 from 1995-2020 dataset, for prey density summary
#ex.strata <- c(13,10,8,6,3,1) # map IBMR strata to 5 CDEC regions


beta_hat <- read.table(file='data/beta_hat.txt') # MC filtering coefficients
#this is model uncertainty
#beta_hat <- read.table(file='BEM consumption rate project/beta_hat_IND.txt') # MC filtering coefficients


### 3. Model bioenergetics and reference points


all.VBL<-all.VBWt<- array(NA,dim=c((n.days), n.yrs,nrow(beta_hat)))
all.L<- data.frame(Day = NA, Tempave = NA,Temp_05 = NA,Temp_1 = NA,Temp_15 = NA,
                   Temp05 = NA,Temp1 = NA,Temp15 = NA, s = NA,level = NA)
all.Wt<- data.frame(Day = NA, Tempave = NA,Temp_05 = NA,Temp_1 = NA,Temp_15 = NA,
                    Temp05 = NA,Temp1 = NA,Temp15 = NA, s = NA,level = NA)
all.delta.BEM <- array(NA,dim=c((n.days),n.strata, n.yrs,nrow(beta_hat)))
all.delta.VB <- array(NA,dim=c((n.days), n.yrs,nrow(beta_hat)))
mn.pred.L<-mn.pred.Wt<- array(NA,dim=c((n.days+1), (n.yrs-1),nrow(beta_hat)))

#d is days, i is strata.

for (s in 1:nrow(beta_hat)) {
#for (s in 1:50) {
  VB.L <- VB.Wt <- delta.VB <- array(NA,dim=c(n.days,n.yrs))
  FOD.fx<-TMP.fx <- NTU.fx <- Cmax1 <- realized.Cmax1 <- C <- Limno <- energy <-delta.BEM <- array(NA,dim=c(n.days,n.strata,n.yrs)) # day x region x yr
  Wt <- L <- Foodsum <- array(NA,dim=c((n.days+1),n.strata,n.yrs)) # day x region x yr
  Food <- array(NA,dim=c(n.days,n.strata,n.yrs,n.prey)) # day x region x yr x prey type
  
  #logit.temp.b[1] <- runif(1,(0.43-0.25*0.43),(0.43+0.25*0.43))
  #logit.temp.c[1] <- runif(1,(-2.7-0.25*2.47),(-2.7+0.25*2.7))
  #logit.temp.c[2] <- runif(1,(24.44-0.25*24.44),(24.44+0.25*24.44))
  #R.Q <- rep(runif(1,(0.036-0.25*0.036),(0.036+0.25*0.036)),5)
  #min1 <- runif(1,(0.64-0.25*0.64),(0.64+0.25*0.64))
  #b.turb <- runif(1,(0.1-0.25*0.1),(0.1+0.25*0.1))
  
  a.c <- rep(beta_hat[s,1],5)
  b.c <- rep(beta_hat[s,2],5)
  T.M1 <- beta_hat[s,3]
  R.Q <- rep(beta_hat[s,4],5)
  min1 <- beta_hat[s,5]
  
  
  # start model
  #FIguring outthe June 1st starting length is prety iffy if you 
  for (h in 2:n.yrs) { # year
    L[1,,h] <- start.L[h] # June start
    VB.L[1,h] <- start.L[h] #  1st value is 2017-2020 mean
    VB.Wt[1,h] <- ln.a[2]*(VB.L[1,h]^ln.b[2])
    Wt[1,,h] <- ln.a[2]*(L[1,,h]^ln.b[2])
    delta.VB[1] <- NA
    
    for (i in 1:n.strata) { # region = Conf, NESuisun, LSac, Sac-Cache, SJR
      Food[1,i,h,] <- PD.mn.array[1,,ex.strata[i],h]*t(V[,4]/K[,4]) #this is selecting V and k for the fourth life stage
      Foodsum[1,i,h] <- sum(Food[1,i,h,c(3:5,9,11,12)])
      FOD.fx[1,i,h] <- sum(Food[1,i,h,]/(1+sum(Food[1,i,h,])))
      
      # Temp-consumption model
      L.1<- exp((1/(T.0[stage]-CQ[stage]))*log(0.98*(1-CK.1[stage])/(CK.1[stage]*0.02))*(as.numeric(obs.temp.dat[1,i,h])-CQ[stage]))
      L.2 <- exp((1/(T.L[stage]-T.M1))*log(0.98*(1-CK.4[stage])/(CK.4[stage]*0.02))*(T.L[stage]-as.numeric(obs.temp.dat[1,i,h])))
      K.A <- CK.1[stage]*L.1/(1+CK.1[stage]*(L.1-1))
      K.B <- CK.4[stage]*L.2/(1+CK.4[stage]*(L.2-1))
      TMP.fx[1,i,h] <- K.A*K.B
      
      NTU.fx[1,i,h] <- min1+(1-min1)/(1+exp(-(b.turb*(as.numeric(obs.turb.dat[1,i,h])-a.turb)))) # Turbidity-consumption
      Cmax1[1,i,h] <- a.c[4]*(Wt[1,i,h]^b.c[4])
      realized.Cmax1[1,i,h] <- Cmax1[1,i,h]*TMP.fx[1,i,h]*NTU.fx[1,i,h]*LT.fx[1]
      C.prey <- (realized.Cmax1[1,i,h]*Food[1,i,h,]) /(1+sum(Food[1,i,h,])) # C.prey : total consumption of each prey type. function of max consumption, food 
      energy.prey <- e.d*C.prey
      Limno[1,i,h] <- energy.prey[1]/sum(energy.prey) # fraction energy from Limno prey
      energy[1,i,h] <- e.d[1]*Limno[1,i,h] + e.d[2]*(1-Limno[1,i,h]) #sum(energy.prey)/n.prey
      C[1,i,h] <- sum(C.prey) # C = realized consumption rate
      R <- a.r[4]*(Wt[1,i,h]^b.r[4])*exp(R.Q[stage]*as.numeric(obs.temp.dat[1,i,h]))
      Feg <- F.a[4]*C[1,i,h]
      U <- U.a[4]*(C[1,i,h]-Feg)
      SDA <- S.d[4]*(C[1,i,h]-Feg)
      Wt.gain <- Wt[1,i,h]*(energy[1,i,h]/e.s)*(C[1,i,h]-R-Feg-U-SDA)
      Wt[2,i,h] <- Wt[1,i,h]+Wt.gain # add the growth
      L[2,i,h] <- (Wt[2,i,h]/ln.a[2])^(1/ln.b[2])
      delta.BEM[1,i,h] <- Wt[2,i,h]-Wt[1,i,h]
    }
    
    for (d in 2:n.days) {
      VB.L[d,h] <- VB.L[(d-1),h]+(Linf-VB.L[(d-1),h])*(1-exp(-VB.k*(1/365)))
      VB.Wt[d,h] <- ln.a[2]*(VB.L[d,h]^ln.b[2])
      delta.VB[d,h] <- VB.Wt[d,h]-VB.Wt[(d-1),h]
    }}
  
  #this is runnign the rest of the model
  for (i in 1:n.strata) { # temperature scenarios

    for (h in 2:n.yrs) { # zoop scenarios
      for (d in 2:n.days) { # day
        Food[d,i,h,] <- PD.mn.array[d,,ex.strata[i],h]*t(V[,4]/K[,4])
        Foodsum[d,i,h] <- sum(Food[d,i,h,c(3:5,9,11,12)])
        FOD.fx[d,i,h] <- sum(Food[d,i,h,]/(1+sum(Food[d,i,h,])))
        
        # Temp-consumption model
        L.1<- exp((1/(T.0[stage]-CQ[stage]))*log(0.98*(1-CK.1[stage])/(CK.1[stage]*0.02))*(as.numeric(obs.temp.dat[d,i,h])-CQ[stage]))
        L.2 <- exp((1/(T.L[stage]-T.M1))*log(0.98*(1-CK.4[stage])/(CK.4[stage]*0.02))*(T.L[stage]-as.numeric(obs.temp.dat[d,i,h])))
        K.A <- CK.1[stage]*L.1/(1+CK.1[stage]*(L.1-1))
        K.B <- CK.4[stage]*L.2/(1+CK.4[stage]*(L.2-1))
        TMP.fx[d,i,h] <- K.A*K.B
        
        NTU.fx[d,i,h] <- min1+(1-min1)/(1+exp(-(b.turb*(as.numeric(obs.turb.dat[d,i,h])-a.turb)))) # Turbidity-consumption
        Cmax1[d,i,h] <- a.c[4]*(Wt[d,i,h]^b.c[4])
        realized.Cmax1[d,i,h] <- Cmax1[d,i,h]*TMP.fx[d,i,h]*NTU.fx[d,i,h]*LT.fx[d]
        C.prey <- (realized.Cmax1[d,i,h]*Food[d,i,h,]) /(1+sum(Food[d,i,h,])) # C.prey : total consumption of each prey type. function of max consumption, food 
        energy.prey <- e.d*C.prey
        Limno[d,i,h] <- energy.prey[1]/sum(energy.prey) # fraction energy from Limno prey
        energy[d,i,h] <- e.d[1]*Limno[d,i,h] + e.d[2]*(1-Limno[d,i,h]) #sum(energy.prey)/n.prey
        C[d,i,h] <- sum(C.prey) # C = realized consumption rate
        R <- a.r[4]*(Wt[d,i,h]^b.r[4])*exp(R.Q[stage]*as.numeric(obs.temp.dat[d,i,h]))
        Feg <- F.a[4]*C[d,i,h]
        U <- U.a[4]*(C[d,i,h]-Feg)
        SDA <- S.d[4]*(C[d,i,h]-Feg)
        Wt.gain <- Wt[d,i,h]*(energy[d,i,h]/e.s)*(C[d,i,h]-R-Feg-U-SDA)
        Wt[(d+1),i,h] <- Wt[d,i,h]+Wt.gain # add the growth
        L[(d+1),i,h] <- (Wt[(d+1),i,h]/ln.a[2])^(1/ln.b[2])
        delta.BEM[d,i,h] <- Wt[(d+1),i,h]-Wt[d,i,h]
        if(is.na(Wt[(d+1),i,h])) {
          print(paste("NA's in weight~!", "d=", d, "temp scenario", i, "zoop scenario",h))
        }
      }}}
  
  ### 4. Get summaries
  L[153,,] #terminal lengths
  
  # for (d in 1:(n.days+1)) {
  #   mn.pred.L[d,s,] <- c(mean(L[d,3:6,2],na.rm=T),
  #                        mean(L[d,2:6,3],na.rm=T),
  #                        mean(L[d,3:6,4],na.rm=T),
  #                        mean(L[d,1:6,5],na.rm=T),
  #                        mean(L[d,3:6,6],na.rm=T))
  #   mn.pred.Wt[d,s,] <- c(mean(Wt[d,3:6,2],na.rm=T),
  #                         mean(Wt[d,2:6,3],na.rm=T),
  #                         mean(Wt[d,3:6,4],na.rm=T),
  #                         mean(Wt[d,1:6,5],na.rm=T),
  #                         mean(Wt[d,3:6,6],na.rm=T))
  # }
  
  
  all.delta.BEM[,,,s] <- delta.BEM
  all.delta.VB[,,s] <- delta.VB
  all.L <- bind_rows(all.L, data.frame(Day = c(152:305), 
                                       `Temp_05` = L[,1,], 
                                       `Temp_1` = L[,2,],
                                       `Temp_15` = L[,3,],
                                       `Temp05` = L[,4,],
                                       `Temp1` = L[,5,],
                                       `Temp15` = L[,6,],
                                       `Tempave` = L[,7,],
                                       level = h, s = s))
  all.VBL[,,s] <- VB.L
  all.Wt <- bind_rows(all.Wt, data.frame(Day = c(152:305), 
                                         `Temp_05` = Wt[,1,], 
                                         `Temp_1` = Wt[,2,],
                                         `Temp_15` = Wt[,3,],
                                         `Temp05` = Wt[,4,],
                                         `Temp1` = Wt[,5,],
                                         `Temp15` = Wt[,6,],
                                         `Tempave` = Wt[,7,],
                                         level = h, s = s))
  all.VBWt[,,s] <- VB.Wt
  
}


#well, that didn't quite work.
all.Wt2_ps = all.Wt %>%
  pivot_longer(cols = c(Temp_05.1:Tempave.11), names_to = "Temp", values_to = "Weight") %>%
  filter(!is.na(Weight)) %>%
  select(s,Temp, Day, Weight) %>%
  mutate(Temp2 = str_sub(Temp, start = 1, end = -3),
         Temp2 = str_remove_all(Temp2, "\\.+"),
         zoopX = str_sub(Temp, start = -2, end = -1))

# ggplot(all.Wt2, aes(x = Day, y = Weight, color = Region2)) +
#   geom_smooth()+
#   facet_wrap(~Year)

#Uh, this is taking too long
pseudosave$pdiapfor = round(pseudosave$pdiapfor, digits =2)
Wt.summary_ps = all.Wt2_ps %>%
  group_by(Temp2, zoopX, Day) %>%
  summarize(MWeight = mean(Weight, na.rm =T), sdWeight = sd(Weight, na.rm =T)) %>%
  mutate(zoop = case_match(zoopX,"13" ~ pseudosave$pdiapfor[12], ".2" ~ pseudosave$pdiapfor[1], ".3" ~pseudosave$pdiapfor[2],
                           ".4" ~ pseudosave$pdiapfor[3], ".5" ~ pseudosave$pdiapfor[4], ".6" ~ pseudosave$pdiapfor[5],
                           ".7" ~ pseudosave$pdiapfor[6], ".8" ~ pseudosave$pdiapfor[7], ".9" ~ pseudosave$pdiapfor[8], 
                           "10" ~ pseudosave$pdiapfor[9], "11" ~ pseudosave$pdiapfor[10], "12" ~pseudosave$pdiapfor[11]))

#some how my temperature labels got scrambled. Sigh
Wt.summary_ps = mutate(Wt.summary_ps, Temp = factor(Temp2, levels = c("Temp_15", "Temp_1", "Temp_05", "Tempave",
                                                                       "Temp05", "Temp1", "Temp15"),
                                                    labels = c(21.86, 21.36, 20.86, 20.36, 19.8646, 19.3646, 18.86)),
                       Temp = as.numeric(as.character(Temp)))



ggplot(Wt.summary_ps, aes(x = Day, y = MWeight, color = as.factor(Temp))) +
  geom_smooth()+
  facet_wrap(~zoop)

ggplot(Wt.summary_ps, aes(x = Day, y = MWeight, color = as.factor(zoop))) +
  geom_smooth()+
  facet_wrap(~Temp)



#OK, how much did they grow over the whole summer (or fall)?
Wtsum2_ps = group_by(Wt.summary_ps, zoop, Temp) %>%
  summarize(startweight = first(MWeight), endwieght = last(MWeight), diffweight = endwieght-startweight)


ggplot(Wtsum2_ps, aes(x = zoop, y = diffweight, fill = as.factor(Temp))) + geom_col(position = "dodge")+
  facet_wrap(~zoop, scales = "free_x")+
  theme(axis.text.x = element_blank())

ggplot(Wtsum2_ps, aes(x = zoop, y = diffweight, color = as.factor(Temp))) + geom_line()+
  ylab("Total summer/fall growth (g)")+ xlab("average pseudodiaptomus biomass")+
  scale_color_brewer(palette = "Dark2", name = "Temperature\nscenario")+ theme_bw()

ggplot(Wtsum2_ps, aes(x = log(zoop), y = diffweight, color = as.factor(Temp))) + geom_line()+ geom_point()

ggplot(Wtsum2_ps, aes(x = Temp, y = diffweight, color = as.factor(zoop))) + geom_line()+
  ylab("total summer-fall growth (g)") + xlab("Temperature regime")+
  scale_color_brewer(palette = "Set3", name = "Zooplankton \nScenario")+
  theme_bw()

#######################################
#change in weight per change in zoops. 

test = filter(Wtsum2_ps, Temp == first(unique(Wt.summary_ps$Temp))) %>%
  ungroup() %>%
  arrange(zoop) %>%
  mutate(lagzoop = lag(zoop), diffzoop = zoop-lagzoop, lagweight = lag(endwieght), diffweigth = endwieght-lagweight,
         ratio = diffweight/(diffzoop*1000))


ggplot(test, aes(x = zoop, y = ratio)) + geom_line()+
  ylab("Increase in weight per increase in zooplankton biomass")+
  xlab("total zooplankton biomass")

#so, I think this is telling me that increasing pseudodiaptomus from 2 to 20 makes a big dofferemce.
#but increasing it from 20 to 200 doesn't do that much for you. There is a threshold effect
#around about 30 or 40. 
save(Wtsum2_ps, Wt.summary_ps, all.Wt2_ps, file = "outputs/increasedzoops.Rdata")
load("outputs/increasedzoops.Rdata")


#So, how mcuh zooplankton do we actually see? And how mcuh does flow actions increase it?
#this is juveniles and adults together, so not quite right, but still helpful
load("C:/Users/rhartman/OneDrive - California Department of Water Resources/salinity control gates/SFHA_synthesis/data/pseudo_data_mass.RData")

#need to divide mass by 1000 to convert to mg.

pseudo_recentx = mutate(pseudo_data_mass2, BPUEmg = BPUE/1000, doy = yday(Date))

ggplot(filter(pseudo_recentx, doy %in% c(150:243), Year >2009), 
       aes(x = as.factor(Year), y = BPUEmg)) + geom_boxplot()+
  geom_hline(yintercept = 30)

#OK, so we are rarely over the point of satiation. I think.

ggplot(pseudo_recent, aes(x = Month, fill = as.factor(Year), y = exp(BPUE)/1000))+ geom_col(position = "dodge")+
  facet_wrap(~Region)+
  coord_cartesian(xlim = c(4,8))
#Egh, i don't think i back-transformed this properly

#OK, things to do:
## figure out a better way of growth rate at less than 6
## Put variability back in there
## What do i do with this zoop threshold analysis?
