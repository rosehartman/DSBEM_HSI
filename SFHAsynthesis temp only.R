
#Version with constant zoops and turbidity

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

### 1. Read model parameters and data
n.strata <- 7 #"Confluence"   "NE Suisun"    "NW Suisun"    "SE Suisun"    "Suisun Marsh" "SW Suisun" "Lower Sacramento RIver"
ex.strata = sort(c( "NE Suisun", "SE Suisun", "NW Suisun", "Lower Sacramento River", "Confluence", "SW Suisun",  "Suisun Marsh"))
n.prey <- 12 #may adjust this later
n.days<-153 # summer and fall
ln.a<-c(0.000005,0.00000183) # length-weight parameters
ln.b<-c(3,3.38)
startY <- 2010
endY <- 2022 #expand through 2024 when all data are in
yr.seq<- seq(startY-2010,endY-2010,by=1)
n.yrs<-13
stage =4

#I'm not going to change any of these
VB.k<-2.72 #2.98 # VBGM parameters estimated in DSM TN 35
Linf<-78.39 #76.1
t0<--0.026 #-0.014

#This is what Will origionally used, may want something else, based on observed values in each year. 
#OK, I need the june first lengths from skt and edsm
#but if i'm looking at growth potential, maybe just start them all off the same is better?
#start.L = c(24.3,29.2,22.2,19.9,21.9,28.1) # assumed length on Jun 1, 1st value = EDSM+SKT observed, (5/15-6/15),2017,2018,2019
start.L = c(rep(23, 13))
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


source("BEM consumption functions_v3.2.R") # load consumption functions

#next step is cdec turbitiy and temperature
load("data/WaterQuality20102022.RData")
X2 <- read.table(file="data/X2_daily.txt",header=T)
obs.temp.dat <- Tempx2
obs.turb.dat <- turbx2_constant
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


####Sensativity analysis ##########################
#First hold zooplankton constant and run everything else as-is.
PD.mn.array = zoopx2_constant #days by prey by strata



all.VBL<-all.VBWt<- array(NA,dim=c((n.days), n.yrs,nrow(beta_hat)))
all.L<- data.frame(Day = NA, Confluence = NA, `Lower Sacramento River` = NA, `NE Suisun` = NA,
                   `NW Suisun` = NA,`SE Suisun` = NA,
                   `Suisun Marsh` = NA, `SW Suisun` = NA, s = NA,Yr = NA)
all.Wt<- data.frame(Day = NA, Confluence = NA, `Lower Sacramento River` = NA, `NE Susiun` = NA,
                    `NW Susiun` = NA,`SE Susiun` = NA,
                    `Suisun Marsh` = NA, `SW Suisun` = NA, s = NA, Yr = NA)
all.delta.BEM <- array(NA,dim=c((n.days),n.strata, n.yrs,nrow(beta_hat)))
all.delta.VB <- array(NA,dim=c((n.days), n.yrs,nrow(beta_hat)))
mn.pred.L<-mn.pred.Wt<- array(NA,dim=c((n.days+1), (n.yrs-1),nrow(beta_hat)))

#d is days, i is strata.

for (s in 1:nrow(beta_hat)) {
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
  for (i in 1:n.strata) { # region = Conf, NESuisun, LSac
    for (h in 2:n.yrs) { # year
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
                                       `NE Suisun` = L[,3,], 
                                       `SE Suisun` = L[,5,],
                                       `NW Suisun` = L[,4,],
                                       `Lower Sacramento River` = L[,2,],
                                       `Confluence` = L[,1,],
                                       `SW Suisun` = L[,7,],
                                       `Suisun Marsh` = L[,6,],
                                       Year = h, s = s))
  all.VBL[,,s] <- VB.L
  all.Wt <- bind_rows(all.Wt, data.frame(Day = c(152:305), 
                                         `NE Suisun` = Wt[,3,], 
                                         `SE Suisun` = Wt[,5,],
                                         `NW Suisun` = Wt[,4,],
                                         `Lower Sacramento River` = Wt[,2,],
                                         `Confluence` = Wt[,1,],
                                         `SW Suisun` = Wt[,7,],
                                         `Suisun Marsh` = Wt[,6,],
                                         Yr = h, s = s))
  all.VBWt[,,s] <- VB.Wt
  
}


#well, that didn't quite work.
all.Wt2_constantzoops = all.Wt %>%
  pivot_longer(cols = c(NE.Suisun.1:Suisun.Marsh.13), names_to = "Region", values_to = "Weight") %>%
  filter(!is.na(Weight)) %>%
  select(s, Region, Yr, Day, Weight) %>%
  mutate(Region2 = str_sub(Region, start = 1, end = -3),
         Region2 = str_remove_all(Region2, "\\.+"),
         YearX = str_sub(Region, start = -2, end = -1))

# ggplot(all.Wt2, aes(x = Day, y = Weight, color = Region2)) +
#   geom_smooth()+
#   facet_wrap(~Year)

#Uh, this is taking too long
Wt.summary_constantzoops = all.Wt2_constantzoops %>%
  group_by(Region2, YearX, Day) %>%
  summarize(MWeight = mean(Weight, na.rm =T), sdWeight = sd(Weight, na.rm =T)) %>%
  mutate(Year = case_match(YearX,"13" ~ 2022, ".2" ~ 2011, ".3" ~2012, ".4" ~ 2013, ".5" ~ 2014, ".6" ~ 2015,
                           ".7" ~ 2016, ".8" ~ 2017, ".9" ~ 2018, "10" ~ 2019, "11" ~ 2020, "12" ~2021))




Wt.summary_constantzoops = left_join(Wt.summary_constantzoops, sal) %>%
  mutate(MWeight6 = case_when(GoodSalinity ~ MWeight,
                              TRUE ~ 0),
         MWeight2 = case_when(GoodSalinity2 ~ MWeight,
                              TRUE ~ 0))

ggplot(Wt.summary_constantzoops, aes(x = Day, y = MWeight, color = Region2)) +
  geom_smooth()+
  facet_wrap(~Year)

ggplot(Wt.summary_constantzoops, aes(x = Day, y = MWeight, color = as.factor(Year))) +
  geom_smooth()+
  facet_wrap(~Region)


Wt.sumx_constantzoops = left_join(Wt.summary_constantzoops, DF)

ggplot(Wt.sumx_constantzoops, aes(x = Day, y = MWeight, color = X2, group = as.factor(Year))) +
  geom_smooth()+
  facet_wrap(~Region2)+
  ylab("Predicted Weight in Grams")


#OK, how much did they grow over the whole summer (or fall) and is it related to X2?
Wtsum2_constantzoops = group_by(Wt.summary_constantzoops, Year, Region2) %>%
  summarize(startweight = first(MWeight), endwieght = last(MWeight), diffweight = endwieght-startweight)


ggplot(Wtsum2_constantzoops, aes(x = Region2, y = diffweight, fill = Region2)) + geom_col(position = "dodge")+
  facet_wrap(~Year, scales = "free_x")+
  theme(axis.text.x = element_blank())



Wtsum2_constantzoops = left_join(Wtsum2_constantzoops, DF)

ggplot(Wtsum2_constantzoops, aes(x = X2, y = diffweight)) + geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Region2, scales = "free_x")+
  
  ylab("Predicted growth (g), June-OCtober")+
  xlab("Mean June-September X2")
#huh,  not seeing any effect of X2, maybe because most of the X2 effect was due to zoops?


Wt.summary_constantzoops = Wt.summary_constantzoops %>%
  group_by(Region2, Year) %>%
  mutate(WperD = MWeight - lag(MWeight),
         WperD = replace_na(WperD, 0)) %>%
  ungroup()

Meangrowth6_constantzoops = Wt.summary_constantzoops %>%
  group_by(Year, Day) %>%
  filter(GoodSalinity) %>%
  summarise(Growth = mean(WperD, na.rm =T)) %>%
  group_by(Year) %>%
  mutate(Weight = cumsum(Growth)+0.073298) %>%
  ungroup()

Meangrowth2_constantzoops = Wt.summary_constantzoops %>%
  group_by(Year, Day) %>%
  filter(GoodSalinity2) %>%
  summarise(Growth2 = mean(WperD, na.rm =T)) %>%
  group_by(Year) %>%
  mutate(Weight2 = cumsum(Growth2)+0.073298) %>%
  ungroup()

Meangrowth_constantzoops = left_join(Meangrowth6_constantzoops, Meangrowth2_constantzoops)

ggplot(Meangrowth_constantzoops, aes(x = Day, y = Weight,color = as.factor(Year))) + geom_line()
ggplot(Meangrowth_constantzoops, aes(x = Day, y = Weight2,color = as.factor(Year))) + geom_line()

#Does this make sense? need to check with the group. 

Meangrowthsum_constantzoops = Meangrowth_constantzoops %>%
  group_by(Year) %>%
  summarise(Weight = max(Weight, na.rm =T), Weight2 = max(Weight2, na.rm =T),) %>%
  left_join( DF) %>%
  left_join(Meantemps)

#x2 versus growth 6 psu
ggplot(Meangrowthsum_constantzoops, aes(x = X2, y = Weight)) + geom_point()+ geom_smooth(method = "lm")+
  geom_text(aes(label = Year))+ ggtitle("Only temperature changes")

ggplot(Meangrowthsum_constantzoops, aes(x = X2, y = Weight)) + geom_point()+ geom_smooth()+
  geom_text(aes(label = Year))



ggplot(Meangrowthsum_constantzoops, aes(x = Temperature, y = Weight)) + geom_point()+ geom_smooth(method = "lm")+
  geom_text(aes(label = Year))


#x2 versus growth 2 psu
ggplot(Meangrowthsum_constantzoops, aes(x = X2, y = Weight2)) + geom_point()+ geom_smooth(method = "lm")+
  geom_text(aes(label = Year))


ggplot(Meangrowthsum_constantzoops, aes(x = Temperature, y = Weight2)) + geom_point()+ geom_smooth(method = "lm")+
  geom_text(aes(label = Year))

