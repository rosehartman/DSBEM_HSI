#predicted growth from the smelt cages - with zooplankton held constant across time
#rosemary hartman 7/1/2024

# 1. Read model parameters and data
# 2. Make a plot of Cmax, R, Cmax-temp.fx, and Cmax-turb.fx - probably don't need to do this
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
n.strata <- 2 #just Rio Vista and Belden's landing
n.prey <- 13 #may adjust this later
n.days<-46 # cages were out for six weeks
ln.a<-c(0.000005,0.00000183) # length-weight parameters
ln.b<-c(3,3.38)
startY <- 2023
endY <- 2023
yr.seq<- 1 # just one year
n.yrs<-length(yr.seq)
stage =4

#I'm not going to change any of these
VB.k<-2.72 #2.98 # VBGM parameters estimated in DSM TN 35
Linf<-78.39 #76.1
t0<--0.026 #-0.014

#need to get these from the pre-deployment lengths. He origionally had four weights, is that importnat?
start.L = read_csv("SmeltLengths.csv") 

#I had more end weights than start weights, which might be important
terminal.L <- read_csv("SmeltLengths_end.csv")
terminal.L %>%
  group_by(Location) %>%
  summarize(FL = mean(ForkLength))
#terminal.L <- c((start.L[1]+30.025),(start.L[2]+30.025),49.6,53.5,52.1,57) # EDSM observed, (10/15-11/15),2017,2018,2019

#I can probably use real weights rather than length-weight equations
start.Wt <- read_csv("SmeltLengths.csv")$Weight 
terminal.Wt <- read_csv("SmeltLengths_end.csv")$Weight 

#i can also use real standard errors, I guess
terminal.Wtse <- c(0.299,0.356,0.367)
mn.terminal.Wtse <- mean(terminal.Wtse)
terminal.Wtse <- c(mn.terminal.Wtse ,mn.terminal.Wtse ,terminal.Wtse,mn.terminal.Wtse)

#i'm not sure where these goes
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

#these could be updated
e.d <- c(2590,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590)		# energy density of prey items

#I kinda want to make thiese into data frames so they make sense
V <- matrix(NA,n.prey,5)					# V.ij is vulnerability of prey type j to fish i. 
#Set to 1 for all life stages eating all zooplankton types; except DS larvae values were 0 except for Limnoithona

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
V[13,] <- c(0,0,1,1,1) #amphipod

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
K[13,] <- c(NA,250,120,4.5,1) #amphipod (coppied from other cladocera)
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

#source("BEM consumption functions_v3.2.R") # load consumption functions

#OK, here is where I replace this table with the one from Belden's and Rio Vista
source("data organization.r")

# I need to get daily X2's from last summer - though it didn't look like thiswas used?
X2 <- read.table(file="X2_daily.txt",header=T)
obs.temp.dat <- select(ungroup(CDEC_wide), Temp.BDL, Temp.RVB)
obs.turb.dat <- select(ungroup(CDEC_wide), NTU.BDL, NTU.RVB)
daylength<- CDEC_wide$daylight
day <- CDEC_wide$Seq.Day
LT.fx <- daylength/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace
wtr.yr<-1  # Sacto WY type wet = 1,above normal = 2, below normal = 3, dry = 4, critical = 5. 2023 was wet


beta_hat <- read.table(file='beta_hat.txt') # MC filtering coefficients
#this is model uncertainty
#beta_hat <- read.table(file='BEM consumption rate project/beta_hat_IND.txt') # MC filtering coefficients


### 3. Model bioenergetics and reference points

#OK, this just fills in prey densities. I can do this with the zoop data. 
PD.mn.array = zoopsx2 #days by prey by strata


all.VBL<-all.VBWt<- array(NA,dim=c((n.days),nrow(beta_hat)))
all.L<- data.frame(Day = NA, Montezuma = NA, RioVista = NA, s = NA)
all.Wt<- data.frame(Day = NA, Montezuma = NA, RioVista = NA, s = NA)
all.delta.BEM <- array(NA,dim=c((n.days),n.strata,nrow(beta_hat)))
all.delta.VB <- array(NA,dim=c((n.days),nrow(beta_hat)))
mn.pred.L<-mn.pred.Wt<- array(NA,dim=c((n.days+1),nrow(beta_hat)))

#d is days, i is strata.

for (s in 1:nrow(beta_hat)) {
  VB.L <- VB.Wt <- delta.VB <- array(NA,dim=c(n.days))
  FOD.fx<-TMP.fx <- NTU.fx <- Cmax1 <- realized.Cmax1 <- C <- Limno <- energy <-delta.BEM <- array(NA,dim=c(n.days,n.strata)) # day x region
  Wt <- L <- Foodsum <- array(NA,dim=c((n.days+1),n.strata)) # day x region 
  Food <- array(NA,dim=c(n.days, n.prey, n.strata)) # day x  prey typexregion  x
  
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
  
  # start model with length/weights as fish went into the cages
  L[1,] <- mean(start.L) # all fish are mean of starting fish lengths
  VB.L[1] <- mean(start.L)
  VB.Wt[1] <- mean(start.Wt)
  Wt[1,] <- mean(start.Wt)
  delta.VB[1] <- NA
  
  for (i in 1:n.strata) { 
    Food[1,,i] <- PD.mn.array[1,,i]*t(V[,4]/K[,4])
    Foodsum[1,i] <- sum(Food[1,c(3:5,9,11,12),i])
    FOD.fx[1,i] <- sum(Food[1,,i]/(1+sum(Food[1,,i])))
    
    # Temp-consumption model
    L.1<- exp((1/(T.0[stage]-CQ[stage]))*log(0.98*(1-CK.1[stage])/(CK.1[stage]*0.02))*(as.numeric(obs.temp.dat[1,i])-CQ[stage]))
    L.2 <- exp((1/(T.L[stage]-T.M1))*log(0.98*(1-CK.4[stage])/(CK.4[stage]*0.02))*(T.L[stage]-as.numeric(obs.temp.dat[1,i])))
    K.A <- CK.1[stage]*L.1/(1+CK.1[stage]*(L.1-1))
    K.B <- CK.4[stage]*L.2/(1+CK.4[stage]*(L.2-1))
    TMP.fx[1,i] <- K.A*K.B
    
    NTU.fx[1,i] <- min1+(1-min1)/(1+exp(-(b.turb*(as.numeric(obs.turb.dat[1,i])-a.turb)))) # Turbidity-consumption
    Cmax1[1,i] <- a.c[4]*(Wt[1,i]^b.c[4])
    realized.Cmax1[1,i] <- Cmax1[1,i]*TMP.fx[1,i]*NTU.fx[1,i]*LT.fx[1]
    C.prey <- (realized.Cmax1[1,i]*Food[1,,i]) /(1+sum(Food[1,,i])) # C.prey : total consumption of each prey type. function of max consumption, food 
    energy.prey <- e.d*C.prey
    Limno[1,i] <- energy.prey[1]/sum(energy.prey) # fraction energy from Limno prey
    energy[1,i] <- e.d[1]*Limno[1,i] + e.d[2]*(1-Limno[1,i]) #sum(energy.prey)/n.prey
    C[1,i] <- sum(C.prey) # C = realized consumption rate
    R <- a.r[4]*(Wt[1,i]^b.r[4])*exp(R.Q[stage]*as.numeric(obs.temp.dat[1,i]))
    Feg <- F.a[4]*C[1,i]
    U <- U.a[4]*(C[1,i]-Feg)
    SDA <- S.d[4]*(C[1,i]-Feg)
    Wt.gain <- Wt[1,i]*(energy[1,i]/e.s)*(C[1,i]-R-Feg-U-SDA)
    Wt[2,i] <- Wt[1,i]+Wt.gain # add the growth
    L[2,i] <- (Wt[2,i]/ln.a[2])^(1/ln.b[2])
    delta.BEM[1,i] <- Wt[2,i]-Wt[1,i]
  }
  
  for (d in 2:n.days) {
    VB.L[d] <- VB.L[(d-1)]+(Linf-VB.L[(d-1)])*(1-exp(-VB.k*(1/365)))
    VB.Wt[d] <- ln.a[2]*(VB.L[d]^ln.b[2])
    delta.VB[d] <- VB.Wt[d]-VB.Wt[(d-1)]
  }
  
  #now this is the rest of the days
  
  for (i in 1:n.strata) { # region = Montezuma, Rio Vista
    for (d in 2:n.days) { # day
      Food[d,,i] <- PD.mn.array[d,,i]*t(V[,4]/K[,4])
      Foodsum[d,i] <- sum(Food[d,c(3:5,9,11,12),i])
      FOD.fx[d,i] <- sum(Food[d,,i]/(1+sum(Food[d,,i])))
      
      # Temp-consumption model
      L.1<- exp((1/(T.0[stage]-CQ[stage]))*log(0.98*(1-CK.1[stage])/(CK.1[stage]*0.02))*(as.numeric(obs.temp.dat[d,i])-CQ[stage]))
      L.2 <- exp((1/(T.L[stage]-T.M1))*log(0.98*(1-CK.4[stage])/(CK.4[stage]*0.02))*(T.L[stage]-as.numeric(obs.temp.dat[d,i])))
      K.A <- CK.1[stage]*L.1/(1+CK.1[stage]*(L.1-1))
      K.B <- CK.4[stage]*L.2/(1+CK.4[stage]*(L.2-1))
      TMP.fx[d,i] <- K.A*K.B
      
      NTU.fx[d,i] <- min1+(1-min1)/(1+exp(-(b.turb*(as.numeric(obs.turb.dat[d,i])-a.turb)))) # Turbidity-consumption
      Cmax1[d,i] <- a.c[4]*(Wt[d,i]^b.c[4])
      realized.Cmax1[d,i] <- Cmax1[d,i]*TMP.fx[d,i]*NTU.fx[d,i]*LT.fx[d]
      C.prey <- (realized.Cmax1[d,i]*Food[d,,i]) /(1+sum(Food[d,,i])) # C.prey : total consumption of each prey type. function of max consumption, food 
      energy.prey <- e.d*C.prey
      Limno[d,i] <- energy.prey[1]/sum(energy.prey) # fraction energy from Limno prey
      energy[d,i] <- e.d[1]*Limno[d,i] + e.d[2]*(1-Limno[d,i]) #sum(energy.prey)/n.prey
      C[d,i] <- sum(C.prey) # C = realized consumption rate
      R <- a.r[4]*(Wt[d,i]^b.r[4])*exp(R.Q[stage]*as.numeric(obs.temp.dat[d,i]))
      Feg <- F.a[4]*C[d,i]
      U <- U.a[4]*(C[d,i]-Feg)
      SDA <- S.d[4]*(C[d,i]-Feg)
      Wt.gain <- Wt[d,i]*(energy[d,i]/e.s)*(C[d,i]-R-Feg-U-SDA)
      Wt[(d+1),i] <- Wt[d,i]+Wt.gain # add the growth
      L[(d+1),i] <- (Wt[(d+1),i]/ln.a[2])^(1/ln.b[2])
      delta.BEM[d,i] <- Wt[(d+1),i]-Wt[d,i]
    }}
  
  ### 4. Get summaries
  L[47,] #terminal lengths
  
  all.delta.BEM[,,s] <- delta.BEM
  all.delta.VB[,s] <- delta.VB
  all.L <- bind_rows(all.L, data.frame(Day = c(1:47), Montezuma = L[,1], RioVista = L[,2], s = s))
  all.VBL[,s] <- VB.L
  all.Wt <- bind_rows(all.Wt, data.frame(Day = c(1:47), Montezuma = Wt[,1], RioVista = Wt[,2], s = s))
  all.VBWt[,s] <- VB.Wt
  
}

# ratio of predicted to reference growth
ratio<-array(NA,dim=c(n.days,n.strata))
for (i in 1:n.strata) {
  for (t in 1:n.days-1) {
    mn <- mean(all.delta.BEM[t,i,],na.rm=T)
    ratio[t,i] <- mn/delta.VB[t+1]
  }}


all.Llong = pivot_longer(all.L, cols = c(RioVista, Montezuma), names_to = "Location", values_to = "ForkLength") %>%
  filter(!is.na(ForkLength))

all.Wtlong =  pivot_longer(all.Wt, cols = c(RioVista, Montezuma), names_to = "Location", values_to = "Weight") %>%
  filter(!is.na(Weight))

# #put the actual lengths on there
# actuallenghts = data.frame(Day = c(rep(1, length(start.L)), rep(47, nrow(terminal.L)) ), ForkLength = c(start.L, terminal.L$ForkLength),
#                            Location = c(rep("Start", length(start.L)), terminal.L$Location)) %>%
#   filter(Location != "FCCL") %>%
#   mutate(Location = case_when(Location == "Rio Vista" ~ "RioVista",
#                               TRUE ~ Location))

ggplot(all.Llong) + #geom_point( alpha = 0.05) + 
 # geom_boxplot(data = actuallenghts, aes(x = Day, y = ForkLength, fill = Location), inherit.aes = FALSE)+
  geom_smooth(aes(x = Day, y = ForkLength, color = Location))+
  scale_color_manual(values = c("blue", "grey30", "red"))+
  theme_bw()+
  scale_fill_manual(values = c("blue", "grey30", "red"))

write.csv(all.Llong, "outputs/smeltcagelengths_same zoops.csv")

ggplot(all.Wtlong) + #geom_point( alpha = 0.05) + 
  #geom_boxplot(data = actuallenghts, aes(x = Day, y = Weight, fill = Location), inherit.aes = FALSE)+
  geom_smooth(aes(x = Day, y = Weight, color = Location))+
  scale_color_manual(values = c("blue", "grey30", "red"))+
  theme_bw()+
  scale_fill_manual(values = c("blue", "grey30", "red"))
