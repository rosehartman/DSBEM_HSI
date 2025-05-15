#predicted growth from the smelt cages
#rosemary hartman 6/19/2024

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
n.prey <- 12 #may adjust this later
#n.days<-46 # cages were out for six weeks
ln.a<-c(0.000005,0.00000183) # length-weight parameters
ln.b<-c(3,3.38)
#startY <- 2023
#endY <- 2023
yr.seq<- 1 # just one year
n.yrs<-length(yr.seq)
stage =4

#I'm not going to change any of these
VB.k<-2.72 #2.98 # VBGM parameters estimated in DSM TN 35
Linf<-78.39 #76.1
t0<--0.026 #-0.014

#need to get these from the pre-deployment lengths. He origionally had four weights, is that importnat?
start.L = read_csv("data/SmeltLengths.csv") 

#I had more end weights than start weights, which might be important
terminal.L <- read_csv("data/SmeltLengths_end.csv")
terminal.L %>%
  group_by(Location) %>%
  summarize(FL = mean(ForkLength))
#terminal.L <- c((start.L[1]+30.025),(start.L[2]+30.025),49.6,53.5,52.1,57) # EDSM observed, (10/15-11/15),2017,2018,2019

#I can probably use real weights rather than length-weight equations
start.Wt <- read_csv("data/SmeltLengths.csv")$Weight 
terminal.Wt <- read_csv("data/SmeltLengths_end.csv")$Weight 

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
T.M <- c(20,20,23,21.6,23) # Temperature at 0.98 of maximum (deg C) Updated to 21.6 based on Smith and Nobriga 2023
# Metabolism (R) parameters
a.r <- c(0.0027,0.0027,0.0027,0.0027,0.0027)	# weight multiplier
b.r <- c(-0.216,-0.216,-0.216,-0.216,-0.216)	# weight exponent
R.Q <- c(0.036,0.036,0.036,0.064,0.036)			# exponent for temperature effect Updated to 0.64 based on Smith and Nobriga 2023
S.d <- c(0.175,0.175,0.175,0.175,0.175)			# Fraction of assimilated food lost to SDA (specific dynamic action)

# Egestion (F) and excretion (U) parameters
F.a <- c(0.16,0.16,0.16,0.16,0.16)				# Fraction of consumed food lost to egestion
U.a <- c(0.1,0.1,0.1,0.1,0.1)				# Fraction of assimilated food lost to excretion

e.s <- 4814									# J/g: convert g(prey)/g(delta smelt) to g(smelt)/g(smelt) - fixed

#these could be updated
e.d <- c(1823,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590)		# energy density of prey items

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
V[10,] <- c(0,0,1,1,1) #other - This was set at zero, but maybe it should be 1?
V[11,] <- c(0,0,1,1,1) #eurytem ***
V[12,] <- c(0,0,1,1,1) #pdiapfor ***
#V[13,] <- c(0,0,1,1,1) #amphipod

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
#K[13,] <- c(NA,250,120,4.5,1) #amphipod (coppied from other cladocera)


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
#source("data organization.r")
load("data/cagezoops_byyear.Rdata")
load("data/CDEC_wide_cages.RData")
beta_hat <- read.table(file='data/beta_hat.txt') # MC filtering coefficients
#this is model uncertainty
#beta_hat <- read.table(file='BEM consumption rate project/beta_hat_IND.txt') # MC filtering coefficients


### 3. Model bioenergetics and reference points



#d is days, i is strata.

cagegrowth = function(PD.mn.array, CDEC_wide, start.Wt, beta_hat) {
  obs.temp.dat <- select(ungroup(CDEC_wide), Temp.BDL, Temp.RVB)
  obs.turb.dat <- select(ungroup(CDEC_wide), NTU.BDL, NTU.RVB)
  Locations = c("BDL", "RVB")
  CDEC_wide$Seq.Day = c(1:nrow(CDEC_wide))
  n.days = nrow(CDEC_wide)
  
  #percent maximum daylight (based on daylight during the solstacve)
  CDEC_wide$Frc.Max.Daylight = CDEC_wide$daylight/(daylength(38.15, 173)$Daylength*60)
  
  
  daylength<- CDEC_wide$daylight
  day <- CDEC_wide$Seq.Day
  LT.fx <- daylength/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace
  
  all.VBL<-all.VBWt<- array(NA,dim=c((n.days),nrow(beta_hat)))
  all.L<- data.frame(Day = NA, Montezuma = NA, RioVista = NA, s = NA)
  all.Wt<- data.frame(Day = NA, Montezuma = NA, RioVista = NA, s = NA)
  all.delta.BEM <- array(NA,dim=c((n.days),n.strata,nrow(beta_hat)))
  all.delta.VB <- array(NA,dim=c((n.days),nrow(beta_hat)))
  mn.pred.L<-mn.pred.Wt<- array(NA,dim=c((n.days+1),nrow(beta_hat)))
  energyout = data.frame(Energy = NA, Location = NA, Day = NA)

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
    L[1,] <- mean(start.L$ForkLength) # all fish are mean of starting fish lengths
    VB.L[1] <- mean(start.L$ForkLength)
    VB.Wt[1] <- mean(start.Wt)
    Wt[1,] <- mean(start.Wt)
    delta.VB[1] <- NA
    
    for (i in 1:n.strata) { 
      stratum = Locations[i]
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
      
      energyout = bind_rows(energyout, data.frame(Energy = energy.prey, 
                                                  Energypred = energy[1,i], 
                                                  Limno = Limno[1,i], Location = stratum, 
                                                  Day = 1))
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

for (i in 1:n.strata) { # region = Conf, NESuisun, LSac
      for (d in 2:n.days) { # day
        stratum = Locations[i]
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
        energyout = bind_rows(energyout, data.frame(Energy = energy[1,i], 
                                                    Energypred = energy[d,i],Limno = Limno[d,i],
                                                    Location = stratum, Day = d,
                                                    Consumption = C[d,i],
                                                    Temp = obs.temp.dat[d,i],
                                                    Turb = obs.turb.dat[d,i],
                                                    Tempeffect =  TMP.fx[d,i],
                                                    Turbeffect = NTU.fx[d,i],
                                                    allfood = sum(energy.prey)))
      }}
  
  ### 4. Get summaries
  #L[47,] #terminal lengths
  
  all.delta.BEM[,,s] <- delta.BEM
  all.delta.VB[,s] <- delta.VB
  all.L <- bind_rows(all.L, data.frame(Day = c(1:nrow(L)), Montezuma = L[,1], RioVista = L[,2], s = s))
  all.VBL[,s] <- VB.L
  all.Wt <- bind_rows(all.Wt, data.frame(Day = c(1:nrow(L)), Montezuma = Wt[,1], RioVista = Wt[,2], s = s))
  all.VBWt[,s] <- VB.Wt

}
  return(list(all.Wt, energyout))
}
  

cages2019 = cagegrowth(PD.mn.array = zoop19, CDEC_wide = filter(CDEC_wide, year(Date) == 2019 ), 
                       start.Wt, beta_hat = beta_hat[1:200,])
View(cages2019[[2]])

cages2023 = cagegrowth(PD.mn.array = zoop23, CDEC_wide = filter(CDEC_wide, year(Date) == 2023 ), 
                       start.Wt, beta_hat = beta_hat[1:200,])

cages2024 = cagegrowth(PD.mn.array = zoop24, CDEC_wide = filter(CDEC_wide, year(Date) == 2024 ), 
                       start.Wt, beta_hat = beta_hat[1:200,])


allcages = bind_rows(mutate(cages2019[[1]], Year = "2019"), mutate(cages2023[[1]], Year = "2023"),
                     mutate(cages2024[[1]], Year = "2024"))

allparams = bind_rows(mutate(cages2019[[2]], Year = "2019"), mutate(cages2023[[2]], Year = "2023"),
                     mutate(cages2024[[2]], Year = "2024"))


all.Wtlong =  pivot_longer(allcages, cols = c(RioVista, Montezuma), names_to = "Location", values_to = "Weight") %>%
  filter(!is.na(Weight))

ggplot(allparams, aes(x = Day, y = Consumption, color = Location)) +
  geom_smooth()+
  facet_wrap(~Year)

ggplot(allparams, aes(x = Day, y = allfood, color = Location)) +
  geom_smooth()+
  facet_wrap(~Year)


ggplot(allparams, aes(x = Day, y = Tempeffect, color = Location)) +
  geom_smooth()+
  facet_wrap(~Year)

ggplot(allparams, aes(x = Day, y = Turbeffect, color = Location)) +
  geom_smooth()+
  facet_wrap(~Year)



ggplot(all.Wtlong, aes(x = Day, y = Weight, color = Location)) + facet_wrap(~Year)+ geom_smooth()#+
  #geom_point(alpha = 0.2)
#This is odd, I'd expect Rio Vista to be better thanmontezuma

CDEC_wide = mutate(CDEC_wide, Year = year(Date))
ggplot(CDEC_wide, aes(x = Date, y = Temp.BDL)) + geom_line()+
  geom_line(aes(y = Temp.RVB), color = "blue") +
  facet_wrap(~Year, scales = "free")+
  ylab("Temperature")
#Rio vista is cooler in all three years

ggplot(CDEC_wide, aes(x = Date, y = NTU.BDL)) + geom_line()+
  geom_line(aes(y = NTU.RVB), color = "blue") +
  facet_wrap(~Year, scales = "free") +ylab("Turbidity")
#but Rio Vista is also much clearer, is that enough to make a difference?


############################################################
#Try replacing the cage zoops with FMWT zoops


cages2019fmwt = cagegrowth(PD.mn.array = zoop19fmwt, CDEC_wide = filter(CDEC_wide, year(Date) == 2019 ), 
                       start.Wt, beta_hat = beta_hat[1:200,])

cages2023fmwt = cagegrowth(PD.mn.array = zoop23fmwt, CDEC_wide = filter(CDEC_wide, year(Date) == 2023 ), 
                       start.Wt, beta_hat = beta_hat[1:200,])

cages2024 = cagegrowth(PD.mn.array = zoop24fmwt, CDEC_wide = filter(CDEC_wide, year(Date) == 2024 ), 
                       start.Wt, beta_hat = beta_hat[1:200,])


allcagesfmwt = bind_rows(mutate(cages2019fmwt[[1]], Year = "2019"), mutate(cages2023fmwt[[1]], Year = "2023"),
                     mutate(cages2024[[1]], Year = "2024"))

allparamsfmwt = bind_rows(mutate(cages2019fmwt[[2]], Year = "2019"), mutate(cages2023fmwt[[2]], Year = "2023"),
                      mutate(cages2024[[2]], Year = "2024"))


all.Wtlongfmwt =  pivot_longer(allcagesfmwt, cols = c(RioVista, Montezuma), names_to = "Location", values_to = "Weight") %>%
  filter(!is.na(Weight))


ggplot(all.Wtlongfmwt, aes(x = Day, y = Weight, color = Location)) + 
  facet_wrap(~Year)+ geom_smooth()#+



ggplot(allparamsfmwt, aes(x = Day, y = Consumption, color = Location)) +
  geom_smooth()+
  facet_wrap(~Year)

ggplot(allparamsfmwt, aes(x = Day, y = allfood, color = Location)) +
  geom_smooth()+
  facet_wrap(~Year)


ggplot(allparamsfmwt, aes(x = Day, y = Tempeffect, color = Location)) +
  geom_smooth()+
  facet_wrap(~Year)

ggplot(allparamsfmwt, aes(x = Day, y = Turbeffect, color = Location)) +
  geom_smooth()+
  facet_wrap(~Year)

#############################################################################################
#put the actual lengths on there
actuallenghts = data.frame(Day = c(rep(1, length(start.L$Date)), rep(47, nrow(terminal.L)) ), 
                           ForkLength = c(start.L$ForkLength, terminal.L$ForkLength),
                           Weight = c(start.Wt, terminal.Wt),
                           Location = c(rep("Start", length(start.L$Date)), terminal.L$Location)) %>%
  filter(Location != "FCCL") %>%
  mutate(Location = case_when(Location == "Rio Vista" ~ "RioVista",
                              TRUE ~ Location),
           Location = factor(Location, levels = c("Start", "RioVista", "Montezuma"))) 

#write.csv(all.Llong, "outputs/smeltcagelegnths_realinputs.csv")
oldresults = read.csv("outputs/smeltcagelegnths_realinputs.csv")


ggplot(oldresults, aes(x = Day, y = ForkLength, color = Location))+ geom_smooth()

ggplot(oldresults) + #geom_point( alpha = 0.05) + 
 geom_boxplot(data = actuallenghts, aes(x = Day, y = ForkLength, fill = Location), inherit.aes = FALSE)+
  
  stat_summary(geom = "point", data = actuallenghts, aes(x = Day, y = ForkLength, fill = Location), 
             inherit.aes = FALSE)+
geom_smooth(aes(x = Day, y = ForkLength, color = Location))+
  scale_color_manual(values = c("blue", "grey30", "red"))+
  theme_bw()+
   scale_fill_manual(values = c("red", "grey30", "blue"))


ggplot(all.Wtlong) + #geom_point( alpha = 0.05) + 
  geom_boxplot(data = actuallenghts, aes(x = Day, y = Weight, fill = Location), inherit.aes = FALSE)+
  geom_smooth(aes(x = Day, y = Weight, color = Location))+
  scale_color_manual(values = c("blue", "grey30", "red"))+
  theme_bw()+
  scale_fill_manual(values = c("blue", "grey30", "red"))


ggplot(actuallenghts, aes(x = Location, y = Weight)) + geom_boxplot()+
  stat_summary(geom = "point", color = "black", shape =1, size = 2)+
  stat_summary(data = filter(all.Wtlong, Day == 47), geom = "point", color = "blue", shape =10, size = 2)


ggplot(actuallenghts, aes(x = Location, y = ForkLength)) + geom_boxplot()
ggplot(actuallenghts, aes(x = ForkLength, fill = Location)) + geom_density(alpha = 0.5)
ggplot(actuallenghts, aes(x = Weight, fill = Location)) + geom_density(alpha = 0.5)

#OK, plot the model and the actual as box plots next to each other
lenghts_realandactual = bind_rows(filter(mutate(all.Wtlong, Type = "Modeled"),
                                         Day == 47),mutate(actuallenghts, Type = "Actual"))


ggplot(lenghts_realandactual, aes(x = Location, y = Weight, fill = Type)) + geom_boxplot()

actuallenghts %>%
  group_by(Location) %>%
  summarize(meanweight = mean(Weight),  weight_max = max(Weight),
            weight_min = min(Weight),meanlength = mean(ForkLength), FLmax = max(ForkLength), FLmin = min(ForkLength))
