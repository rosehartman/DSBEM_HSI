### Plots and models for bioenergetics model, consumption report ###
# William Smith, 19 Feb 2021
# 1. Read model parameters and data
# 2. Make a plot of Cmax, R, Cmax-temp.fx, and Cmax-turb.fx
# 3. Model bioenergetics and reference points
# 4. Get summaries
# 5. Graphs

library(abind)
library(zoo)
library(lubridate)
library(lattice)
library(ggplot2)
library(gridExtra)
library(paletteer)
library(latticeExtra)
library(gsw)

#setwd("D:/FWS/R code/Rose BEM/Peterson-Smith-Rose BEM, 2020 revisions")

### 1. Read model parameters and data
n.strata <- 13
n.prey <- 12
n.days<-153
ln.a<-c(0.000005,0.00000183) # length-weight parameters
ln.b<-c(3,3.38)
startY <- 2015
endY <- 2020
yr.seq<-seq(startY-1994,endY-1994,by=1) # year 1 = 1995
n.yrs<-length(yr.seq)

VB.k<-2.72 #2.98 # VBGM parameters estimated in DSM TN 35
Linf<-78.39 #76.1
t0<--0.026 #-0.014
start.L <- c(24.3,29.2,22.2,19.9,21.9,28.1) # assumed length on Jun 1, 1st value = EDSM+SKT observed, (5/15-6/15),2017,2018,2019
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

source("BEM consumption functions_v3.2.R") # load consumption functions
CDEC <- read.table(file="CDEC_temp_turb.txt",header=T)
X2 <- read.table(file="X2_daily.txt",header=T)
obs.temp.dat <- make.temp()
obs.turb.dat <- make.turb()
daylength<-read.table(file="daylength.txt")
day <- CDEC[,2]
LT.fx <- daylength/max(daylength)
wtr.yr<-c(1,1,1,1,1,1,4,4,2,3,3,1,4,5,4,3,1,3,4,5,5) # Sacto WY type wet = 1,above normal = 2, below normal = 3, dry = 4, critical = 5
rep.yrs <- seq(21,26,by=1) # years 2015-2020 from 1995-2020 dataset, for prey density summary
ex.strata <- c(13,10,8,6,3,1) # map IBMR strata to 5 CDEC regions

beta_hat <- read.table(file='beta_hat.txt') # MC filtering coefficients
#this is model uncertainty
#beta_hat <- read.table(file='BEM consumption rate project/beta_hat_IND.txt') # MC filtering coefficients

### 2. Make a plot of Cmax-L, R-temp, Cmax-temp.fx, and Cmax-turb.fx
# Cmax and R
Len1 <- seq(25,75,by=1)
Len2 <- c(20,75)
temp.dat <- seq(14,27,length.out=100) # make a temp vector for plotting predictions (18C is approx. minimum of CDEC temps)
NTU.plot<-seq(1,100,length.out=100) # NTU=100 is approx. max of revised CDEC data

W1 <- ln.a[2]*(Len1^ln.b[2])
W2 <- ln.a[2]*(Len2^ln.b[2])
Cmax <- vector()
R.Rose<-R <- matrix(NA,length(W2),length(temp.dat))
for (i in 1:length(W1)) {
 Cmax[i] <- mean(beta_hat[,1]*(W1[i]^beta_hat[,2]))
 }
for (i in 1:length(W2)) {
for (j in 1:length(temp.dat)) {
 R[i,j] <- mean(a.r[4]*(W2[i]^b.r[4])*exp(beta_hat[,4]*temp.dat[j]))
 R.Rose[i,j] <- 0.0027*(W2[i]^b.r[4])*exp(0.036*temp.dat[j])
 }}
 
# Cmax-temp.fx and -turb.fx
temp.fx.out<-temp.fx.out2<-turb.fx.out <- vector()
for (i in 1:length(temp.dat)) {
 temp.fx.out[i] <- mean(temp.fx(new.temp.dat<-temp.dat[i],stage<-4,T.M1<-beta_hat[,3]))
 temp.fx.out2[i] <- temp.fx(new.temp.dat<-temp.dat[i],stage<-4,T.M1<-23)
 turb.fx.out[i] <- mean(beta_hat[,5]+(1-beta_hat[,5])/(1+exp(-(b.turb*(NTU.plot[i]-a.turb)))))
 } 

# plots
par(mfrow=c(3,2))
par(family='serif',mar = c(4,4,2,2), oma = c(1,1.1,1,1),cex.lab=1.5,cex.axis=1.3)
plot(W1,Cmax,ylab="Cmax",xlab="Weight (g)",type="l",ylim=c(0,0.30))
lines(W1,a.c[stage]*W1^b.c[stage],col='red')
legend('topright',legend=c("Current model","Rose et al. (2013a) model"),lty=c("solid","solid"),col=c("black","red"),cex=1.15)
text(0.25,0.293,'A',cex=1.5)

plot(temp.dat,R[1,],ylab="Respiration",xlab="Temperature (C)",ylim=c(0,1.5*max(R)),type="l")
lines(temp.dat,R[2,],lty="dashed")
lines(temp.dat,R.Rose[1,],col='red')
lines(temp.dat,R.Rose[2,],col='red',lty="dashed")
#rug(obs.temp.dat, ticksize = 0.02, side = 1, lwd = 0.01)
legend('topright',legend=c("25 mm FL","45 mm FL"),lty=c("solid","dashed","solid","dashed"),col=c("black","black"),cex=1.15)
text(1.025*min(temp.dat),1.5*max(R)-0.0007,'B',cex=1.5)

plot(temp.dat,temp.fx.out,ylab='Temperature effect',xlab='Temperature (C)',type='l',ylim=c(0,1))
lines(temp.dat,temp.fx.out2,col='red')
#rug(obs.temp.dat, ticksize = 0.02, side = 1, lwd = 0.01)
text(1.025*min(temp.dat),0.955,'C',cex=1.5)

plot(NTU.plot,turb.fx.out,type='l',xlab='Turbidity (NTU)',ylab='Turbidity effect',ylim=c(0,1))
#rug(obs.turb.dat, ticksize = 0.02, side = 1, lwd = 0.01)
text(3+min(NTU.plot),0.955,'D',cex=1.5)

plot(seq(1,n.days,by=1),LT.fx[(1:n.days),1],ylab="Daylight effect",xlab="Days since May 31",type="l",ylim=c(0,1))
text(5,0.95,'E',cex=1.5)

# Prey selection plots
Food.mn <- c(0.502,0.132,2.10,1.26,1.74,0.141,0,0.258,0.186,0.0048,0.759,8.94)
Food.lo <- rep(0,12) #c(0.0043573129,0.0150948662,0.0704961354,0.1889142921,0.207291051,0.0036499878,0.0000000000,0.0038854292,0.0205220011,0.0009358132,0.0311836170,0.6078908852)
Food.hi <- c(0.702396111,0.216322882,3.170230321,1.644875925,2.576003187,0.137358171,0.000000000,0.134133713,0.183127854,0.004932837,0.819951747,12.690013802 )

EUT <- seq(0.003,5,length.out=100)
DAP <- seq(0.003,5,length.out=100)
food.plot.EUT<-food.plot.DAP <- matrix(NA,length(EUT),n.prey)

for (i in 1:length(EUT)) {
 for (j in c(1:7,9:n.prey)) {
  food.plot.DAP[i,j] <- (Food.lo[j]*V[j,4]/K[j,4])/(1+sum(c(Food.lo[1:7],DAP[i],Food.lo[9:n.prey])*V[,4]/K[,4]))
  }
 food.plot.DAP[i,8] <- (DAP[i]*V[8,4]/K[8,4])/(1+sum(c(Food.lo[1:7],DAP[i],Food.lo[9:n.prey])*V[,4]/K[,4]))
 for (j in c(1:10,n.prey)) {
  food.plot.EUT[i,j] <- (Food.lo[j]*V[j,4]/K[j,4])/(1+sum(c(Food.lo[1:10],EUT[i],Food.lo[n.prey])*V[,4]/K[,4]))
  }
 food.plot.EUT[i,11] <- (EUT[i]*V[11,4]/K[11,4])/(1+sum(c(Food.lo[1:10],EUT[i],Food.lo[n.prey])*V[,4]/K[,4]))
 }
C.DAP <- rowSums(food.plot.DAP)
C.EUT <- rowSums(food.plot.EUT)

xlab.text = expression(paste("Food (mg Carbon/m³)"))

plot(EUT,C.EUT,type='l',ylab='Prey effect',xlab=xlab.text,ylim=c(0,1))
lines(DAP,C.DAP,lty='dashed')
legend('bottomright',legend=c("Eurytemora","Daphnia"),lty=c("solid","dashed"),col=c("black","black"),cex=1.15)
text(0.09,0.955,'F',cex=1.5)

### 3. Model bioenergetics and reference points
PD.array<-array(NA,dim=c(length(rep.yrs),365,n.strata,n.prey)) # get prey densities
for (i in 1:length(rep.yrs)) {
 PD.array[i,,,] <- make.food(rep.yrs[i]) # year x day x strata x prey type
 }
 
PD.mn.array<-array(NA,dim=c(n.yrs,n.days,n.strata,n.prey)) # Make into smaller matrix for modeling
for (y in 1:n.yrs) {
for (j in 1:n.strata) {
for (k in 1:n.prey) {
 PD.mn.array[y,,j,k] <- PD.array[y,151:(151+n.days-1),j,k]
 }}}
for (d in 139:n.days) {
 PD.mn.array[6,d,ex.strata[2],] <- PD.mn.array[6,138,ex.strata[2],] # Honker Bay prey data are missing for the end of 2020
}

#what is the average zooplankton density here?

PD.array[1,100,1,]


all.VBL<-all.VBWt<- array(NA,dim=c((n.days),n.yrs,nrow(beta_hat)))
all.L<-all.Wt<- array(NA,dim=c((n.days+1),length(ex.strata),n.yrs,nrow(beta_hat)))
all.delta.BEM <- array(NA,dim=c((n.days),length(ex.strata),n.yrs,nrow(beta_hat)))
all.delta.VB <- array(NA,dim=c((n.days),n.yrs,nrow(beta_hat)))
mn.pred.L<-mn.pred.Wt<- array(NA,dim=c((n.days+1),nrow(beta_hat),(n.yrs-1)))

for (s in 1:nrow(beta_hat)) {
VB.L <- VB.Wt <- delta.VB <- array(NA,dim=c(n.days,n.yrs))
FOD.fx<-TMP.fx <- NTU.fx <- Cmax1 <- realized.Cmax1 <- C <- Limno <- energy <-delta.BEM <- array(NA,dim=c(n.days,length(ex.strata),n.yrs)) # day x region x yr
Wt <- L <- Foodsum <- array(NA,dim=c((n.days+1),length(ex.strata),n.yrs)) # day x region x yr
Food <- array(NA,dim=c(n.days,length(ex.strata),n.yrs,n.prey)) # day x region x yr x prey type

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

for (i in 1:length(ex.strata)) { # region = Conf, NESuisun, LSac, Sac-Cache, SJR
 Food[1,i,h,] <- PD.mn.array[h,1,ex.strata[i],]*t(V[,4]/K[,4]) #this is selecting V and k for the fourth life stage
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
 realized.Cmax1[1,i,h] <- Cmax1[1,i,h]*TMP.fx[1,i,h]*NTU.fx[1,i,h]*LT.fx[1,1]
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
for (i in 1:length(ex.strata)) { # region = Conf, NESuisun, LSac
for (h in 2:n.yrs) { # year
for (d in 2:n.days) { # day
 Food[d,i,h,] <- PD.mn.array[h,d,ex.strata[i],]*t(V[,4]/K[,4])
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
 realized.Cmax1[d,i,h] <- Cmax1[d,i,h]*TMP.fx[d,i,h]*NTU.fx[d,i,h]*LT.fx[d,1]
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

for (d in 1:(n.days+1)) {
mn.pred.L[d,s,] <- c(mean(L[d,3:6,2],na.rm=T),
					mean(L[d,2:6,3],na.rm=T),
					mean(L[d,3:6,4],na.rm=T),
					mean(L[d,1:6,5],na.rm=T),
					mean(L[d,3:6,6],na.rm=T))
mn.pred.Wt[d,s,] <- c(mean(Wt[d,3:6,2],na.rm=T),
					mean(Wt[d,2:6,3],na.rm=T),
					mean(Wt[d,3:6,4],na.rm=T),
					mean(Wt[d,1:6,5],na.rm=T),
					mean(Wt[d,3:6,6],na.rm=T))
 }

all.delta.BEM[,,,s] <- delta.BEM
all.delta.VB[,,s] <- delta.VB
all.L[,,,s] <- L
all.VBL[,,s] <- VB.L
all.Wt[,,,s] <- Wt
all.VBWt[,,s] <- VB.Wt
}

# ratio of predicted to reference growth
ratio<-array(NA,dim=c(n.days,length(ex.strata),n.yrs))
for (i in 1:length(ex.strata)) {
for (j in 1:n.yrs) {
 for (t in 1:n.days-1) {
  mn <- mean(all.delta.BEM[t,i,j,],na.rm=T)
  ratio[t,i,j] <- mn/delta.VB[t+1,j]
  }
 }}

### 5. Graphs
# plot food, NTU.fx, TMP.fx, and Light.fx time series
thresh<-0.75 # get foodsum at threshold of Food effect for plot
i<-1
j<-1
Foodref <- Foodsum[which(round(FOD.fx[,i,j],digits=2)==thresh),i,j]
for (j in 2:n.yrs) {
 Foodref <- c(Foodref,Foodsum[which(round(FOD.fx[,i,j],digits=2)==thresh),i,j])
 }
for (i in 2:n.yrs) {
for (j in 2:n.yrs) {
 Foodref <- c(Foodref,Foodsum[which(round(FOD.fx[,i,j],digits=2)==thresh),i,j])
 }}
Foodref.mn <- mean(Foodref)


graph.yr<-c(2,3,4,5,6)
par(mfrow=c(3,5))
par(family='serif',mar = c(1,1,0.92,0.1), oma = c(3,3.2,0.9,3.5),cex.lab=1.75,cex.axis=1.75)
colscheme  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Food plot
plot(seq(1,n.days,by=1),Foodsum[(1:n.days),1,graph.yr[1]],ylab="",xlab="",type="l",ylim=c(0.75,86),xaxt="n",lty="solid",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),2,graph.yr[1]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),3,graph.yr[1]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),4,graph.yr[1]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),5,graph.yr[1]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),6,graph.yr[1]],lty="twodash",col=colscheme[4])
abline(Foodref.mn,0,lty='dotted')
legend(50,83,legend=c("NW Suisun","NE Suisun","Confluence","Lower Sac","Sac-Cache","DWSC"),lty=c("solid","dotted","dashed","dotdash","longdash","twodash"),lwd=rep(1.25,6),col=c(colscheme[7],colscheme[7],colscheme[3],colscheme[3],colscheme[4],colscheme[4]),cex=1.5,bty='n')
text(65,86,'2016 - Below Normal',cex=1.8)
ylab.text = expression(paste("Food (mg Carbon/m³)"))
mtext(ylab.text, side=2,line=2.07,cex=1.3)

plot(seq(1,n.days,by=1),Foodsum[(1:n.days),1,graph.yr[2]],ylab="",xlab="",type="l",ylim=c(0.75,86),xaxt="n",yaxt="n",lty="solid",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),2,graph.yr[2]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),3,graph.yr[2]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),4,graph.yr[2]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),5,graph.yr[2]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),6,graph.yr[2]],lty="twodash",col=colscheme[4])
abline(Foodref.mn,0,lty='dotted')
text(33,86,'2017 - Wet',cex=1.8)
plot(seq(1,n.days,by=1),Foodsum[(1:n.days),1,graph.yr[3]],ylab="",xlab="",type="l",ylim=c(0.75,86),xaxt="n",yaxt="n",lty="solid",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),2,graph.yr[3]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),3,graph.yr[3]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),4,graph.yr[3]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),5,graph.yr[3]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),6,graph.yr[3]],lty="twodash",col=colscheme[4])
abline(Foodref.mn,0,lty='dotted')
text(65,86,'2018 - Below Normal',cex=1.8)
plot(seq(1,n.days,by=1),Foodsum[(1:n.days),1,graph.yr[4]],ylab="",xlab="",type="l",ylim=c(0.75,86),xaxt="n",yaxt="n",lty="solid",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),2,graph.yr[4]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),3,graph.yr[4]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),4,graph.yr[4]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),5,graph.yr[4]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),6,graph.yr[4]],lty="twodash",col=colscheme[4])
abline(Foodref.mn,0,lty='dotted')
text(33,86,'2019 - Wet',cex=1.8)
plot(seq(1,n.days,by=1),Foodsum[(1:n.days),1,graph.yr[5]],ylab="",xlab="",type="l",ylim=c(0.75,86),xaxt="n",yaxt="n",lty="solid",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),2,graph.yr[5]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),3,graph.yr[5]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),4,graph.yr[5]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),Foodsum[(1:n.days),5,graph.yr[5]],lty="longdash",col=colscheme[4])
lines(seq(1,75,by=1),Foodsum[(1:75),6,graph.yr[5]],lty="twodash",col=colscheme[4])
abline(Foodref.mn,0,lty='dotted')
text(32,85,'2020 - Dry',cex=1.8)

# Turbdiity plot
plot(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),1,graph.yr[1]],ylab="",xlab="",type="l",xaxt="n",lty="solid",col=colscheme[7],ylim=c(0,165))
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),2,graph.yr[1]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),3,graph.yr[1]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),4,graph.yr[1]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),5,graph.yr[1]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),6,graph.yr[1]],lty="twodash",col=colscheme[4])
#mtext("Turbidity effect", side=2,line=2.3,cex=1.1)
mtext("Turbidity (NTU)", side=2, line=2.37,cex=1.3)
abline(9,0,lty='dotted')

plot(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),1,graph.yr[2]],ylab="",xlab="",type="l",xaxt="n",yaxt="n",lty="solid",col=colscheme[7],ylim=c(0,165))
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),2,graph.yr[2]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),3,graph.yr[2]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),4,graph.yr[2]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),5,graph.yr[2]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),6,graph.yr[2]],lty="twodash",col=colscheme[4])
abline(9,0,lty='dotted')

plot(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),1,graph.yr[3]],ylab="",xlab="",type="l",xaxt="n",yaxt="n",lty="solid",col=colscheme[7],ylim=c(0,165))
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),2,graph.yr[3]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),3,graph.yr[3]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),4,graph.yr[3]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),5,graph.yr[3]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),6,graph.yr[3]],lty="twodash",col=colscheme[4])
abline(9,0,lty='dotted')

plot(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),1,graph.yr[4]],ylab="",xlab="",type="l",xaxt="n",yaxt="n",lty="solid",col=colscheme[7],ylim=c(0,165))
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),2,graph.yr[4]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),3,graph.yr[4]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),4,graph.yr[4]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),5,graph.yr[4]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),6,graph.yr[4]],lty="twodash",col=colscheme[4])
abline(9,0,lty='dotted')

plot(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),1,graph.yr[5]],ylab="",xlab="",type="l",xaxt="n",yaxt="n",lty="solid",col=colscheme[7],ylim=c(0,165))
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),2,graph.yr[5]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),3,graph.yr[5]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),4,graph.yr[5]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),5,graph.yr[5]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.turb.dat[(1:n.days),6,graph.yr[5]],lty="twodash",col=colscheme[4])
abline(9,0,lty='dotted')

# Temp plot
plot(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),1,graph.yr[1]],ylab="",xlab="",type="l",lty="solid",col=colscheme[7],ylim=c(14,26))
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),2,graph.yr[1]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),3,graph.yr[1]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),4,graph.yr[1]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),5,graph.yr[1]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),6,graph.yr[1]],lty="twodash",col=colscheme[4])
mtext("Temperature (C)", side=2, line=2.37,cex=1.3)
abline(22.8,0,lty='dotted')

plot(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),1,graph.yr[2]],ylab="",xlab="",type="l",yaxt="n",lty="solid",col=colscheme[7],ylim=c(14,26))
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),2,graph.yr[2]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),3,graph.yr[2]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),4,graph.yr[2]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),5,graph.yr[2]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),6,graph.yr[2]],lty="twodash",col=colscheme[4])
abline(22.8,0,lty='dotted')

plot(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),1,graph.yr[3]],ylab="",xlab="",type="l",yaxt="n",lty="solid",col=colscheme[7],ylim=c(14,26))
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),2,graph.yr[3]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),3,graph.yr[3]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),4,graph.yr[3]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),5,graph.yr[3]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),6,graph.yr[3]],lty="twodash",col=colscheme[4])
abline(22.8,0,lty='dotted')

plot(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),1,graph.yr[4]],ylab="",xlab="",type="l",yaxt="n",lty="solid",col=colscheme[7],ylim=c(14,26))
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),2,graph.yr[4]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),3,graph.yr[4]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),4,graph.yr[4]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),5,graph.yr[4]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),6,graph.yr[4]],lty="twodash",col=colscheme[4])
abline(22.8,0,lty='dotted')

plot(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),1,graph.yr[5]],ylab="",xlab="",type="l",yaxt="n",lty="solid",col=colscheme[7],ylim=c(14,26))
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),2,graph.yr[5]],lty="dotted",col=colscheme[7])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),3,graph.yr[5]],lty="dashed",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),4,graph.yr[5]],lty="dotdash",col=colscheme[3])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),5,graph.yr[5]],lty="longdash",col=colscheme[4])
lines(seq(1,n.days,by=1),obs.temp.dat[(1:n.days),6,graph.yr[5]],lty="twodash",col=colscheme[4])
abline(22.8,0,lty='dotted')

mtext("Days since May 31", side=1,line=1.5,outer=T,cex=1.3)


# Make heatmap plots of growth potential ratio
area<-cbind(rep(1,(n.days-1)),
			rep(2,(n.days-1)),
			rep(3,(n.days-1)),
			rep(4,(n.days-1)),
			rep(5,(n.days-1)),
			rep(6,(n.days-1)))
day.seq<-seq(1,(n.days-1),by=1)
day<-day.seq
for (i in 2:n.yrs) {
 day <- cbind(day,day.seq)
 }

ckey <- list(labels=list(cex=1))
margin <-  theme(plot.margin=unit(c(0,0,0,0),"cm"))
par.main.text = list(font = 2, # make it bold
                       just = "left", 
                       x = grid::unit(25, "mm"),
					   y = grid::unit(25, "mm"))


regions <- c('NW Suisun','NE Suisun','Confluence','Lower Sac','Sac-Cache','Yolo')


p11 <- levelplot(ratio[(1:(n.days-1)),,2]~day*area,
ylab=list(label='',cex=3),xlab='',zlab=list('Predicted/reference growth',cex=1),
col.regions=c(paletteer_c("grDevices::Sunset", 100)[100],
 paletteer_c("grDevices::Sunset", 100)[100:1],
 paletteer_c("grDevices::Sunset", 100)[1]),at=c(1.84,seq(1,0.5,length.out=100),-0.26),
#colorkey=list(ckey,title=expression('Predicted/reference growth'), row=3, column=1, hjust=1, vjust=1,cex=0.5),
main=list('2016',side=1,line=5),
scales=list(z=list(cex=1),x=list(at=NULL,cex=1),y=list(cex=1)),
par.settings=list(layout.heights=list(top.padding=-5, bottom.padding=-5.5),superpose.line = list(col=c("black"))))
p12 <- lattice::xyplot(X2[1:(n.days-1),4] ~ day[,1], type="l",ylim=c(52,123),col='black',
                       lwd=2,par.settings = list(superpose.line = list(col=c("black"))))
p13 <- lattice::xyplot(sal.bound[1:(n.days-1),2,2] ~ day[,1], type="l",ylim=c(0.5,6.5),col='black',
                       lwd=2,par.settings = list(superpose.line = list(col=c("black"))),
 scales=list(y = list(draw = FALSE)))

p21 <- levelplot(ratio[(1:(n.days-1)),,3]~day*area,
ylab=list(label='',cex=3),xlab='',zlab=list('Predicted/reference growth',cex=1),
col.regions=c(paletteer_c("grDevices::Sunset", 100)[100],
 paletteer_c("grDevices::Sunset", 100)[100:1],
 paletteer_c("grDevices::Sunset", 100)[1]),at=c(1.84,seq(1,0.5,length.out=100),-0.26),
colorkey=F,scales=list(z=list(cex=1),x=list(at=NULL,cex=1),y=list(cex=1)),
par.settings=list(layout.heights=list(top.padding=-6, bottom.padding=-3),layout.widths=list(right.padding=7.7),superpose.line = list(col=c("black"))))
p22 <- lattice::xyplot(X2[1:(n.days-1),5] ~ day[,1], type="l",ylim=c(52,123),col='black',lwd=2,par.settings = list(superpose.line = list(col=c("black"))))
p23 <- lattice::xyplot(sal.bound[1:(n.days-1),2,3] ~ day[,1], type="l",ylim=c(0.5,6.5),col='black',lwd=2,par.settings = list(superpose.line = list(col=c("black"))),
 scales=list(y = list(draw = FALSE)))

p31 <- levelplot(ratio[(1:(n.days-1)),,4]~day*area,
ylab=list(label='',cex=3),xlab='',zlab=list('Predicted/reference growth',cex=1),
col.regions=c(paletteer_c("grDevices::Sunset", 100)[100],
 paletteer_c("grDevices::Sunset", 100)[100:1],
 paletteer_c("grDevices::Sunset", 100)[1]),at=c(1.84,seq(1,0.5,length.out=100),-0.26),
colorkey=F,scales=list(z=list(cex=1),x=list(at=NULL,cex=1),y=list(cex=1)),
par.settings=list(layout.heights=list(top.padding=-8.5, bottom.padding=-2.3),
                  layout.widths=list(right.padding=7.7),superpose.line = list(col=c("black"))))
p32 <- lattice::xyplot(X2[1:(n.days-1),6] ~ day[,1], type="l",ylim=c(52,123),col='black',lwd=2,
                       par.settings = list(superpose.line = list(col=c("black"))))
p33 <- lattice::xyplot(sal.bound[1:(n.days-1),2,4] ~ day[,1], type="l",ylim=c(0.5,6.5),
                       col='black',lwd=2,par.settings = list(superpose.line = list(col=c("black"))),
 scales=list(y = list(draw = FALSE)))

p41 <- levelplot(ratio[(1:(n.days-1)),,5]~day*area,
ylab=list(label='',cex=3),xlab='',zlab=list('Predicted/reference growth',cex=1),
col.regions=c(paletteer_c("grDevices::Sunset", 100)[100],
 paletteer_c("grDevices::Sunset", 100)[100:1],
 paletteer_c("grDevices::Sunset", 100)[1]),at=c(1.84,seq(1,0.5,length.out=100),-0.26),
colorkey=F,scales=list(z=list(cex=1),x=list(at=NULL,cex=1),y=list(cex=1)),
par.settings=list(layout.heights=list(top.padding=-9.075, bottom.padding=-2.25),layout.widths=list(right.padding=7.7),superpose.line = list(col=c("black"))))
p42 <- lattice::xyplot(X2[1:(n.days-1),7] ~ day[,1], type="l",ylim=c(52,123),col='black',lwd=2,par.settings = list(superpose.line = list(col=c("black"))))
p43 <- lattice::xyplot(sal.bound[1:(n.days-1),2,5] ~ day[,1], type="l",ylim=c(0.5,6.5),col='black',lwd=2,par.settings = list(superpose.line = list(col=c("black"))),
 scales=list(y = list(draw = FALSE)))

p51 <- levelplot(ratio[(1:(n.days-1)),,6]~day*area,
ylab=list(label='',cex=3),xlab='',zlab=list('Predicted/reference growth',cex=1),
col.regions=c(paletteer_c("grDevices::Sunset", 100)[100],
 paletteer_c("grDevices::Sunset", 100)[100:1],
 paletteer_c("grDevices::Sunset", 100)[1]),at=c(1.84,seq(1,0.5,length.out=100),-0.26),
colorkey=F,scales=list(z=list(cex=1),x=list(cex=1),y=list(cex=1)),
par.settings=list(layout.heights=list(top.padding=-9.15, bottom.padding=-0.75),layout.widths=list(right.padding=7.7),superpose.line = list(col=c("black"))))
p52 <- lattice::xyplot(X2[1:(n.days-1),8] ~ day[,1], type="l",ylim=c(52,123),col='black',lwd=2,par.settings = list(superpose.line = list(col=c("black"))))
p53 <- lattice::xyplot(sal.bound[1:(n.days-1),2,6] ~ day[,1], type="l",ylim=c(0.5,6.5),col='black',lwd=2,par.settings = list(superpose.line = list(col=c("black"))),
 scales=list(y = list(draw = FALSE)))

grid.arrange(doubleYScale(p11,p13),doubleYScale(p21,p23),doubleYScale(p31,p33),doubleYScale(p41,p43),doubleYScale(p51,p53),nrow=5)

mtext('Region',side=2,line=1,cex=1.5,outer=T)
mtext('Days since May 31',side=1,line=1.5,cex=1.5,outer=T)

# Make observed length vs predicted BEM growth trajectory plot
med.L<-lo.L<-hi.L<-med.Wt<-lo.Wt<-hi.Wt <- array(NA,dim=c((n.days+1),5))

for (j in 1:5) {
for (t in 1:(n.days+1)) {
 med <- mean(mn.pred.L[t,,j],na.rm=T)
 lo <- quantile(mn.pred.L[t,,j],0.025,na.rm=T)
 hi <- quantile(mn.pred.L[t,,j],0.975,na.rm=T)
 med.L[t,j] <- med
 lo.L[t,j] <- lo
 hi.L[t,j] <- hi
 }}
for (j in 1:5) {
for (t in 1:(n.days+1)) {
 med <- mean(mn.pred.Wt[t,,j],na.rm=T)
 lo <- quantile(mn.pred.Wt[t,,j],0.025,na.rm=T)
 hi <- quantile(mn.pred.Wt[t,,j],0.975,na.rm=T)
 med.Wt[t,j] <- med
 lo.Wt[t,j] <- lo
 hi.Wt[t,j] <- hi
 }}

par(mfrow=c(2,3))
par(family='serif',mar = c(1,1,0.5,0.5), oma = c(3,3.5,1,1),cex.lab=1.5,cex.axis=1.3)

for (h in 1:(n.yrs-1)) {
plot(seq(1,n.days,by=1),med.Wt[(2:(n.days+1)),h],type='l',ylim=c(0.04,2.25),xaxt='n',lwd=2,xlab='',ylab='',cex.axis=1.5)
lines(seq(1,n.days,by=1),rowMeans(all.VBWt[(1:n.days),(h+1),],na.rm=T),lty='dotted',col='red')
lines(seq(1,n.days,by=1),lo.Wt[(2:(n.days+1)),h],col='gray')
lines(seq(1,n.days,by=1),hi.Wt[(2:(n.days+1)),h],col='gray')
abline(terminal.Wt[h+1],0,lty='dotted')
title(paste0('20',h+15),line=-1.5,cex.main=2)
}

frame()
legend('bottomright',legend=c("Reference weight","Terminal observed weight","Observed weight +/- se","Mean BEM prediction","95% CI BEM prediction"),lty=c("dotted","dotted","dotted","solid","solid"),col=c('red','black','gray','black','gray'),lwd=c(1,1,1,2,1),cex=1.75)


mtext('Weight (g)',side=2.1,outer=T,line=1.3,cex=1.25)
mtext('Days since May 31',side=1,outer=T,line=1.5,cex=1.25)
