#BEM funciton
library(tidyverse)
library(arrayhelpers)

smelt_bioenergetics = function(PD.mn.array, obs.temp.dat, obs.turb.dat, start.L,
                               ex.strata, beta_hat, start.year, Season = "Both"){

  ln.a<-c(0.000005,0.00000183) # length-weight parameters
  ln.b<-c(3,3.38)
  stage =4
  
  
  n.days = dim(PD.mn.array)[1]
  n.yrs = dim(PD.mn.array)[4]
  n.prey = dim(PD.mn.array)[2]
  n.strata = length(ex.strata)
  start.Wt <- ln.a[2]*start.L^ln.b[2]
  
  #I'm not going to change any of these
  VB.k<-2.72 #2.98 # VBGM parameters estimated in DSM TN 35
  Linf<-78.39 #76.1
  t0<--0.026 #-0.014
  
  
  #prey stuff
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
  
  # bioenergetics model parameters, given by Rose et al. 2013a
  a.c <- c(0.18,0.18,0.18,0.1,0.1) # Rose
  b.c <- c(-0.275,-0.275,-0.275,-0.54,-0.54)
  CQ <- c(7,7,10,10,10)						# Temperature at CK1 of maximum (deg C)
  T.L <- c(28,28,27,27,27)					# Temperature at CK4 of maximum (deg C)
  CK.1 <- c(0.4,0.4,0.4,0.4,0.4)				# effect at temperature CQ
  CK.4 <- c(0.01,0.01,0.01,0.001,0.001)			# effect at temperature T.L
  T.0 <- c(17,17,20,20,20) # Temperature at 0.98 of maximum (deg C)
  T.M <- 21.6 # Temperature above which things start delining.  smith and nobriga had this at 21.6
  #T.M <- c(20,20,21.6,21.6,21.6) #
  # Metabolism (R) parameters
  a.r <- c(0.0027,0.0027,0.0027,0.0027,0.0027)	# weight multiplier
  b.r <- c(-0.216,-0.216,-0.216,-0.216,-0.216)	# weight exponent
  #R.Q <- c(0.036,0.036,0.036,0.036,0.036)			# exponent for temperature effect
  R.Q <- c(0.064,0.064,0.064,0.064,0.064)			# updating based on Smith and Nobriga
  S.d <- c(0.175,0.175,0.175,0.175,0.175)			# Fraction of assimilated food lost to SDA (specific dynamic action)
  
  # Egestion (F) and excretion (U) parameters
  F.a <- c(0.16,0.16,0.16,0.16,0.16)				# Fraction of consumed food lost to egestion
  U.a <- c(0.1,0.1,0.1,0.1,0.1)				# Fraction of assimilated food lost to excretion
  
  e.s <- 4814									# J/g: convert g(prey)/g(delta smelt) to g(smelt)/g(smelt) - fixed
  #e.d <- c(2590,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590)		# energy density of prey items
  #e.d = c(2226, rep(3180, 11)) #updated energy densities based on Brock's analysis.
  e.d = c(1823, rep(2590, 11)) #go back to what is in Rose et al. 
  
  
  


all.L<- data.frame(NA)
all.Wt<- data.frame(NA)

daylight = daylength(38.15, c(153:310))$Daylength*60

LT.fx <- daylight/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace
  if(Season == "Summer"){
    LT.fx = LT.fx[1:60]
    daylabs = c(1:61)
  } else {
    if(Season == "Fall"){LT.fx = LT.fx[92:153]
    daylabs = c(1:63)} else daylabs = c(1:n.days)}
for(s in 1:100) {
    FOD.fx<-TMP.fx <- NTU.fx <- Cmax1 <- realized.Cmax1 <- C <- Limno <- energy <-delta.BEM <- array(NA,dim=c(n.days,n.strata,n.yrs)) # day x region x yr
  Wt <- L <- Foodsum <- array(NA,dim=c((n.days+1),n.strata,n.yrs)) # day x region x yr
  Food <- array(NA,dim=c(n.days,n.strata,n.yrs,n.prey)) # day x region x yr x prey type

  #We still want variability, but have it based on new parameters  
  logit.temp.b[1] <- runif(1,(0.43-0.25*0.43),(0.43+0.25*0.43))
  logit.temp.c[1] <- runif(1,(-2.7-0.25*2.47),(-2.7+0.25*2.7))
  logit.temp.c[2] <- runif(1,(24.44-0.25*24.44),(24.44+0.25*24.44))
  R.Q <- rep(runif(1,(0.064-0.25*0.064),(0.064+0.25*0.064)),1)
  min1 <- runif(1,(0.64-0.25*0.64),(0.64+0.25*0.64))
  b.turb <- runif(1,(0.12-0.25*0.12),(0.12+0.25*0.12))
  
 # a.c <- rep(beta_hat[s,1],5)
 # b.c <- rep(beta_hat[s,2],5)
  #T.M1 <- beta_hat[s,3]
 # R.Q <- rep(beta_hat[s,4],5)
 # min1 <- beta_hat[s,5]
  
  
  # start model
  #
  for (h in 1:n.yrs) { # year
    L[1,,h] <- start.L[h] # June start

    Wt[1,,h] <- ln.a[2]*(L[1,,h]^ln.b[2])
print(dimnames(PD.mn.array)[[4]][h])
print(h)
    
    for (i in 1:n.strata) { # region = Conf, NESuisun, LSac, Sac-Cache, SJR
      Food[1,i,h,] <- PD.mn.array[1,,ex.strata[i],h]*t(V[,4]/K[,4]) #this is selecting V and k for the fourth life stage
      Foodsum[1,i,h] <- sum(Food[1,i,h,c(3:5,9,11,12)]) #this doesn't appear to be used again....
      FOD.fx[1,i,h] <- sum(Food[1,i,h,]/(1+sum(Food[1,i,h,]))) #This is also not used gain
      
      # Temp-consumption model
      L.1<- exp((1/(T.0[stage]-CQ[stage]))*log(0.98*(1-CK.1[stage])/(CK.1[stage]*0.02))*(as.numeric(obs.temp.dat[1,i,h])-CQ[stage]))
      L.2 <- exp((1/(T.L[stage]-T.M))*log(0.98*(1-CK.4[stage])/(CK.4[stage]*0.02))*(T.L[stage]-as.numeric(obs.temp.dat[1,i,h])))
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
      
      R <- a.r[4]*(Wt[1,i,h]^b.r[4])*exp(R.Q*as.numeric(obs.temp.dat[1,i,h]))
      Feg <- F.a[4]*C[1,i,h]
      U <- U.a[4]*(C[1,i,h]-Feg)
      SDA <- S.d[4]*(C[1,i,h]-Feg)
      Wt.gain <- Wt[1,i,h]*(energy[1,i,h]/e.s)*(C[1,i,h]-R-Feg-U-SDA)
      Wt[2,i,h] <- Wt[1,i,h]+Wt.gain # add the growth
      L[2,i,h] <- (Wt[2,i,h]/ln.a[2])^(1/ln.b[2])
      delta.BEM[1,i,h] <- Wt[2,i,h]-Wt[1,i,h]
    }}
    
  
  #this is runnign the rest of the model
  for (i in 1:n.strata) { # region = Conf, NESuisun, LSac
    for (h in 1:n.yrs) { # year
      for (d in 2:n.days) { # day
        Food[d,i,h,] <- PD.mn.array[d,,ex.strata[i],h]*t(V[,4]/K[,4])
        Foodsum[d,i,h] <- sum(Food[d,i,h,c(3:5,9,11,12)])
        FOD.fx[d,i,h] <- sum(Food[d,i,h,]/(1+sum(Food[d,i,h,])))
        
        # Temp-consumption model
        L.1<- exp((1/(T.0[stage]-CQ[stage]))*log(0.98*(1-CK.1[stage])/(CK.1[stage]*0.02))*(as.numeric(obs.temp.dat[d,i,h])-CQ[stage]))
        L.2 <- exp((1/(T.L[stage]-T.M))*log(0.98*(1-CK.4[stage])/(CK.4[stage]*0.02))*(T.L[stage]-as.numeric(obs.temp.dat[d,i,h])))
        K.A <- CK.1[stage]*L.1/(1+CK.1[stage]*(L.1-1))
        K.B <- CK.4[stage]*L.2/(1+CK.4[stage]*(L.2-1))
        TMP.fx[d,i,h] <- K.A*K.B
        
        NTU.fx[d,i,h] <- min1+(1-min1)/(1+exp(-(b.turb*(as.numeric(obs.turb.dat[d,i,h])-a.turb)))) # Turbidity-consumption
        Cmax1[d,i,h] <- a.c[4]*(Wt[d,i,h]^b.c[4])
        realized.Cmax1[d,i,h] <- Cmax1[d,i,h]*TMP.fx[d,i,h]*NTU.fx[d,i,h]*LT.fx[d]
        #  
        #  #OKk percent max consumption times each food source, but why is it devided by 1+total food?
        #  C.prey <- (realized.Cmax1[d,i,h]*Food[d,i,h,]) /(1+sum(Food[d,i,h,])) # C.prey : total consumption of each prey type. function of max consumption, food 
        #  energy.prey <- e.d*C.prey
        #  #  Limno[1,i,h] <- energy.prey[1]/sum(energy.prey) # fraction energy from Limno prey
        #  #  energy[1,i,h] <- e.d[1]*Limno[1,i,h] + e.d[2]*(1-Limno[1,i,h]) #sum(energy.prey)/n.prey
        #  energyx = energy.prey/sum(energy.prey) # fraction energy from each prey
        #  energy[1,i,h] = sum(e.d*energyx)
        #  
        #  #This is a really odd way to do this. Why not just give limno a different entergy density and multiply it on the line above?
        # # Limno[d,i,h] <- energy.prey[1]/sum(energy.prey) # fraction energy from Limno prey
        # # energy[d,i,h] <- e.d[1]*Limno[d,i,h] + e.d[2]*(1-Limno[d,i,h]) #sum(energy.prey)/n.prey
        
        realized.Cmax1[d,i,h] <- Cmax1[d,i,h]*TMP.fx[d,i,h]*NTU.fx[d,i,h]*LT.fx[d]
        C.prey <- (realized.Cmax1[d,i,h]*Food[d,i,h,]) /(1+sum(Food[d,i,h,])) # C.prey : total consumption of each prey type. function of max consumption, food 
        energy.prey <- e.d*C.prey
        Limno[d,i,h] <- energy.prey[1]/sum(energy.prey) # fraction energy from Limno prey
        energy[d,i,h] <- e.d[1]*Limno[d,i,h] + e.d[2]*(1-Limno[d,i,h]) #sum(energy.prey)/n.prey
        C[d,i,h] <- sum(C.prey) # C = realized consumption rate
        
        #respiration
        R <- a.r[4]*(Wt[d,i,h]^b.r[4])*exp(R.Q*as.numeric(obs.temp.dat[d,i,h]))
        #egestion
        Feg <- F.a[4]*C[d,i,h]
        #excretion
        U <- U.a[4]*(C[d,i,h]-Feg)
        #specific dynamic action
        SDA <- S.d[4]*(C[d,i,h]-Feg)
        #weght gain
        Wt.gain <- Wt[d,i,h]*(energy[d,i,h]/e.s)*(C[d,i,h]-R-Feg-U-SDA)
        Wt[(d+1),i,h] <- Wt[d,i,h]+Wt.gain # add the growth
      #  if(Wt.gain <0) {
      #    L[(d+1),i,h]  = L[(d),i,h]
      #  } else L[(d+1),i,h] <- (Wt[(d+1),i,h]/ln.a[2])^(1/ln.b[2])
        
        delta.BEM[d,i,h] <- Wt[(d+1),i,h]-Wt[d,i,h]
      }}}
  
  ### 4. Get summaries

  AllL = array2df(L, levels = list(NA,ex.strata,c(1:n.yrs)), label.x = "Length", na.rm =T) %>%
    mutate(s = s)
  all.L <- bind_rows(all.L, AllL) 
  
  all.Wt <- bind_rows(all.Wt, mutate(array2df(Wt,levels = list(NA,ex.strata,c(1:n.yrs)),  label.x = "Weight", na.rm =T),s = s))
}
lengthweight = left_join(all.Wt, all.L)

Outputs = select(lengthweight, -NA.) %>%
  filter(!is.na(d1)) %>%
  rename(Day = d1 ,
         Stratum = d2,
         Year = d3) %>%
  mutate(Day = as.numeric(as.character(factor(Day))),
         Year = as.numeric(as.character(factor(Year, labels = c((start.year):(start.year+n.yrs-1))))))
        

return(Outputs)
}


