#BEM funciton
library(tidyverse)
library(arrayhelpers)

smelt_bioenergetics = function(PD.mn.array, obs.temp.dat, obs.turb.dat, start.L,
                               ex.strata, beta_hat, start.year){
  n.days = dim(PD.mn.array)[1]
  n.years = dim(PD.mn.array)[4]
  n.prey = dim(PD.mn.array)[2]
  n.strata = length(ex.strata)
  start.Wt <- ln.a[2]*start.L^ln.b[2]

all.VBL<-all.VBWt<- array(NA,dim=c((n.days), n.yrs,nrow(beta_hat)))
all.L<- data.frame(NA)
all.Wt<- data.frame(NA)
all.delta.BEM <- array(NA,dim=c((n.days),n.strata, n.yrs,nrow(beta_hat)))
all.delta.VB <- array(NA,dim=c((n.days), n.yrs,nrow(beta_hat)))
mn.pred.L<-mn.pred.Wt<- array(NA,dim=c((n.days+1), (n.yrs-1),nrow(beta_hat)))

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
      Foodsum[1,i,h] <- sum(Food[1,i,h,c(3:5,9,11,12)]) #this doesn't appear to be used again....
      FOD.fx[1,i,h] <- sum(Food[1,i,h,]/(1+sum(Food[1,i,h,]))) #This is also not used gain
      
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
        R <- a.r[4]*(Wt[d,i,h]^b.r[4])*exp(R.Q[stage]*as.numeric(obs.temp.dat[d,i,h]))
        #egestion
        Feg <- F.a[4]*C[d,i,h]
        #excretion
        U <- U.a[4]*(C[d,i,h]-Feg)
        #specific dynamic action
        SDA <- S.d[4]*(C[d,i,h]-Feg)
        #weght gain
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
  # all.L <- bind_rows(all.L, data.frame(Day = c(152:305), 
  #                                      `NE Suisun` = L[,3,], 
  #                                      `SE Suisun` = L[,5,],
  #                                      `NW Suisun` = L[,4,],
  #                                      `Lower Sacramento River` = L[,2,],
  #                                      `Confluence` = L[,1,],
  #                                      `SW Suisun` = L[,7,],
  #                                      `Suisun Marsh` = L[,6,],
  #                                      Year = h, s = s))
  # all.VBL[,,s] <- VB.L
  # all.Wt <- bind_rows(all.Wt, data.frame(Day = c(152:305), 
  #                                        `NE Suisun` = Wt[,3,], 
  #                                        `SE Suisun` = Wt[,5,],
  #                                        `NW Suisun` = Wt[,4,],
  #                                        `Lower Sacramento River` = Wt[,2,],
  #                                        `Confluence` = Wt[,1,],
  #                                        `SW Suisun` = Wt[,7,],
  #                                        `Suisun Marsh` = Wt[,6,],
  #                                        Yr = h, s = s))
  
  AllL = array2df(L, levels = list(NA,ex.strata,c(2:n.yrs)), label.x = "Length", na.rm =T) %>%
    mutate(s = s)
  all.L <- bind_rows(all.L, AllL) 
  all.VBL[,,s] <- VB.L
  all.Wt <- bind_rows(all.Wt, mutate(array2df(Wt,levels = list(NA,ex.strata,c(2:n.yrs)),  label.x = "Weight", na.rm =T), s = s))
  
  all.VBWt[,,s] <- VB.Wt
  
}
lengthweight = left_join(all.Wt, all.L)

Outputs = select(lengthweight, -NA.) %>%
  filter(!is.na(d1)) %>%
  rename(Day = d1 ,
         Stratum = d2,
         Year = d3) %>%
  mutate(Day = as.numeric(as.character(factor(Day, labels = c(0:n.days)))),
         Year = as.numeric(as.character(factor(Year, labels = c((start.year+1):(start.year+n.yrs-1))))))

return(Outputs)
}


