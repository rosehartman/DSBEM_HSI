### Functions for bioenergetics model, consumption report ###
# William Smith, 9 Aug 2022
# this script contains several functions to process data and estimate consumption and respiration rates of delta smelt
# the bioenergetics model was described by Rose et al. 2013a
# 1. Model temperature effect on Cmax # Rose and Eder models
# 3. Summarize CDEC data, fill-in missing values and outliers
# 4. Get food data and make prey density estimates

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
 
### 2. Summarize CDEC data
CDEC15<-read.table(file="data/CDEC15.txt",header=T) # CDEC data

#OK, X is just a matrix of all the water quality variables
X<-cbind(with(CDEC15, tapply(CDEC15$TempF_HON, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_HON, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$TempF_ANH, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_ANH, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$TempF_SSI, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_SSI, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$TempF_RVB, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_RVB, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$TempF_TWI, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_TWI, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$TempF_GZL, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_GZL, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$TempF_BDL, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_BDL, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$TempF_NSL, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_NSL, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$TempF_LIB, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_LIB, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$TempF_DWS, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$NTU_DWS, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$Cond_GZL, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$Cond_HON, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$Cond_CSE, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$Cond_SSI, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$Cond_RVB, list(CDEC15$date), FUN=mean, na.rm='T')),
with(CDEC15, tapply(CDEC15$Cond_DWS, list(CDEC15$date), FUN=mean, na.rm='T')))

X[,c(1,3,5,7,9,11,13,15,17,19)] <- (X[,c(1,3,5,7,9,11,13,15,17,19)]-32)*5/9 # convert fahrenheit to celcius
X[,(21:26)] <- gsw_SP_from_C(X[,(21:26)]/1000, X[,c(11,1,3,5,7,19)], 0) # convert conductivity to PSS
X[,c(18,20)] <- (1.1*X[,c(18,20)]-7)/0.74 # convert FNU to NTU
# Lee, C. M., Hestir, E. L., Tufillaro, N., Palmieri, B., Acuña, S., Osti, A., Bergamaschi, B. and Sommer, T. (2021). 
# Monitoring turbidity in San Francisco Estuary and Sacramento–San Joaquin delta using satellite remote sensing. 
# JAWRA Journal of the American Water Resources Association, 57(5), 737-751.

for (i in c(1,3,5,7,9,11,13,15,17,19)) { # censor outliers
 X[,i] <- ifelse(X[,i]>29,X[,i]==NA,X[,i]<-X[,i])
 X[,i] <- ifelse(X[,i]<1,X[,i]==NA,X[,i]<-X[,i])
 }
for (i in c(2,4,6,8,10,12,14,16,18,20)) {
 X[,i] <- ifelse(X[,i]>200,X[,i]==NA,X[,i]<-X[,i])
 X[,i] <- ifelse(X[,i]<0.03,X[,i]==NA,X[,i]<-X[,i])
 }
for (i in c(21:26)) {
 X[,i] <- ifelse(X[,i]>33,X[,i]==NA,X[,i]<-X[,i])
 X[,i] <- ifelse(X[,i]<0.03,X[,i]==NA,X[,i]<-X[,i])
 }

X <- X[order(as.Date(rownames(X), format="%m/%d/%Y")),]

turb.flag<-c(sd(X[,2],na.rm=T),sd(X[,4],na.rm=T),sd(X[,2],na.rm=T),sd(X[,8],na.rm=T),sd(X[,10],na.rm=T),sd(X[,12],na.rm=T),sd(X[,14],na.rm=T),sd(X[,16],na.rm=T),sd(X[,18],na.rm=T),sd(X[,20],na.rm=T))
for (h in 4:nrow(X)) {
for (i in c(2)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>2.3*turb.flag[1],NA,X[h,i])
 }
for (i in c(4)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>2.3*turb.flag[2],NA,X[h,i])
 }
for (i in c(6)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>2.3*turb.flag[3],NA,X[h,i])
 }
for (i in c(8)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>2.3*turb.flag[4],NA,X[h,i])
 }
for (i in c(10)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>2.3*turb.flag[5],NA,X[h,i])
 }
for (i in c(12)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>2.3*turb.flag[6],NA,X[h,i])
 }
for (i in c(14)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>2.3*turb.flag[7],NA,X[h,i])
 }
for (i in c(16)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>2.3*turb.flag[8],NA,X[h,i])
 }
for (i in c(18)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>1.25*turb.flag[9],NA,X[h,i])
 }
for (i in c(20)) {
 X[h,i] <- ifelse(abs(X[h,i]-mean(X[,i],na.rm=T))>2*turb.flag[10],NA,X[h,i])
 }}

glm.X <- X
jday <- yday(as.Date(dimnames(X)[[1]], format="%m/%d/%Y"))
jday2 <- jday^2

#I'm not sure what these models do
# NE Suisun
mod1<-glm(glm.X[,1]~glm.X[,c(3,4,5,6)]+jday+jday2) #,family=Gamma(link="log"))
mod2<-glm(glm.X[,2]~glm.X[,c(3,4,5)]+jday+jday2,family=Gamma(link="log"))
# Confluence
mod3<-glm(glm.X[,3]~glm.X[,c(1,2,5,6)]+jday+jday2) #,family=Gamma(link="log")) # fit models to predict missing values
mod4<-glm(glm.X[,4]~glm.X[,c(2,5,6)]+jday+jday2,family=Gamma(link="log"))
# Lower Sac
mod5<-glm(glm.X[,5]~glm.X[,c(1,3,4)]+jday+jday2) #,family=Gamma(link="log"))
mod6<-glm(glm.X[,6]~glm.X[,c(1,3,4)]+jday, family=Gamma(link="log"))
mod61<-glm(glm.X[,8]~glm.X[,c(2,18,20)]+jday, family=Gamma(link="log"))
# SJR
mod7<-glm(glm.X[,9]~glm.X[,c(1,3,5,7)]+jday) #,family=Gamma(link="log"))
mod8<-glm(glm.X[,10]~glm.X[,c(1,3,6,7,8)]+jday+jday2, family=Gamma(link="log"))
# NW Suisun
mod9<-glm(glm.X[,11]~glm.X[,c(1,13,15)]+jday+jday2) #,family=Gamma(link="log"))
mod10<-glm(glm.X[,12]~glm.X[,c(2,16)]+jday+jday2, family=Gamma(link="log"))
# Beldon
mod11<-glm(glm.X[,13]~glm.X[,c(1,3,11,15)]+jday+jday2) #,family=Gamma(link="log"))
mod12<-glm(glm.X[,14]~glm.X[,c(2,4,12,16)]+jday+jday2, family=Gamma(link="log"))
# Steel
mod13<-glm(glm.X[,15]~glm.X[,c(1,3,11,13)]+jday+jday2) #,family=Gamma(link="log"))
mod14<-glm(glm.X[,16]~glm.X[,c(2,4,12,14)]+jday2, family=Gamma(link="log"))
#Liberty
mod15<-glm(glm.X[,17]~glm.X[,c(5,7,15)]+jday+jday2) #,family=Gamma(link="log"))
mod16<-glm(glm.X[,18]~glm.X[,c(8)]+jday2, family=Gamma(link="log"))
#DWSC
mod17<-glm(glm.X[,19]~glm.X[,c(5,7,13)]+jday+jday2) #,family=Gamma(link="log"))
mod18<-glm(glm.X[,20]~glm.X[,c(8,10)]+jday+jday2, family=Gamma(link="log"))

pred.X <- matrix(NA,nrow(X),ncol(X))
new.X <- X

#OK, so we modeled all the parameters based on other parameters? I'm not sure the point of that. 
pred.X[,1] <- (mod1$coefficients[1]+mod1$coefficients[2]*X[,3]+mod1$coefficients[3]*X[,4]+mod1$coefficients[4]*X[,5]+mod1$coefficients[5]*X[,6]+mod1$coefficients[6]*jday+mod1$coefficients[7]*jday2)
pred.X[,2] <- exp(mod2$coefficients[1]+mod2$coefficients[2]*X[,3]+mod2$coefficients[3]*X[,4]+mod2$coefficients[4]*X[,5]+mod2$coefficients[5]*jday+mod2$coefficients[6]*jday2)
pred.X[,3] <- (mod3$coefficients[1]+mod3$coefficients[2]*X[,1]+mod3$coefficients[3]*X[,2]+mod3$coefficients[4]*X[,5]+mod3$coefficients[5]*X[,6]+mod3$coefficients[6]*jday+mod3$coefficients[7]*jday2)
pred.X[,4] <- exp(mod4$coefficients[1]+mod4$coefficients[2]*X[,2]+mod4$coefficients[3]*X[,5]+mod4$coefficients[4]*X[,6]+mod4$coefficients[5]*jday+mod4$coefficients[6]*jday2)
pred.X[,5] <- (mod5$coefficients[1]+mod5$coefficients[2]*X[,1]+mod5$coefficients[3]*X[,3]+mod5$coefficients[4]*X[,4]+mod5$coefficients[5]*jday+mod5$coefficients[6]*jday2)
pred.X[,6] <- exp(mod6$coefficients[1]+mod6$coefficients[2]*X[,1]+mod6$coefficients[3]*X[,3]+mod6$coefficients[4]*X[,4]+mod6$coefficients[5]*jday)
pred.X[,8] <- exp(mod61$coefficients[1]+mod61$coefficients[2]*X[,2]+mod61$coefficients[3]*X[,18]+mod61$coefficients[4]*X[,20]+mod61$coefficients[5]*jday)
pred.X[,9] <- (mod7$coefficients[1]+mod7$coefficients[2]*X[,1]+mod7$coefficients[3]*X[,3]+mod7$coefficients[4]*X[,5]+mod7$coefficients[5]*X[,7]+mod7$coefficients[6]*jday)
pred.X[,10] <- exp(mod8$coefficients[1]+mod8$coefficients[2]*X[,1]+mod8$coefficients[3]*X[,3]+mod8$coefficients[4]*X[,6]+mod8$coefficients[5]*X[,7]+mod8$coefficients[6]*X[,8]+mod8$coefficients[7]*jday+mod8$coefficients[8]*jday)
pred.X[,11] <- (mod9$coefficients[1]+mod9$coefficients[2]*X[,1]+mod9$coefficients[3]*X[,13]+mod9$coefficients[4]*X[,15]+mod9$coefficients[5]*jday+mod9$coefficients[6]*jday2)
pred.X[,12] <- exp(mod10$coefficients[1]+mod10$coefficients[2]*X[,2]+mod10$coefficients[3]*X[,16]+mod10$coefficients[4]*jday+mod10$coefficients[5]*jday2)
pred.X[,13] <- (mod11$coefficients[1]+mod11$coefficients[2]*X[,1]+mod11$coefficients[3]*X[,3]+mod11$coefficients[4]*X[,11]+mod11$coefficients[5]*X[,15]+mod11$coefficients[6]*jday+mod11$coefficients[7]*jday2)
pred.X[,14] <- exp(mod12$coefficients[1]+mod12$coefficients[2]*X[,2]+mod12$coefficients[3]*X[,4]+mod12$coefficients[4]*X[,12]+mod12$coefficients[5]*X[,16]+mod12$coefficients[6]*jday+mod12$coefficients[7]*jday2)
pred.X[,15] <- (mod13$coefficients[1]+mod13$coefficients[2]*X[,1]+mod13$coefficients[3]*X[,3]+mod13$coefficients[4]*X[,11]+mod13$coefficients[5]*X[,13]+mod13$coefficients[6]*jday+mod13$coefficients[7]*jday2)
pred.X[,16] <- exp(mod14$coefficients[1]+mod14$coefficients[2]*X[,2]+mod14$coefficients[3]*X[,4]+mod14$coefficients[4]*X[,12]+mod14$coefficients[5]*X[,14]+mod14$coefficients[6]*jday2)
pred.X[,17] <- (mod15$coefficients[1]+mod15$coefficients[2]*X[,5]+mod15$coefficients[3]*X[,7]+mod15$coefficients[4]*X[,15]+mod15$coefficients[5]*jday+mod15$coefficients[6]*jday2)
pred.X[,18] <- exp(mod16$coefficients[1]+mod16$coefficients[2]*X[,8]+mod16$coefficients[3]*jday2)
pred.X[,19] <- (mod17$coefficients[1]+mod17$coefficients[2]*X[,5]+mod17$coefficients[3]*X[,7]+mod17$coefficients[4]*X[,13]+mod17$coefficients[5]*jday+mod17$coefficients[6]*jday2)
pred.X[,20] <- exp(mod18$coefficients[1]+mod18$coefficients[2]*X[,8]+mod18$coefficients[3]*X[,10]+mod18$coefficients[4]*jday+mod18$coefficients[5]*jday2)

for (h in 1:nrow(new.X)) {
for (i in c(1:6,8:20)) {
 if(is.na(new.X[h,i])==T) { new.X[h,i] <-pred.X[h,i] } # replace missing values
 }}
for (h in 1:3) {
for (i in c(6)) {
 new.X[h,i] <- ifelse(is.na(pred.X[h,i])==F & abs(new.X[h,i]-pred.X[h,i])>turb.flag[1],pred.X[h,i],new.X[h,i])
 }}
for (h in 4:nrow(new.X)) {
for (i in c(6)) {
 new.X[h,i] <- ifelse(is.na(pred.X[h,i])==F & abs(new.X[h,i]-pred.X[h,i])>turb.flag[1] & abs(new.X[h,i]-mean(new.X[c((h-3),(h-2),(h-1)),i],na.rm=T))>turb.flag[1],pred.X[h,i],new.X[h,i])
 }} 
new.X[1956:1969,6] <- pred.X[1956:1969,6]
for (h in 1228:1245) {
for (i in c(2)) {
 new.X[h,i] <- ifelse(is.na(pred.X[h,i])==F & abs(new.X[h,i]-pred.X[h,i])>turb.flag[2] & abs(new.X[h,i]-mean(new.X[c((h-3),(h-2),(h-1)),i],na.rm=T))>turb.flag[2],pred.X[h,i],new.X[h,i])
 }}
new.X[490:493,20] <- NA
#new.X[1603,9] <- (new.X[1602,9]+new.X[1604,9])/2 # linear interpolation for some values
#new.X[1603,10] <- (new.X[1602,10]+new.X[1604,10])/2
#new.X[372,12] <- (new.X[371,12]+new.X[373,12])/2
#new.X[372,2] <- (new.X[371,2]+new.X[373,2])/2
#new.X[372,4] <- (new.X[371,4]+new.X[373,4])/2
#new.X[479,12] <- (new.X[478,4]+new.X[481,4])/2
#new.X[480,12] <- (new.X[478,4]+new.X[481,4])/2

for (i in 1:26) {
 new.X[,i] <- approx(seq(1,nrow(new.X),by=1),new.X[,i],xout=seq(1,nrow(new.X),by=1))$y
 }
low.bound<-matrix(NA,nrow(X),2)
for (i in 1:nrow(X)) {
 low.bound[i,1] <- min(which(new.X[i,(21:26)]<2))
 low.bound[i,2] <- min(which(new.X[i,(21:26)]<5.6))
 }

# graphs to check distributions of amended CDEC data
par(mfrow=c(1,2))
plot(new.X[,1],ylim=c(4,30))
points(new.X[,3],col='red')
points(new.X[,5],col='blue')

plot(new.X[,2],ylim=c(0,500))
points(new.X[,4],col='red')
points(new.X[,6],col='blue')

par(mfrow=c(3,4))
plot(new.X[1:n.days,2],type="l",ylim=c(0,100))
lines(new.X[(1+365):(n.days+365),2])
lines(new.X[(1+365+366):(n.days+365+366),2])
lines(new.X[(1+365+366+365):(n.days+365+366+365),2])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),2])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),2])
title("Grizzly Bay",line=-1)

plot(new.X[1:n.days,4],type="l",ylim=c(0,50))
lines(new.X[(1+365):(n.days+365),4])
lines(new.X[(1+365+366):(n.days+365+366),4])
lines(new.X[(1+365+366+365):(n.days+365+366+365),4])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),4])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),4])
title("Antioch",line=-1)

plot(new.X[1:n.days,6],type="l",ylim=c(0,50))
lines(new.X[(1+365):(n.days+365),6])
lines(new.X[(1+365+366):(n.days+365+366),6])
lines(new.X[(1+365+366+365):(n.days+365+366+365),6])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),6])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),6])
title("Decker Island",line=-1)

plot(new.X[1:n.days,8],type="l",ylim=c(0,50))
lines(new.X[(1+365):(n.days+365),8])
lines(new.X[(1+365+366):(n.days+365+366),8])
lines(new.X[(1+365+366+365):(n.days+365+366+365),8])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),8])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),8])
title("Rio Vista",line=-1)

plot(new.X[1:n.days,10],type="l",ylim=c(0,50))
lines(new.X[(1+365):(n.days+365),10])
lines(new.X[(1+365+366):(n.days+365+366),10])
lines(new.X[(1+365+366+365):(n.days+365+366+365),10])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),10])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),10])
title("Twitchell",line=-1)

plot(new.X[1:n.days,12],type="l",ylim=c(0,150))
lines(new.X[(1+365):(n.days+365),12])
lines(new.X[(1+365+366):(n.days+365+366),12])
lines(new.X[(1+365+366+365):(n.days+365+366+365),12])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),12])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),12])
title("NW Suisun",line=-1)

plot(new.X[1:n.days,14],type="l",ylim=c(0,50))
lines(new.X[(1+365):(n.days+365),14])
lines(new.X[(1+365+366):(n.days+365+366),14])
lines(new.X[(1+365+366+365):(n.days+365+366+365),14])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),14])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),14])
title("Beldon",line=-1)

plot(new.X[1:n.days,16],type="l",ylim=c(0,50))
lines(new.X[(1+365):(n.days+365),16])
lines(new.X[(1+365+366):(n.days+365+366),16])
lines(new.X[(1+365+366+365):(n.days+365+366+365),16])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),16])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),16])
title("Steel",line=-1)

plot(new.X[1:n.days,18],type="l",ylim=c(0,50))
lines(new.X[(1+365):(n.days+365),18])
lines(new.X[(1+365+366):(n.days+365+366),18])
lines(new.X[(1+365+366+365):(n.days+365+366+365),18])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),18])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),18])
title("Liberty",line=-1)

plot(new.X[1:n.days,20],type="l",ylim=c(0,50))
lines(new.X[(1+365):(n.days+365),20])
lines(new.X[(1+365+366):(n.days+365+366),20])
lines(new.X[(1+365+366+365):(n.days+365+366+365),20])
lines(new.X[(1+365+366+365+365):(n.days+365+366+365+365),20])
lines(new.X[(1+365+366+365+365+365):(n.days+365+366+365+365+365),20])
title("DWSC",line=-1)

new.X<-new.X[,c(11:12,1:8,19:20)] # only keep NW Suisun, NE Suisun, Conf, LowSac, SacCache, DWSC

obs.temp.dat <- abind(new.X[1:n.days,c(1,3,5,7,9,11)], # select columns where temp data are stored
 new.X[(1+366):(n.days+366),c(1,3,5,7,9,11)],
 new.X[(1+366+365):(n.days+366+365),c(1,3,5,7,9,11)],
 new.X[(1+366+365+365):(n.days+366+365+365),c(1,3,5,7,9,11)],
 new.X[(1+366+365+365+365):(n.days+366+365+365+365),c(1,3,5,7,9,11)],
 new.X[(1+366+365+365+365+366):(n.days+366+365+365+365+366),c(1,3,5,7,9,11)],along=3)
obs.turb.dat <- abind(new.X[1:n.days,c(2,4,6,8,10,12)], # select columns where turbidity data are stored
 new.X[(1+366):(n.days+366),c(2,4,6,8,10,12)],
 new.X[(1+365+366):(n.days+365+366),c(2,4,6,8,10,12)],
 new.X[(1+365+366+365):(n.days+365+366+365),c(2,4,6,8,10,12)],
 new.X[(1+365+366+365+365):(n.days+365+366+365+365),c(2,4,6,8,10,12)],
 new.X[(1+365+366+365+365+366):(n.days+365+366+365+365+366),c(2,4,6,8,10,12)],along=3)
sal.bound <- abind(low.bound[1:n.days,1:2], 
 low.bound[(1+366):(n.days+366),1:2],
 low.bound[(1+365+366):(n.days+365+366),1:2],
 low.bound[(1+365+366+365):(n.days+365+366+365),1:2],
 low.bound[(1+365+366+365+365):(n.days+365+366+365+365),1:2],
 low.bound[(1+365+366+365+365+366):(n.days+365+366+365+365+366),1:2],along=3)

make.temp<-function() {
 return(obs.temp.dat)
 }
make.turb<-function() {
 return(obs.turb.dat)
 }

### 3. Get food data and make prey density estimates
#raw.zoop.dat11 <- read.table(file='raw_zoop_dat.txt',header=T)
load("data/zoopsmwide.RData") # i made this file inthe 'ZoopSalinityGams' file in teh SFHA_synthesis project
# zoop1 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X1, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop2 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X2, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop3 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X3, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop4 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X4, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop5 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X5, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop6 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X6, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop7 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X7, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop8 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X8, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop9 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X9, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop10 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X10, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop11 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X11, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))
# zoop12 <- with(raw.zoop.dat11, tapply(raw.zoop.dat11$X12, list(raw.zoop.dat11$Year,raw.zoop.dat11$Month,raw.zoop.dat11$IBMR_strata2), FUN=mean, na.rm='T'))

zoop1 <- with(zoopsmwide, tapply(zoopsmwide$limno, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop2 <- with(zoopsmwide, tapply(zoopsmwide$othcaljuv, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop3 <- with(zoopsmwide, tapply(zoopsmwide$pdiapjuv, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop4 <- with(zoopsmwide, tapply(zoopsmwide$othclad, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop5 <- with(zoopsmwide, tapply(zoopsmwide$acartela, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop6 <- with(zoopsmwide, tapply(zoopsmwide$othclad, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop7 <- with(zoopsmwide, tapply(zoopsmwide$allcopnaup, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop8 <- with(zoopsmwide, tapply(zoopsmwide$daphnia, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop9 <- with(zoopsmwide, tapply(zoopsmwide$othcyc, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop10 <- with(zoopsmwide, tapply(zoopsmwide$other, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop11 <- with(zoopsmwide, tapply(zoopsmwide$eurytem, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))
zoop12 <- with(zoopsmwide, tapply(zoopsmwide$pdiapfor, list(zoopsmwide$Year,zoopsmwide$Month,zoopsmwide$Region), FUN=mean, na.rm='T'))



zoop <- abind(zoop1,zoop2,zoop3,zoop4,zoop5,zoop6,zoop7,zoop8,zoop9,zoop10,zoop11,zoop12,along=4)

# linear interpolation
zoop.flat <- zoop[1,,,]
for (y in 2:nrow(zoop)) {
 zoop.flat <- abind(zoop.flat,zoop[y,,,],along=1)
 }

#Monthly
#zoop.flat.filled <- zoop.flat
#for (i in 1:n.strata) {
#for (j in 1:n.prey) {
# zoop.flat.filled[,i,j] <- approx(seq(1,nrow(zoop.flat),by=1),zoop.flat[,i,j],xout=seq(1,nrow(zoop.flat),by=1))$y
# zoop.flat.filled[,i,j] <-ifelse(is.na(zoop.flat.filled[,i,j])==T,0,zoop.flat.filled[,i,j])
# }}
#yr.seq <- seq(1,nrow(zoop),by=1)
#zoop.filled <- array(NA,dim=c(nrow(zoop),12,n.strata,n.prey))
#for (y in yr.seq) {
#for (i in 1:n.strata) {
#for (j in 1:n.prey) {
# zoop.filled[y,1:12,i,j] <- zoop.flat.filled[(1+12*(y-1)):(12*y),i,j]
# }}}
 
#Daily
d15.seq <- c(15,46,74,105,135,166,196,227,258,288,319,349) # julian day of the 15th of each month in non-leap yrs
zoop.daily <- array(NA,dim=c(nrow(zoop),365,n.strata,n.prey))
for (y in yr.seq) {
for (m in 1:12) {
for (i in 1:n.strata) {
for (j in 1:n.prey) {
 zoop.daily[y,d15.seq[m],i,j] <- zoop[y,m,i,j] # Make a daily-scaled file, put monthly means on the 15th
 }}}}
 
zoop.flat.daily <- zoop.daily[1,,,] # Make y x m file into long string of daily values, for interpolation
for (y in 2:nrow(zoop)) {
 zoop.flat.daily <- abind(zoop.flat.daily,zoop.daily[y,,,],along=1)
 }
zoop.flat.filled <- zoop.flat.daily # Do interpolation
for (i in 1:n.strata) {
for (j in 1:n.prey) {
 zoop.flat.filled[,i,j] <- approx(seq(1,nrow(zoop.flat.daily),by=1),zoop.flat.daily[,i,j],xout=seq(1,nrow(zoop.flat.daily),by=1))$y
 zoop.flat.filled[,i,j] <-ifelse(is.na(zoop.flat.filled[,i,j])==T,0,zoop.flat.filled[,i,j]) # Fill in 0s at the ends
 }} 
yr.seq <- seq(1,nrow(zoop),by=1)
zoop.filled <- array(NA,dim=c(nrow(zoop),365,n.strata,n.prey)) # reshape into yr x day x strata x prey array
for (y in yr.seq) {
for (i in 1:n.strata) {
for (j in 1:n.prey) {
 zoop.filled[y,1:365,i,j] <- zoop.flat.filled[(1+365*(y-1)):(365*y),i,j]
 }}}

make.food<-function(yr) {
#Monthly
#hldr<-array(NA,dim=c(12,n.strata,n.prey)) # day x strata x prey type
#for (h in 1:12) {
#for (i in 1:n.strata) {
# hldr[h,i,] <- zoop.filled[yr,h,i,]/1000 # convert from ug/m3 to mg/cm3
# }}

hldr<-array(NA,dim=c(365,n.strata,n.prey)) # day x strata x prey type
for (h in 1:365) {
for (i in 1:n.strata) {
 hldr[h,i,] <- zoop.filled[yr,h,i,]/1000 # convert from ug/m3 to mg/cm3
 }}
 return(hldr)
 }
