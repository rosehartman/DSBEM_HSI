#Run the model for summer and fall 2011-2024
#R_DSIBM


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
source("BEMfunction.R")
### 1. Read model parameters and data
n.strata <- 7 #"Confluence"   "NE Suisun"    "NW Suisun"    "SE Suisun"    "Suisun Marsh" "SW Suisun" "Lower Sacramento RIver"
ex.strata = sort(c( "NE Suisun", "SE Suisun", "NW Suisun", "Lower Sacramento River", "Confluence", "SW Suisun",  "Suisun Marsh"))
n.prey <- 12 #may adjust this later
n.days<-153 # summer and fall
ln.a<-c(0.000005,0.00000183) # length-weight parameters
ln.b<-c(3,3.38)
startY <- 2010
endY <- 2024 #expand through 2024 when all data are in
yr.seq<- seq(startY-2010,endY-2010,by=1)
n.yrs<-15
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
T.M <- c(20,20,23,23,23) # Temperature above which things start delining. I'm bumping this down from 23 to 22
#T.M <- c(20,20,21.6,21.6,21.6) # smith and nobriga had this at 21.6
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


# The ep and es terms (J/g) were used
# to convert grams of prey per gram of Delta Smelt to grams of
# smelt per gram of smelt, which was then multiplied by weight
# (W) to yield the weight change in grams of Delta Smelt per
# individual per day. The value of es was fixed at 4,814 J/g, while ep
# was computed each day based on the fraction of Limnoithona in
# Downloaded by [University of California Davis] at 10:03 19 February 2014
# DELTA SMELT MODEL BASELINE 1245
# the diet. All zooplankton groups had an energy density of 2,590
# J/g; the exception was Limnoithona, for which energy density
# was assumed to be 30% lower (1,823 J/g) 

source("BEM consumption functions_v3.2.R") # load consumption functions
load("WaterQuality20102024.RData")

obs.temp.dat <- Tempx2
obs.turb.dat <- Turbx2
daylight = daylength(38.15, c(130:310))$Daylength*60
LT.fx <- daylight/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace

#day <- CDEC[,2]
#LT.fx <- daylength/max(daylength)
wtr.yr = c("BN", "W", "BN", "D", "C", "C", "BN", "W", "BN", "W", "D", "C", "C", "W", "AN") #water year types 2010-2024
#wtr.yr<-c(1,1,1,1,1,1,4,4,2,3,3,1,4,5,4,3,1,3,4,5,5) # Sacto WY type wet = 1,above normal = 2, below normal = 3, dry = 4, critical = 5
#rep.yrs <- seq(21,26,by=1) # years 2015-2020 from 1995-2020 dataset, for prey density summary
#ex.strata <- c(13,10,8,6,3,1) # map IBMR strata to 5 CDEC regions




beta_hat <- read.table(file='data/beta_hat.txt') # MC filtering coefficients
#this is model uncertainty
#beta_hat <- read.table(file='BEM consumption rate project/beta_hat_IND.txt') # MC filtering coefficients


### 3. Model bioenergetics and reference points
load("zoopsmwide.RData")
#OK, this just fills in prey densities. I can do this with the zoop data. 
PD.mn.array = zoopx2 #days by prey by strata

sfharun = smelt_bioenergetics(PD.mn.array = zoopx2, obs.temp.dat = Tempx2, obs.turb.dat = Turbx2, 
                    start.L = rep(23, 15), ex.strata = ex.strata, beta_hat = beta_hat[1:200,], start.year = 2010)



summary(sfharun)

ggplot(sfharun, aes(x = Day, y = Weight, color = Stratum))+
  facet_wrap(~Year)+
  geom_smooth()
#something is screwy here. Missing some data and probably have it assigned to the wrong years. SIgh. 

ggplot(sfharun, aes(x = Day, y = Length, color = as.factor(Year)))+
  facet_wrap(~Stratum)+
  geom_smooth()

load("data/sal.Rdata")
load("data/DF.RData")

Wt.summary = sfharun %>%
  rename(Region = Stratum) %>%
  group_by(Region, Year, Day) %>%
  summarize(MWeight = mean(Weight, na.rm =T), sdWeight = sd(Weight, na.rm =T),
            MLength = mean(Length, na.rm =T), sdLength = sd(Length, na.rm =T)) 


Wt.summary = left_join(Wt.summary, sal) 

Wt.sumx = left_join(Wt.summary, DF)

 ggplot(Wt.sumx, aes(x = Day, y = MWeight, color = X2, group = as.factor(Year))) +
   geom_smooth()+
   facet_wrap(~Region2)+
   ylab("Predicted Weight in Grams")
 
 
 ggplot(Wt.summary, aes(x = Day, y = MWeight, color = as.factor(Year))) +
   geom_smooth()+
   facet_wrap(~Region2)+
   ylab("Predicted Weight in Grams")
 

 ggplot(filter(Wt.summary, Year == 2014), aes(x = Day, y = MWeight, color = Region2)) +
   geom_smooth()+ ylab("Weight (g)")+ xlab("Day of Year")+ theme_bw()+
   theme(legend.position = "bottom")
 

 #I need to check and make sure this is what i'd expect based on food, temp, turbidity
 

 #OK, how much did they grow over the whole summer (or fall) and is it related to X2?
 Wtsum2 = group_by(Wt.summary, Year, Region) %>%
   summarize(startweight = first(MWeight), endwieght = last(MWeight), diffweight = endwieght-startweight)
 
 
ggplot(Wtsum2, aes(x = Region, y = diffweight, fill = Region)) + geom_col(position = "dodge")+
  facet_wrap(~Year, scales = "free_x")+
  theme(axis.text.x = element_blank())
 

DFmeans = filter(DF, Month %in% c(6:9)) %>%
  group_by(Year) %>%
  summarize(OUT = mean(OUT, na.rm =T), X2 = mean(X2, na.rm =T))

Wtsum2 = left_join(Wtsum2, DFmeans)

ggplot(Wtsum2, aes(x = X2, y = diffweight)) + geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Region, scales = "free_x")+

  ylab("Predicted growth (g), June-OCtober")+
  xlab("Mean June-September X2")

x2grwothlm = lm(diffweight ~ X2*Region, data = Wtsum2)
summary(x2grwothlm)
library(emmeans)
trends = emtrends(x2grwothlm, pairwise ~ Region, var = "X2")
trends
#OK, so maybe I don't worry about the weird things in DW suisun that much because it isn't significant


save(Wt.summary, sfharun, DF, Wtsum2, file = "outputs/summerflowweights.RData")
############################################################
#OK, now I need to throw out growht when it gets too salty

#Maybe use growth per day, distribute smelt based on where they can grow, and figure out an average?


Wt.summary = Wt.summary %>%
  group_by(Region, Year) %>%
  mutate(WperD = MWeight - lag(MWeight),
         WperD = replace_na(WperD, 0)) %>%
  ungroup()

#I should use the full dataset and get the varience. Somehow. 
#also weight gain should be g/g.
salx = select(sal, Region2, Day, Year, GoodSalinity) %>%
  distinct()
Wtsal = sfharun %>%
  rename(Region = Stratum) %>%
  group_by(s, Region, Year) %>%
  arrange(Day) %>%
  mutate(gpergperd = (Weight -lag(Weight))/lag(Weight),
         gpergperd = replace_na(gpergperd, 0)) %>%
  ungroup() %>%
  left_join(salx)




Meangrowth6a = Wtsal %>%
  group_by(Year, Day) %>%
  filter(GoodSalinity) %>%
  summarise(Growth = mean(gpergperd, na.rm =T), Weight = mean(Weight, na.rm =T)) %>%
  group_by(Year) %>%
  mutate(Weight = cumsum(Growth*Weight)+0.073298) %>%
  ungroup()



ggplot(Meangrowth6a, aes(x = Day, y = Weight,color = as.factor(Year))) + geom_line()

ggplot(Meangrowth6a, aes(x = Day, y = Growth,color = as.factor(Year))) + geom_line()
#why did my results change so much with th enew data?


#Does this make sense? need to check with the group. 

Meangrowthsum = Meangrowth6a %>%
  group_by(Year) %>%
  summarise(Weight = max(Weight, na.rm =T)) %>%
  left_join( DFmeans) %>%
  left_join(Meantemps6)

#x2 versus growth 6 psu
ggplot(Meangrowthsum, aes(x = X2, y = Weight)) + geom_point()+ geom_smooth(method = "lm")+
  geom_text(aes(label = Year))

x2lm6 = lm(Weight ~ X2, Meangrowthsum)
summary(x2lm6)
#barely significant

ggplot(Meangrowthsum, aes(x = X2, y = Weight)) + geom_point()+ geom_smooth()+
  geom_text(aes(label = Year))



ggplot(Meangrowthsum, aes(x = Temperature, y = Weight)) + geom_point()+ geom_smooth(method = "lm")+
  geom_text(aes(label = Year))


#x2 versus growth 2 psu
#ggplot(Meangrowthsum, aes(x = X2, y = Weight2)) + geom_point()+ geom_smooth(method = "lm")+
#  geom_text(aes(label = Year))

#x2lm2 = lm(Weight2 ~ X2, Meangrowthsum)

#summary(x2lm2)
#much more significant

#ggplot(Meangrowthsum, aes(x = Temperature, y = Weight2)) + geom_point()+ geom_smooth(method = "lm")+
#  geom_text(aes(label = Year))


########################
#x2 versus good temperature

TempsandX2 = left_join(Meantemps, DF)

ggplot(TempsandX2, aes(x = X2, y = Temperature)) + geom_point()+ geom_smooth(method = "lm")

#temperature versus X2 by region
tempregion = AllWQmean2 %>%
  filter(Month %in% c(6:9), Parameter == "watertemperature") %>%
  group_by(Region, Year, DOY, Date)  %>%
  summarize(Temp = mean(Value, na.rm =T)) %>%
  left_join(mutate(Dayflow, Date = mdy(Date)))


ggplot(tempregion, aes(x = Region, y = Temp)) + geom_boxplot()

tempregionmean = group_by(tempregion, Region, Year) %>%
  summarize(Temp = mean(Temp), X2 = mean(X2))

ggplot(tempregionmean, aes(x = Year, y = Temp, color = Region)) + geom_point()+ geom_line()+
  ylab("average Summer temperature")+
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020, 2022))

ggplot(tempregionmean, aes(x = X2, y = Temp, color = Region)) + geom_point() + geom_smooth(method = "lm", se = F)

library(deltamapr)
ggplot()+
  geom_sf(data = WW_Delta)+
  geom_sf(data = Regions, aes(fill = Region), alpha = 0.5)+
  coord_sf(xlim = c(-122.2, -121.6), ylim = c(38.0, 38.3))

#zoops versus X2 by region and year and stuff

zoopsave = zoopsallm %>%
  filter(Month %in% c(6:9)) %>%
  group_by(Year, IBMR, Region) %>%
  summarize(BPUE = mean(BPUE, na.rm =T))

ggplot(zoopsave, aes(x = Year, y = BPUE, fill = IBMR)) +
  facet_wrap(~Region)+ geom_col()

##########################################################
