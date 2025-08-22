#seperate summer and fall


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
source("BEMfunction.R")
### 1. Read model parameters and data
n.strata <- 7 #"Confluence"   "NE Suisun"    "NW Suisun"    "SE Suisun"    "Suisun Marsh" "SW Suisun" "Lower Sacramento RIver"
ex.strata = sort(c( "NE Suisun", "SE Suisun", "NW Suisun", "Lower Sacramento River", "Confluence", "SW Suisun",  "Suisun Marsh"))
n.prey <- 12 #may adjust this later
n.days<-92 # just one season
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

#source("BEM consumption functions_v3.2.R") # load consumption functions
load("data/WaterQuality20102024.RData")
### 3. Model bioenergetics and reference points
load("data/zoopsmwide.RData")
Temp_summer = Tempx2[1:62,,]
Temp_fall = Tempx2[92:153,,]

Turb_summer = Turbx2[1:62,,]
Turb_fall = Turbx2[92:153,,]

zoopX2summer = zoopx2[1:62,,,]
zoopX2fall = zoopx2[92:153,,,]


beta_hat <- read.table(file='data/beta_hat.txt') # MC filtering coefficients
#this is model uncertainty
#beta_hat <- read.table(file='BEM consumption rate project/beta_hat_IND.txt') # MC filtering coefficients



#OK, this just fills in prey densities. I can do this with the zoop data. 

sfharunsummer = smelt_bioenergetics(PD.mn.array = zoopX2summer, obs.temp.dat = Temp_summer, 
                                    obs.turb.dat = Turb_summer, 
                              start.L = rep(23, 15), ex.strata = ex.strata, 
                              beta_hat = beta_hat[1:200,], start.year = 2010, Season = "Summer")
sfharunsummer2 = mutate(sfharunsummer, Season = "Summer")

sfharunfall = smelt_bioenergetics(PD.mn.array = zoopX2fall, obs.temp.dat = Temp_fall, 
                                    obs.turb.dat = Turb_fall, 
                                    start.L = rep(23, 15), ex.strata = ex.strata, 
                                    beta_hat = beta_hat[1:200,], start.year = 2010, Season = "Fall")
sfharunfall2 = mutate(sfharunfall, Day = Day +90, Season = "Fall")


load("data/sal.Rdata")

load("data/DF.RData")

Wt.summary1 = bind_rows(sfharunfall2, sfharunsummer2) %>%
  rename(Region = Stratum) %>%
  group_by(Region, Year, Day, Season) %>%
  summarize(MWeight = mean(Weight, na.rm =T), sdWeight = sd(Weight, na.rm =T),
            MLength = mean(Length, na.rm =T), sdLength = sd(Length, na.rm =T)) 


Wt.summary1 = left_join(Wt.summary1, sal) 

Wt.sumx1 = left_join(Wt.summary1, DF)

ggplot(Wt.sumx1, aes(x = Day, y = MWeight, color = X2, group = as.factor(Year))) +
  geom_smooth()+
  facet_grid(Season~Region2)+
  ylab("Predicted Weight in Grams")



#OK, how much did they grow over the whole summer (or fall) and is it related to X2?
Wtsum2seasonal = group_by(Wt.summary1, Year, Region, Season) %>%
  summarize(startweight = first(MWeight), endwieght = last(MWeight), diffweight = endwieght-startweight)


DFmeans2 = filter(DF, Month %in% c(6:10)) %>%
  mutate(Season = case_when(Month %in% c(6:7) ~ "summer",
                            Month %in% c(9,10) ~ "fall")) %>%
  group_by(Year, Season) %>%
  summarize(OUT = mean(OUT, na.rm =T), X2 = mean(X2, na.rm =T))

Wtsum2seasonal = left_join(Wtsum2seasonal, DFmeans)

ggplot(Wtsum2seasonal, aes(x = X2, y = diffweight, color = Season)) + geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Region, scales = "free_x")+
  
  ylab("Predicted growth (g), seasonal")+
  xlab("Mean seasonal X2")

#now growth at les sthan 6 PSU


############################################################
#OK, now I need to throw out growht when it gets too salty

#Maybe use growth per day, distribute smelt based on where they can grow, and figure out an average?

#I should use the full dataset and get the varience. Somehow. 
#also weight gain should be g/g.
salx = ungroup(sal) %>%
  select(Region2, Day, Year, GoodSalinity) %>%
  distinct()

Wtsalseasonal = bind_rows(sfharunfall2, sfharunsummer2) %>%
  rename(Region = Stratum) %>%
  group_by(s, Region, Year, Season) %>%
  arrange(Day) %>%
  mutate(gpergperd = (Weight -lag(Weight))/lag(Weight),
         gpergperd = replace_na(gpergperd, 0)) %>%
  ungroup() %>%
  left_join(salx)




Meangrowth6asa = Wtsalseasonal %>%
  group_by(Year, Day, s, Season) %>%
  filter(GoodSalinity) %>%
  summarise(Growth = mean(gpergperd, na.rm =T), Weight = mean(Weight, na.rm =T)) %>%
  group_by(Year, s, Season) %>%
  mutate(Weight = cumsum(Growth*Weight)+0.073298) %>%
  ungroup()

Meangrowth6bsa = group_by(Meangrowth6asa, Year, Day, Season) %>%
  summarise(Growth = mean(Growth), Weight = mean(Weight)) %>%
  mutate(Seasonx = factor(Season, levels = c("Summer", "Fall"),
                          labels = c("Summer (Jun-Jul)", "Fall (Sep-Oct)")))

ggplot(Meangrowth6bsa, aes(x = Day, y = Weight,color = as.factor(Year))) + 
  geom_line(linewidth = 1)+
  facet_wrap(~Seasonx, scales = "free_x")+
  scale_x_continuous(breaks = c(1, 31, 62, 92, 123, 153), labels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Noc"))+
  scale_color_viridis_d(option = "turbo", name = NULL)+
  ylab("Predicted weight (g) in regions <6 PSU")+
  theme_bw()


ggsave("plots/sfha_growthlessthan6_seasonal.png", device = "png", width =8, height = 4)

Meangrowthseas = Meangrowth6asa %>%
  group_by(Year, s, Season) %>%
  summarise(Weight = max(Weight, na.rm =T)) %>%
  group_by(Year, Season) %>%
  summarize(sdWeight = sd(Weight), Weight = mean(Weight)) %>%
  left_join( mutate(DFmeans2, Season = case_match(Season, "summer" ~ "Summer",
                                                  "fall" ~ "Fall")))%>%
  left_join(mutate(Meantemps6b, Season = case_match(Season, "summer" ~ "Summer",
                                                   "fall" ~ "Fall"))) %>%
  mutate(Season = factor(Season, levels = c("Summer", "Fall"),
                         labels = c("Summer (Jun-Jul)", "Fall (Sep-Oct)"))) 

ggplot(Meangrowthseas, aes(x = X2, y = Weight)) + geom_point()+ 
  geom_text(aes(label = Year))+
  geom_smooth(method = "lm")+
  facet_wrap(~Season, scales = "free")+
  theme_bw()+
  ylab("Predicted Growth (g)")+
  xlab("Mean X2")

ggsave("plots/GrowthSeasonal.tiff", device = 'tiff', width =8, height = 4.5)

#linear model of final weight versus X2, with interaction of season

seaslm = lm(Weight ~Season*X2, data = Meangrowthseas)
summary(seaslm)
write.csv(summary(seaslm)$coefficients, file = "outputs/seasonalgrowth_lm.csv")

library(emmeans)
emtrends(seaslm, pairwise ~ Season, var = "X2")
#now growth versus temperature
ggplot(Meangrowthseas, aes(x = Temperature, y = Weight)) + geom_point()+ 
  geom_text(aes(label = Year))+
  geom_smooth(method = "lm")+
  facet_wrap(~Season, scales = "free")+
  theme_bw()+
  ylab("Predicted Growth (g)")+
  xlab("Mean Temperature (c)")

ggsave("plots/GrowthSeasonal_temps.tiff", device = 'tiff', width =8, height = 4.5
)
