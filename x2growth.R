#OK, how much change in X2 is required for a change in smelt growth?

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
library(arrayhelpers)
source("BEMfunction.R")
library(RColorBrewer)


mypal = c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))


#let's start with average X2 by month for the whole time period
#it might be better to do just the wetter years, but IDK
load("data/DF.RData")
# X2means = DF %>%
#   filter(Year>2009) %>%
#   group_by(Month) %>%
#   summarise(X2 = mean(X2, na.rm =T)) 
# 
# write.csv(X2means, file = "data/X2means.csv")

X2scenarios = read_csv("data/X2means.csv") %>%
  pivot_longer(cols = c(X2:last_col()), names_to = "Scenario", values_to = "X2") %>%
  filter(Month %in% c(6:10))

#load X2 model
load("C:/Users/rhartman/OneDrive - California Department of Water Resources/salinity control gates/SFHA_synthesis/outputs/zoopx2model.RData")
summary(zm3)

#so i just need x2, region, season, and month. one for each day of the period


alldaysSF = seq(ymd("2024-06-01"), ymd("2024-10-31"), by = "day")
RegsSF = c("River", "Suisun Marsh", "Suisun Bay", "Grizzly Bay")

allvalsSF = data.frame(Date = alldaysSF) %>%
  mutate(Month = month(Date), Season = case_when(Month %in% c(6,7,8) ~ "Summer",
                                                 Month %in% c(9,10) ~ "Fall"),
         Year = year(Date)) %>%
  merge(data.frame(Region = RegsSF)) %>%
  full_join(X2scenarios)

Pesudopredict = predict(zm3, newdata = allvalsSF)

Predictallvalssf = bind_cols(allvalsSF, data.frame(prediction = Pesudopredict))

#OK, now I have to convert that to biomass. The model is adults and copepodites, so its a bit tricky. 
#probably make an assumption about hte ratio of adults to juveniles based on past data.

ggplot(Predictallvalssf, aes(x = Date, y = prediction, color = Region))+ geom_line()+
  facet_wrap(~Scenario)

ggplot(filter(Predictallvalssf, Scenario %in% c("X2", "X2minus2Summer", "X2minus2Fall")), aes(x = Date, y = prediction, color = Scenario))+ geom_line()+
  facet_wrap(~Region)

#this is the average critters as biomass in mg
#pseudo adult has a carbon weight of 2.66 ug, juv is 1.15
load("data/zoopsmwide_constant.RData")
ratioAtoJ = zoops_constant %>%
  ungroup() %>%
  select(Region, Month, doy, pdiapjuv, pdiapfor) %>%
  distinct() %>%
  mutate(Ratio = (pdiapfor/2.66)/(pdiapjuv/1.15+ pdiapfor/2.66))
#so adults are 40-50% of the total count. Nice. 

Predictallvalssf = mutate(Predictallvalssf,pdiapfor = exp(prediction)*0.45*2.66/1000, pdiapjuv = exp(prediction)*0.55*1.15/1000)

ggplot(X2scenarios, aes(x = Month, y = X2, color = Scenario)) + geom_line()

#alright, now attach those values to the rest of the zoops.
#Need ot reorganize the regions to match tho.
load("data/zoopsmwide.RData")

zoopsmwidef2 = filter(zoopsmAve, doy %in% c(153:305))%>%
  arrange(Region) %>%
  mutate(Region = case_when(Region %in% c("Confluence", "Lower Sacramento River") ~ "River",
                            Region %in% c("NE Suisun" , "SE Suisun", "SW Suisun") ~ "Suisun Bay",
                            Region == "Suisun Marsh" ~ "Suisun Marsh", 
                            Region == "NW Suisun" ~ "Grizzly Bay")) %>%
  group_by(Region, Date, Year, doy, Month, Day) %>%
  summarize(across(limno:pdiapfor, mean))

zoopsX2s = zoopsmwidef2 %>%
  filter(Year == 2024) %>%
  select(-pdiapfor, -pdiapjuv)%>%
  left_join(select(Predictallvalssf, Date, Month, Year, Scenario, pdiapfor, pdiapjuv, Region))

#why is this wrong?
ggplot(zoopsX2s, aes(x = doy, y = pdiapfor, color = Scenario))+
  geom_line()+ facet_wrap(~Region)+
  geom_line(aes(y = othcaljuv))


zoopsX2ready = zoopsX2s %>%
  ungroup() %>%
  arrange(Region) %>%
  select(Scenario, Region,doy, 
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) 

test_ps_x2= zoopsX2ready%>%
  split(list(zoopsX2ready$Region, zoopsX2ready$Scenario))

foo = unlist(test_ps_x2)

zoop_ps_x2= array(unlist(test_ps_x2),dim=c(153,15,4,11), 
               dimnames = list(c(153:305), names(zoopsX2ready),
                               unique(zoopsX2ready$Region), sort(unique(zoopsX2ready$Scenario))))




#OK! Now i just need to get rid of the id columns
zoopx_psX2 = zoop_ps_x2[(1:153), c(4:15), c(1:4),c(1:11)]
zoopx2_psX2 = apply(zoopx_psX2, c(2,3,4), as.numeric) #days by taxa by scenario by region
#days by prey by strata by year!

#try to put it back together 
test = array2df(zoopx2_psX2,levels = list(NA, names(zoopsX2ready)[4:15],
                                         unique(zoopsX2ready$Region),unique(zoopsX2ready$Scenario)),  
                label.x = "Weight", na.rm =T) %>%
  rename(Day = d1, IBMR = d2, Region = d3, Scenario = d4)

ggplot(filter(test, IBMR == "pdiapfor"), aes(x = Day, y = Weight, color = Scenario))+
  facet_wrap(~Region)+ geom_line()
#OK!.

#constant temperature, just redo the regions ##################################
load("WaterQuality20102024.RData")
Tempave_x2s = filter(AllWQmean2, Parameter == "watertemperature", DOY %in% c(153:305)) %>%
  mutate(Region = case_when(Region %in% c("Confluence", "Lower Sacramento River") ~ "River",
                            Region %in% c("NE Suisun" , "SE Suisun", "SW Suisun") ~ "Suisun Bay",
                            Region == "Suisun Marsh" ~ "Suisun Marsh", 
                            Region == "NW Suisun" ~ "Grizzly Bay")) %>%
  group_by(Region, DOY) %>%
  summarize(Value = mean(Value, na.rm =T))
  
Scenarios = data.frame(Scenario = unique(Predictallvalssf$Scenario))
  
Tempwide_constantX2 = Tempave_x2s %>%
  mutate(Year = 2024) %>%
  merge(Scenarios) %>%
  arrange(Region) %>%
  pivot_wider(names_from = Region, values_from = Value)


Tempwidef_constantX2 = filter(Tempwide_constantX2, DOY %in% c(153:305))

test_constantX2 = Tempwidef_constantX2  %>%
  split(list(Tempwidef_constantX2$Scenario))

Temp_constantX2 = array(unlist(test_constantX2),dim=c(153,7, 11), 
                      dimnames = list(c(153:305), names(Tempwidef_constantX2),
                                      sort(unique(Tempwidef_constantX2$Scenario))))

#OK! Now i just need to get rid of the id columns
Tempx_constantX2 = Temp_constantX2[c(1:153), c(4:7), c(1:11)]
Tempx2_constantX2 = apply(Tempx_constantX2, c(2,3), as.numeric)

#######now turbidity ###############################################

Turbave_x2s = filter(AllWQmean2, Parameter == "turbidity", DOY %in% c(153:305)) %>%
  mutate(Region = case_when(Region %in% c("Confluence", "Lower Sacramento River") ~ "River",
                            Region %in% c("NE Suisun" , "SE Suisun", "SW Suisun") ~ "Suisun Bay",
                            Region == "Suisun Marsh" ~ "Suisun Marsh", 
                            Region == "NW Suisun" ~ "Grizzly Bay")) %>%
  group_by(Region, DOY) %>%
  summarize(Value = mean(Value, na.rm =T))

Scenarios = data.frame(Scenario = unique(Predictallvalssf$Scenario))

Turbwide_constantX2 = Turbave_x2s %>%
  mutate(Year = 2024) %>%
  merge(Scenarios) %>%
  arrange(Region) %>%
  pivot_wider(names_from = Region, values_from = Value)


Turbwidef_constantX2 = filter(Turbwide_constantX2, DOY %in% c(153:305))

test_constantX2 = Turbwidef_constantX2  %>%
  split(list(Turbwidef_constantX2$Scenario))

Turb_constantX2 = array(unlist(test_constantX2),dim=c(153,7, 11), 
                        dimnames = list(c(153:305), names(Turbwidef_constantX2),
                                        sort(unique(Turbwidef_constantX2$Scenario))))

#OK! Now i just need to get rid of the id columns
Turbx_constantX2 = Turb_constantX2[c(1:153), c(4:7), c(1:11)]
Turbx2_constantX2 = apply(Turbx_constantX2, c(2,3), as.numeric)

##################################################################
#run th emodel################################################
n.strata <- 4 #"Confluence"   "NE Suisun"    "NW Suisun"    "SE Suisun"    "Suisun Marsh" "SW Suisun" "Lower Sacramento RIver"
ex.strata = sort(c( "Grizzly Bay", "River", "Suisun Bay",  "Suisun Marsh"))
n.prey <- 12 #may adjust this later
ln.a<-c(0.000005,0.00000183) # length-weight parameters
ln.b<-c(3,3.38)


#I'm not going to change any of these
VB.k<-2.72 #2.98 # VBGM parameters estimated in DSM TN 35
Linf<-78.39 #76.1
t0<--0.026 #-0.014
stage =4


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


zooprunX2 = smelt_bioenergetics(PD.mn.array = zoopx2_psX2, obs.temp.dat = Tempx2_constantX2, obs.turb.dat = Turbx2_constantX2, 
                              start.L = rep(23, 15), ex.strata =c("Grizzly Bay", "River", "Suisun Bay", "Suisun Marsh"),
                              beta_hat = beta_hat[1:200,], start.year = 2010)


#relable the zooplankton levels correctly
unique(zooprunX2$Year)

zooprunX2a = mutate(zooprunX2, Scenario = case_match(Year, 2010 ~ "X2",2011 ~ "X274fall",
                                                     2012 ~ "X274sum80fall",
                                                     2013 ~ "X274summer",
                                                     2014 ~ "X280fall",
                                                      2015~ "X2minus2all", 2016 ~ "X2minus2Fall", 
                                             2017 ~ "X2minus2Summer", 2018 ~ "X2minus4all",
                                             2019 ~ "X2minus4Fall", 2020 ~  "X2minus4Summer"),
                 Month = case_match(Day, c(0:30) ~ 6,
                                    c(31:62) ~7,
                                    c(63:94) ~8,
                                    c(95:125) ~9,
                                    c(126:153) ~10)) %>%
  rename(Region = Stratum)


zooprunsumX2 = zooprunX2a %>%
  group_by(Day, Month, Scenario, Region) %>%
  summarize(MWeight = mean(Weight))

ggplot(zooprunsumX2, aes(x = Day, y = MWeight, color = Scenario)) +
  geom_smooth()+
  facet_wrap(~Region)+
  scale_color_manual(values = mypal)

#Hm. Not shrinking, but also not much difference in scenarios. 

ggplot(filter(zooprunsumX2, Scenario %in% c("X2", "X274summer", "X274fall", "X280fall", "X274sum80fall")),
              aes(x = Day, y = MWeight, color = Scenario)) +
  geom_smooth()+
  facet_wrap(~Region)+
  scale_color_manual(values = mypal)

ggplot(filter(zooprunsumX2, Scenario %in% c("X2", "X274summer", "X274fall", "X280fall", "X274sum80fall"), Region == "Suisun Bay"),
       aes(x = Day, y = MWeight, color = Scenario)) +
  geom_smooth()+
  scale_color_manual(values = mypal)+
  coord_cartesian(xlim = c(120, 155), ylim = c(1.75, 2.25))


ggplot(filter(zooprunsumX2, Scenario %in% c("X2", "X274summer", "X274fall", "X280fall", "X274sum80fall"), Region == "Suisun Marsh"),
       aes(x = Day, y = MWeight, color = Scenario)) +
  geom_smooth()+
  scale_color_manual(values = mypal)+
  coord_cartesian(xlim = c(120, 155), ylim = c(1.75, 2.5))

#now restrict it to LSZ? ###############################
#UGH! How do i do that with the different X2 scenarios? Maybe base it on X2? Grizzly is only available if X2 is less than 80?

#I"m thinking I don't do this, Or maybe just show Suisun Bay results?

#maybe it makes more sense to just caculate the Pseudo percent of average in each scenario

load("data/pseudosave.RData")

percentpseudo = zoopsX2ready %>%
  filter(Region == "Suisun Bay") %>%
  mutate(Month = case_when(doy %in% c(152:181) ~ 6,
                           doy %in% c(182:212) ~ 7,
                           doy %in% c(213:243) ~8,
                           doy %in% c(244:273) ~ 9,
                           doy > 273 ~ 10)) %>%
  left_join(rename(filter(pseudosave_monthly, level ==1), pdiapAve = pdiapfor)) %>%
  group_by(Scenario) %>%
  summarise(pdiapfor = mean(pdiapfor), pdiapAve = mean(pdiapAve), percent = pdiapfor/pdiapAve)

ggplot(percentpseudo, aes(x = Scenario, y = percent))+ geom_col()
