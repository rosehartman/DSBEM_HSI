#How much do you need to increase Pseudodiaptomus in order to see a difference in smelt?

#In this case I think I"ll use average temp and turbidity, use different zoop scenarios instead of years


#Version with constant Temperature and turbididy
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
library(arrayhelpers)

library(patchwork)
library(RColorBrewer)
source("BEMfunction.R")

mypal = c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))

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
n.yrs<-11
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



####Sensativity analysis ##########################
#See if I can change psdudos in each year, but hold other taxa the same. 

load("data/zoopsmwide.RData")

zoopsmwidef = filter(zoopsmAve, doy %in% c(153:305))%>%
  arrange(Region)

#average biomass of pseudodiaptomus by month in Suisun Bay
# I should probably do this in suisun marsh instead
zoopsSuisun = filter(zoopsmwidef, Region %in% c("NE Suisun", "NW Suisun","SE Suisun")) %>%
  group_by(doy) %>%
  summarize(across(limno:pdiapfor, function(x) mean(x, na.rm = TRUE))) %>%
  mutate(level =1)
zoopsS = zoopsSuisun

#look at the mean and 95% conf interval by day 

zoopsSuisunSD = filter(zoopsmwidef, Region %in% c("NE Suisun", "NW Suisun","SE Suisun")) %>%
  group_by(doy) %>%
  summarize(pdiapforMean = mean(pdiapfor, na.rm =T), pdiapjuvMean = mean(pdiapjuv, na.rm=T),
            MaxPF = max(pdiapfor, na.rm = T), MaxPJ = max(pdiapjuv, na.rm = T),
            MinPF = min(pdiapfor, na.rm = T), MinPJ = min(pdiapjuv, na.rm = T),
            sdPF = sd(pdiapfor, na.rm =T), sdpj = sd(pdiapjuv, na.rm = T)) %>%
  mutate(Totmean = pdiapforMean + pdiapjuvMean, TotMax = MaxPF+MaxPJ,
         Totmin = MinPF + MinPJ, sdtot = (sdPF+sdpj)/2)

ggplot(zoopsSuisunSD, aes(x =doy, y = Totmean)) + geom_point()+
  geom_errorbar(aes(ymin = Totmin, ymax = TotMax))

ggplot(zoopsSuisunSD, aes(x =doy, y = Totmean)) + geom_point()+
  geom_errorbar(aes(ymin = Totmean-sdtot, ymax = Totmean+sdtot))


ggplot(zoopsSuisunSD) + geom_smooth(aes(x =doy, y = Totmean), se = F)+
  geom_smooth(aes(x =doy, y = Totmean+sdtot), se = F, linetype = 2)+
  geom_smooth(aes(x =doy, y = Totmean-sdtot), se = F, linetype = 2)+
  ylab("Pseudodiaptomus biomass (mg/m3)")

annualpseudo = filter(zoopsmwidef, Region %in% c("NE Suisun", "NW Suisun","SE Suisun")) %>%
  mutate(totP = pdiapjuv+pdiapfor, Month = month(Date)) %>%
  group_by(Year, Month) %>% 
  filter(Month != 11)%>%
summarise(totP = mean(totP, na.rm =T))

ggplot(annualpseudo , aes(x = Year, y = totP, fill = as.factor(Month)))+ 
  geom_col(position = "dodge")+
  ylab("Pseudodiaptomus biomass (mg/m3)")

ggplot(annualpseudo , aes(x = Year, y = totP, fill = as.factor(Month)))+ 
  geom_col(position = "dodge")+ facet_wrap(~Month, scales = "free_y")+
  ylab("Pseudodiaptomus biomass (mg/m3)")

#now replicate this for each biomass level
Levels = c(0.01, 0.5, 1.5, 2, 5, 10, 20, 40, 100, 200,500)
for(i in 1:10){
  Level = Levels[i]
  Pseudo = zoopsSuisun %>%
    mutate(level = Level, pdiapfor = level*pdiapfor, pdiapjuv = level*pdiapjuv)
  zoopsS = bind_rows(zoopsS, Pseudo)
}

pseudosave = ungroup(zoopsS) %>%
  group_by(level) %>%
  summarize(pdiapfor = exp(mean(log(pdiapfor), na.rm =T)),
            pdiapjuv = exp(mean(log(pdiapjuv), na.rm =T))) %>%
            mutate(pdiaptot = pdiapfor + pdiapjuv)
yr.seq = pseudosave$level

pseudosave_monthly = ungroup(zoopsS) %>%
  mutate(Month = case_match(doy, c(153:183) ~ 6,
                            c(184:215) ~7,
                            c(216:247) ~8,
                            c(248:278) ~9,
                            c(279:305) ~10)) %>%
  group_by(level, Month) %>%
  summarize(pdiapfor = mean(pdiapfor, na.rm =T))

save(zoopsSuisun, pseudosave, pseudosave_monthly, file = "data/pseudosave.RData")

#now the different temperature regimes

zoopsSS = bind_rows(mutate(zoopsS, Temp = "Temp_05"),
                    mutate(zoopsS, Temp = "Tempave"),
                    mutate(zoopsS, Temp = "Temp1"),
                    mutate(zoopsS, Temp = "Temp15"),
                    mutate(zoopsS, Temp = "Temp05"),
                    mutate(zoopsS, Temp = "Temp_1"),
                    mutate(zoopsS, Temp = "Temp_15"))  %>%
  select(Temp, level,doy, 
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) 

test_ps= zoopsSS%>%
  split(list(zoopsSS$Temp, zoopsSS$level))

zoop_ps= array(unlist(test_ps),dim=c(153,15,7, 11), 
               dimnames = list(c(153:305), names(zoopsSS), unique(zoopsSS$Temp),
                               unique(zoopsSS$level)))


#OK! Now i just need to get rid of the id columns
zoopx_ps = zoop_ps[(1:153), c(4:15), c(1:7), c(1:11)]
zoopx2_ps = apply(zoopx_ps, c(2,3,4), as.numeric)
PD.mn.array= zoopx2_ps #days by prey by strata by year!


#second version with quantiles ################################

PseudoQuantiles = filter(zoopsmwidef, Region %in% c("NE Suisun", "NW Suisun","SE Suisun")) %>%
  filter(!is.na(pdiapfor)) %>%
  group_by(doy) %>%
  summarize(pdiapfor_0.05 = quantile(pdiapfor, 0.05),pdiapjuv_0.05 = quantile(pdiapjuv, 0.05),
            pdiapfor_0.25 = quantile(pdiapfor, 0.25),pdiapjuv_0.25 = quantile(pdiapjuv, 0.25), 
            pdiapfor_0.50 = quantile(pdiapfor, 0.5),pdiapjuv_0.50 = quantile(pdiapjuv, 0.50),
            pdiapfor_0.75 = quantile(pdiapfor, 0.75),pdiapjuv_0.75 = quantile(pdiapjuv, 0.75),
            pdiapfor_0.95 = quantile(pdiapfor, 0.75),pdiapjuv_0.95 = quantile(pdiapjuv, 0.95)) 

Pseuodquantlong = pivot_longer(PseudoQuantiles, cols = c(pdiapfor_0.05:pdiapjuv_0.95), names_to = c("Lifestage", "Quantile"),
                               names_sep = "_") %>%
  pivot_wider(names_from = "Lifestage", values_from = value)

zoopquantiles = left_join(select(zoopsSuisun, -pdiapfor, -pdiapjuv, -level), Pseuodquantlong) 

#now the different temperature regimes

zoopsQS = bind_rows(mutate(zoopquantiles, Temp = "Temp_05"),
                    mutate(zoopquantiles, Temp = "Tempave"),
                    mutate(zoopquantiles, Temp = "Temp1"),
                    mutate(zoopquantiles, Temp = "Temp15"),
                    mutate(zoopquantiles, Temp = "Temp05"),
                    mutate(zoopquantiles, Temp = "Temp_1"),
                    mutate(zoopquantiles, Temp = "Temp_15"))  %>%
select(Temp, Quantile,doy, 
       limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
       daphnia, othcyc, other, eurytem, pdiapfor) 

test_psq= zoopsQS%>%
  split(list(zoopsQS$Temp, zoopsQS$Quantile))

zoop_psq= array(unlist(test_psq),dim=c(153,15,7, 5), 
               dimnames = list(c(153:305), names(zoopsQS), unique(zoopsQS$Temp),
                               unique(zoopsQS$Quantile)))


#OK! Now i just need to get rid of the id columns
zoopx_psq = zoop_psq[(1:153), c(4:15), c(1:7), c(1:5)]
zoopx2_psq = apply(zoopx_psq, c(2,3,4), as.numeric) #days by prey by temperature by quantile!



### 2.

#Different temperature regimes ######################
load("data/WaterQuality20102024.RData")

#switch this to mean daily temp
#Tempave = mean(filter(AllWQmean2, Parameter == "watertemperature", DOY %in% c(153:305))$ValueImputed)

Tempsave = AllWQmean2 %>%
  filter(DOY %in%c(153:305), , Parameter == "watertemperature",
         Region %in% c("NE Suisun", "NW Suisun","SE Suisun")) %>%
  group_by(DOY) %>%
  summarize(Temperature = mean(ValueImputed))
  
Tempwide_x = bind_rows(mutate(Tempsave, Temp = "Temp_15", Temperature = Temperature-1.5),
                       mutate(Tempsave, Temp = "Temp_1", Temperature = Temperature-1),
                       mutate(Tempsave, Temp = "Temp_05", Temperature = Temperature-0.5),
                       mutate(Tempsave, Temp = "Tempave"),
                       mutate(Tempsave, Temp = "Temp05", Temperature = Temperature+.5),
                       mutate(Tempsave, Temp = "Temp1", Temperature = Temperature+1),
                       mutate(Tempsave, Temp = "Temp15", Temperature = Temperature+1.5)) %>%
  pivot_wider(names_from = Temp, values_from = Temperature) 

# Tempwide_x = filter(AllWQmean2, Parameter == "watertemperature")%>%
#   select(Parameter, DOY, Region) %>%
#   distinct() %>%
#   mutate(ValueImputed = case_when(Region == "Lower Sacramento River" ~ Tempave,
#                                   Region == "Confluence" ~ Tempave-.5,
#                                   Region == "Suisun Marsh" ~ Tempave-1,
#                                   Region == "SW Suisun" ~ Tempave-1.5,
#                                   Region == "SE Suisun" ~ Tempave+.5,
#                                   Region == "NE Suisun" ~ Tempave+1,
#                                   Region == "NW Suisun" ~ Tempave+1.5),
#          Temp = case_when(Region == "Lower Sacramento River" ~ "Tempave",
#                           Region == "Confluence" ~ "Temp_05",
#                           Region == "Suisun Marsh" ~ "Temp_1",
#                           Region == "SW Suisun" ~ "Temp_15",
#                           Region == "SE Suisun" ~ "Temp05",
#                           Region == "NE Suisun" ~ "Temp1",
#                           Region == "NW Suisun" ~ "Temp15")) %>%
#   arrange(Temp) %>%
#    ungroup() %>%
#   select(-Region) %>%
#   pivot_wider(names_from = Temp, values_from = ValueImputed) 
 


Tempwidef_x = filter(Tempwide_x, DOY %in% c(153:305))
Tempwidef_x2 = bind_rows(mutate(Tempwidef_x, level = pseudosave$level[1]),
                         mutate(Tempwidef_x, level = pseudosave$level[2]),
                         mutate(Tempwidef_x, level = pseudosave$level[3]),
                         mutate(Tempwidef_x, level = pseudosave$level[4]),
                         mutate(Tempwidef_x, level = pseudosave$level[5]),
                         mutate(Tempwidef_x, level = pseudosave$level[6]),
                         mutate(Tempwidef_x, level = pseudosave$level[7]),
                         mutate(Tempwidef_x, level = pseudosave$level[8]),
                         mutate(Tempwidef_x, level = pseudosave$level[9]),
                         mutate(Tempwidef_x, level = pseudosave$level[10]),
                         mutate(Tempwidef_x, level = pseudosave$level[11]))

test_x = Tempwidef_x2  %>%
  split(list(Tempwidef_x2$level))

Temp_x = array(unlist(test_x),dim=c(153,9, 11), 
                      dimnames = list(c(153:305), names(Tempwidef_x2),
                                      unique(Tempwidef_x2$level)))

#OK! Now i just need to get rid of the id columns
Tempx_x = Temp_x[c(1:153), c(2:8), c(1:11)]
Tempx2_x = apply(Tempx_x, c(2,3), as.numeric)

#Quantile version #######

Tempwidef_x2q = bind_rows(mutate(Tempwidef_x, Quantile = "0.05"),
                         mutate(Tempwidef_x,Quantile = "0.25"),
                         mutate(Tempwidef_x, Quantile = "0.50"),
                         mutate(Tempwidef_x, Quantile = "0.75"),
                         mutate(Tempwidef_x, Quantile = "0.95"))

test_xq = Tempwidef_x2q  %>%
  split(list(Tempwidef_x2q$Quantile))

Temp_xq = array(unlist(test_xq),dim=c(153,9, 5), 
               dimnames = list(c(153:305), names(Tempwidef_x2q),
                               unique(Tempwidef_x2q$Quantile)))

#OK! Now i just need to get rid of the id columns
Tempx_xq = Temp_xq[c(1:153), c(2:8), c(1:5)]
Tempx2_xq = apply(Tempx_x, c(2,3), as.numeric)


#now turbidity (constant)


Turbave = mean(filter(AllWQmean2, Parameter == "turbidity", DOY %in% c(153:305),
                      Region %in% c("NE Suisun", "NW Suisun","SE Suisun"))$ValueImputed)

turbwide_x = filter(AllWQmean2, Parameter == "turbidity")%>%
  select(Parameter, DOY, Region) %>%
  distinct() %>%
  mutate(ValueImputed = Turbave,
         Temp = case_when(Region == "Lower Sacramento River" ~ "Tempave",
                          Region == "Confluence" ~ "Temp_05",
                          Region == "Suisun Marsh" ~ "Temp_1",
                          Region == "SW Suisun" ~ "Temp_15",
                          Region == "SE Suisun" ~ "Temp05",
                          Region == "NE Suisun" ~ "Temp1",
                          Region == "NW Suisun" ~ "Temp15")) %>%
  arrange(Temp) %>%
  ungroup() %>%
  select(-Region, -Parameter) %>%
  pivot_wider(names_from = Temp, values_from = ValueImputed) 


Turbwidef_x = filter(turbwide_x, DOY %in% c(153:305))
Turbwidef_x2 = bind_rows(mutate(Turbwidef_x, level = pseudosave$level[1]),
                         mutate(Turbwidef_x, level = pseudosave$level[2]),
                         mutate(Turbwidef_x, level = pseudosave$level[3]),
                         mutate(Turbwidef_x, level = pseudosave$level[4]),
                         mutate(Turbwidef_x, level = pseudosave$level[5]),
                         mutate(Turbwidef_x, level = pseudosave$level[6]),
                         mutate(Turbwidef_x, level = pseudosave$level[7]),
                         mutate(Turbwidef_x, level = pseudosave$level[8]),
                         mutate(Turbwidef_x, level = pseudosave$level[9]),
                         mutate(Turbwidef_x, level = pseudosave$level[10]),
                         mutate(Turbwidef_x, level = pseudosave$level[11]),)


test_xt = Turbwidef_x2  %>%
  split(list(Turbwidef_x2$level))

Turb_x = array(unlist(test_xt),dim=c(153,9, 11), 
               dimnames = list(c(153:305), names(Turbwidef_x2),
                               unique(Turbwidef_x2$level)))

#OK! Now i just need to get rid of the id columns
turbx_x = Turb_x[c(1:153), c(2:8), c(1:11)]
turbx2_x = apply(turbx_x, c(2,3), as.numeric)

obs.temp.dat <- Tempx2_x 
obs.turb.dat <- turbx2_x
daylight = daylength(38.15, c(130:310))$Daylength*60
LT.fx <- daylight/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace

#quantile version of turbidity ##############
Turbwidef_x2q = bind_rows(mutate(Turbwidef_x, Quantile = "0.05"),
                          mutate(Turbwidef_x,Quantile = "0.25"),
                          mutate(Turbwidef_x, Quantile = "0.50"),
                          mutate(Turbwidef_x, Quantile = "0.75"),
                          mutate(Turbwidef_x, Quantile = "0.95"))


test_xtq = Turbwidef_x2q  %>%
  split(list(Turbwidef_x2q$Quantile))

Turb_xq = array(unlist(test_xt),dim=c(153,9, 5), 
               dimnames = list(c(153:305), names(Turbwidef_x2q),
                               unique(Turbwidef_x2q$Quantile)))

#OK! Now i just need to get rid of the id columns
turbx_xq = Turb_x[c(1:153), c(2:8), c(1:5)]
turbx2_xq = apply(turbx_x, c(2,3), as.numeric)


### 3. Model bioenergetics  #############################

zooprun = smelt_bioenergetics(PD.mn.array = zoopx2_ps, obs.temp.dat = Tempx2_x, obs.turb.dat = turbx2_x, 
                              start.L = rep(23, 15), ex.strata =c("Temp_15","Temp_1","Temp_05","Tempave","Temp05","Temp1","Temp15"),
                              beta_hat = beta_hat[1:200,], start.year = 2010)


#relable the zooplankton levels correctly
unique(zooprun$Year)

zooprun = mutate(zooprun, Level = case_match(Year, 2010 ~ 0.01,2011 ~ 0.5, 2012 ~ 1, 2013 ~ 1.5, 2014 ~ 2, 2015 ~ 5,
                                             2016~ 10, 2017~20, 2018 ~ 40, 2019~ 100, 2020 ~ 200),
                 Month = case_match(Day, c(0:30) ~ 6,
                                    c(31:62) ~7,
                                    c(63:94) ~8,
                                    c(95:125) ~9,
                                    c(126:153) ~10)) %>%
  rename(Temp = Stratum)

zooprunsum = zooprun %>%
  group_by(Day, Month, Level, Temp) %>%
  summarize(MWeight = mean(Weight))

ggplot(zooprunsum, aes(x = Day, y = MWeight, color = as.factor(Temp))) +
  geom_smooth()+
  facet_wrap(~Level)
#weights never decrease, so they aren't starving. 

#OK, how much did they grow over the whole summer (or fall)?
Wtsum2_ps = zooprunsum %>%
  mutate(doy = Day+153) %>%
  left_join(select(Pseudo, doy, pdiapfor)) %>%
  group_by(Level, Temp) %>%
  summarize(startweight = first(MWeight), endwieght = last(MWeight), diffweight = endwieght-startweight) %>%
  left_join(rename(pseudosave, Level = level)) %>%
  mutate(Temp = factor(Temp, levels = c("Temp_15", "Temp_1", "Temp_05", "Tempave", 
                                        "Temp05", "Temp1", 
                                        "Temp15")),
         Tempnum = as.numeric(Temp))

#now do it keeping hte variability
Wtsum2_ps2 = zooprun %>%
  mutate(doy = Day+153) %>%
  left_join(select(Pseudo, doy, pdiapfor)) %>%
  group_by(Level, Temp, s) %>%
  summarize(startweight = first(Weight), endwieght = last(Weight), diffweight = endwieght-startweight) %>%
  group_by(Level, Temp) %>%
  summarize(meadidff = mean(diffweight), sdweight = sd(diffweight), minweight = mean(diffweight) - sdweight, max = mean(diffweight)+ sdweight) %>%
  left_join(rename(pseudosave, Level = level)) %>%
  mutate(Temp = factor(Temp, levels = c("Temp_15", "Temp_1", "Temp_05", "Tempave", 
                                        "Temp05", "Temp1", 
                                        "Temp15")),
         Tempnum = as.numeric(Temp)) 



ggplot(Wtsum2_ps, aes(x = Level, y = diffweight, fill = as.factor(Temp))) + geom_col(position = "dodge")+
  facet_wrap(~Level, scales = "free_x")+
  theme(axis.text.x = element_blank())

ggplot(Wtsum2_ps, aes(x = pdiapfor, y = diffweight, color = as.factor(Temp))) + geom_line(size =1)+
  ylab("Total summer/fall growth (g)")+ xlab("Average Pseudodiaptomus biomass (mg/m3)")+
  scale_color_brewer(palette = "Dark2", name = "Temperature\nscenario",
                     labels = c("Mean Temp - 1.5 C", "Mean Temp - 1 C", "Mean Temp - 0.5 C", 
                                "Mean Temp", "Mean Temp + 0.5 C", "Mean Temp + 1 C", "Mean Temp + 1.5 C"))+ theme_bw()+
  coord_cartesian(xlim = c(0,75))


ggsave("plots/zoop_saturattion_curv.tiff", device = "tiff", width =8, height =5)

#saturation curve with variance
ggplot(Wtsum2_ps, aes(x = pdiapfor, y = diffweight, color = as.factor(Temp))) + 
  geom_ribbon(data = Wtsum2_ps2, aes(x = pdiapfor, ymin = minweight, ymax = max, 
                                     fill = as.factor(Temp), color = as.factor(Temp)), 
              alpha = 0.2, inherit.aes = F, linetype =3, size = 0.25)+
  ylab("Growth (g) from \nJune 1 through Oct 31")+ xlab("Average Pseudodiaptomus biomass (mg/m3)")+ geom_line(size =1)+
  scale_color_viridis_d( name = "Temperature\nscenario",
                     labels = c("Mean Temp - 1.5 C", "Mean Temp - 1 C", "Mean Temp - 0.5 C", 
                                "Mean Temp", "Mean Temp + 0.5 C", "Mean Temp + 1 C", "Mean Temp + 1.5 C"))+ 
  scale_fill_viridis_d(name = "Temperature\nscenario",
                     labels = c("Mean Temp - 1.5 C", "Mean Temp - 1 C", "Mean Temp - 0.5 C", 
                                "Mean Temp", "Mean Temp + 0.5 C", "Mean Temp + 1 C", "Mean Temp + 1.5 C"))+theme_bw()+
  coord_cartesian(xlim = c(0,55), ylim = c(2,3.9))

ggsave("plots/zoop_saturattion_curv2.tiff", device = "tiff", width =8, height =5)

ggplot(Wtsum2_ps, aes(x = Level, y = diffweight, color = as.factor(Temp))) + geom_line()+
  ylab("Total summer/fall growth (g)")+ xlab("average pseudodiaptomus biomass")+
  scale_color_brewer(palette = "Dark2", name = "Temperature\nscenario")+ theme_bw()+
  coord_cartesian(xlim = c(0,100))


ggplot(Wtsum2_ps, aes(x = log(pdiapfor), y = diffweight, color = as.factor(Temp))) + geom_line()+ geom_point()

ggplot(Wtsum2_ps, aes(x = Tempnum, y = diffweight, color = as.factor(round(pdiapfor, 1))))+ geom_line()+
  ylab("total summer-fall growth (g)") + xlab("Temperature regime")+
  scale_color_brewer(palette = "Set3", name = "Zooplankton \nScenario")+
  theme_bw()

#######################################
#change in weight per change in zoops. 

test = filter(Wtsum2_ps, Temp == first(unique(Wtsum2_ps$Temp))) %>%
  ungroup() %>%
  arrange(pdiapfor) %>%
  mutate(lagzoop = lag(pdiapfor), diffzoop = pdiapfor-lagzoop, lagweight = lag(endwieght), 
         diffweigth = endwieght-lagweight,
         ratio = diffweight/(diffzoop*1000))


ggplot(test, aes(x = pdiapfor, y = ratio)) + geom_line()+
  ylab("Increase in weight per increase in zooplankton biomass")+
  xlab("average zooplankton biomass")+
  coord_cartesian(xlim = c(0,80))

#so, I think this is telling me that increasing pseudodiaptomus from 5 to 20 makes a big dofferemce.
#but increasing it from 20 to 1000 doesn't do that much for you. There is a threshold effect
#around about 15 or 20.
save(Wtsum2_ps, zooprun, zooprunsum, file = "outputs/increasedzoops.Rdata")
load("outputs/increasedzoops.Rdata")


#So, how mcuh zooplankton do we actually see? And how mcuh does flow actions increase it?
#this is juveniles and adults together, so not quite right, but still helpful
load("C:/Users/rhartman/OneDrive - California Department of Water Resources/salinity control gates/SFHA_synthesis/data/pseudo_data_mass.RData")

#need to divide mass by 1000 to convert to mg.

pseudo_recentx = mutate(pseudo_data_mass, BPUEmg = BPUE/1000, doy = yday(Date))

ggplot(filter(pseudo_recentx, doy %in% c(150:243), Year >2009), 
       aes(x = as.factor(Year), y = BPUEmg)) + geom_boxplot()+
  geom_hline(yintercept = 18)

ggplot(filter(pseudo_recentx, doy %in% c(150:243), Year >2009, Region == "Suisun Bay"), 
       aes(x = as.factor(Year), y = BPUEmg)) + geom_boxplot()+
  geom_hline(yintercept = 18)


#OK, so we are usually under the point of saturation.

ggplot(pseudo_recentx, aes(x = Month, fill = as.factor(Year), y = BPUEmg))+ geom_col(position = "dodge")+
  facet_wrap(~Region)#+
  #coord_cartesian(xlim = c(4,8))
#Egh, i don't think i back-transformed this properly

#look at average zoops somehow.
avePseudoAnnualSuisun = select(zoopsmAve, Date, Month, doy, Year, pdiapfor, Region, pdiapjuv) %>%
  filter(Region %in% c("NE Suisun", "NW Suisun","SE Suisun"), doy %in% c(153:305)) %>%
  group_by(Year) %>%
  summarize(MeanPseudo = exp(mean(log(pdiapfor), na.rm =T)), MeanPseudojuv = exp(mean(log(pdiapjuv), na.rm =T))) %>%
  mutate(totPseudo = MeanPseudo + MeanPseudojuv)


ggplot(filter(Wtsum2_ps, Temp == "Tempave"), aes(x = pdiapfor, y = diffweight)) + geom_line()+
  ylab("Total summer/fall growth (g)")+ xlab("Average Pseudodiaptomus Biomass (mg/L)")+theme_bw()+
  geom_vline(data = avePseudoAnnualSuisun, aes(xintercept = c(MeanPseudo + MeanPseudojuv), color = as.factor(Year)))+
  coord_cartesian(xlim = c(0,2), ylim = c(2, 3))+
  scale_color_manual(values = mypal)

ggplot(Wtsum2_ps,  aes(x = pdiapfor, y = diffweight, color = Temp)) + geom_line()+
  ylab("Total summer/fall growth (g)")+ xlab("Average Pseudodiaptomus Biomass (mg/L)")+theme_bw()+
  geom_vline(data = avePseudoAnnualSuisun, aes(xintercept = c(MeanPseudo + MeanPseudojuv)))+
  coord_cartesian(xlim = c(0,5), ylim = c(2, 3))#+
  #scale_color_manual(values = mypal)

Wtsum2_ave = filter(Wtsum2_ps, Temp == "Tempave")
Wtsum2_ps = left_join(Wtsum2_ps, pseudosave)
Wtsum2_ps2 = left_join(Wtsum2_ps2, pseudosave)


pcurve = ggplot(filter(Wtsum2_ps, Temp == "Tempave", pdiaptot <5),  aes(x = pdiaptot, y = diffweight)) + 
  geom_ribbon(data = filter(Wtsum2_ps2, Temp == "Tempave", pdiaptot <5),
              aes(x = pdiaptot, ymin = minweight, ymax = max), alpha = 0.2, inherit.aes = F, color = "grey")+
  geom_point()+geom_line(size =1)+
  #geom_text(aes(label = paste(Level*100, "%", sep = "")), nudge_x = 3, nudge_y = 0.05)+
 # geom_segment(data = avePseudoAnnualSuisun, aes(x = totPseudo, y = 2, yend = 2.1, color = as.factor(Year)),
 #              linewidth = 1)+ 
  annotate("text", x =.30, y = 2.38, label = "Average Pseudodiaptomus", vjust =0, size =3)+
  annotate("text", x =.30, y = 2.46, label = "50% more Pseudodiaptomus", vjust =0, size =3)+
  annotate("text", x =.30, y = 2.23, label = "50% less Pseudodiaptomus", vjust =0, size =3)+
  annotate("segment", x = 0, y = 2.44, xend = .972,linetype =2)+
  annotate("segment", x = 0, y = 2.35, xend = .648, linetype =2)+
  annotate("segment", x = 0, y = 2.21, xend = .324,  linetype =2)+
  
  annotate("text", x = 0, y = 2.82, label = "B", size =10)+
  ylab("Total summer/fall growth (g)")+ xlab(NULL)+
  theme_bw()+ coord_cartesian(xlim = c(0,2.5))
 # coord_cartesian(xlim = c(0,250), ylim = c(2, 2.75))+
#  scale_color_viridis_d(option = "turbo")
  
  # scale_x_continuous(breaks = c(0,50,100, 150, 200, 250), labels = c("0", "50%\n0.2mg/m3", "Average\n0.4mg/m3", "150%\n0.6mg/m3", "200%\n0.8mg/m3", "250%"))
pcurve
ggsave("plots/pseudo_percent_change_growth.tiff", device = "tiff", width =7, height =5)

#boxplot of pseudo biomass for context


pave = ggplot(avePseudoAnnualSuisun, aes(x =1, y = totPseudo)) + geom_boxplot()+coord_flip(ylim = c(0,2.5))+
  geom_point(aes(color = as.factor(Year)), size =3)+ xlab(NULL)+ theme_bw()+
  scale_color_viridis_d(option = "turbo", name = NULL)+
  
  annotate("text", x = 1.3, y = 0, label = "A", size =10)+
  theme(axis.text.y = element_blank(), legend.position = "inside", legend.position.inside = c(1, .8),
        legend.direction = "horizontal", legend.key.spacing =unit(.03, "cm"), legend.justification = 1)+ylab(NULL)
 #OK, what is the percent increase in pseudodiaptomus per X2? By month or day of year?
pave
save(Wtsum2_ps, Wtsum2_ps2, pseudo_recentx, avePseudoAnnualSuisun , file = "outputs/foodSatCurveData.RData")
load( "outputs/foodSatCurveData.RData")
########### change in mortality #######################
#Size-dependent mortality.—Mortality in the original baseline
##version was constant within each stage but decreased with successive stages, 
#so penalties in survival for slow growth occurred only through the delay in transition from 
#larvae to postlarvae and from postlarvae to juveniles. Making mortality length dependent 
#reflected the idea that vulnerability to predation mortality decreases with increasing size 
#(Sogard 1997; Bailey and DuffyAnderson 2010; Gislason et al. 2010), so that faster growth would 
#increase cumulative survival regardless of how stage transitions were triggered. We assumed that 
#mortality rate was a function of length (ML;d−1) for larvae through adults; we then fit the function to the 
#constant stage-specific mortality rates from the baseline simulation, associating the rate with the midpoint length of each stage:
#  ML =−0.034 + 0.165 · L−0.322.

#so going from weight to length would be
WttoL = function(Weight) {
  Length = (Weight/0.00000183)^(1/3.38)
return(Length)
}
WttoL(3)

zooprunsum_mort = zooprunsum %>%
  mutate(Length = WttoL(MWeight), Mortality = -0.034+0.165*Length^-0.322, Survived = 1-Mortality)

survival = zooprunsum_mort %>%
  group_by(Level, Temp) %>%
  summarize(Survived = prod(Survived)) %>%
  left_join(Wtsum2_ps2) %>%
  left_join(pseudosave)

#test just one sitution and make sure i'm doing it right
Test = filter(zooprunsum, Level ==1, Temp == "Tempave")%>%
  mutate(Length = WttoL(MWeight), Mortality = -0.034+0.165*Length^-0.322, Survived = 1-Mortality)

prod(Test$Survived)

write.csv(filter(zooprunsum, Temp == "Tempave"), "outputs/growthzoops.csv", row.names =F)

psurv = ggplot(filter(survival, Temp == "Tempave"), aes(x = pdiaptot, y = Survived)) +
  geom_point()+geom_line(size =1)+
  annotate("text", x =.25, y = 0.1245, label = "Average Pseudodiaptomus", vjust =0, size =3)+
  annotate("text", x =.25, y = 0.120, label = "50% less Pseudodiaptomus", vjust =0, size =3)+
  annotate("text", x =.25, y = 0.1275, label = "50% more Pseudodiaptomus", vjust =0, size =3)+
  annotate("segment", x = 0, y = 0.127, xend = .972,linetype =2)+
  annotate("segment", x = 0, y = 0.1241, xend = .648, linetype =2)+
  annotate("segment", x = 0, y = 0.1193, xend = .324,  linetype =2)+
  annotate("text", x = 0, y = 0.132, label = "C", size = 10)+
  ylab("Proportional survival Jun-Oct")+ xlab("Mean Pseudodiapotmus biomass (mg/m3)")+theme_bw()+
  coord_cartesian(xlim = c(0,2.50), ylim = c(0.105, 0.135))
  #scale_color_viridis_d(option = "turbo", name = "Year")+
  #scale_x_continuous(breaks = c(0,50,100, 150, 200, 250), labels = c("0", "50%\n0.2mg/m3", "Average\n0.4mg/m3", "150%\n0.6mg/m3", "200%\n0.8mg/m3", "250%"))
psurv
ggsave("plots/survivalcurve.tiff", device = "tiff", width = 6, height =4)


pave/pcurve/psurv+ plot_layout(heights = c(1,2,2))

ggsave("plots/pseudo_percent_change_wboxplot.tiff", device = "tiff", width =7, height =10)




###### If i increased zoops in suisun to levels in the river, what would happen? ##################
#Is that the right question? or is it looking at difference in pseudo abundance with X2?

#maybe i'll look at difference in growth with the increase in pseudo's that we'd expect with a change in X2 versus not.

#How about, for an average above normal year with average zoops, what's the increase in growth in just suisun bay
#caused by an X2 action, asuming X2 changes from 85 to 75?

#look at decrease in X2 in summer, then decrease in X2 in fall. Just suisun bay. 

#use my X2-pseudo relationship - object zm3
load("C:/Users/rhartman/OneDrive - California Department of Water Resources/salinity control gates/SFHA_synthesis/outputs/zoopx2model.RData")
library(glmmTMB)

#
X2s = data.frame(X2 = c(60:90))

#what was the actual x2 for each month?
X22024 = DF %>%
  filter(Year == 2024) %>%
  group_by(Month, Year) %>%
  summarize(X2 = mean(X2))

X22024 = data.frame(Region = "Suisun Bay",
  Scenario = c(rep("Base",5),  rep("X2 -4 Summer",5), rep("X2 -4 Fall", 5), rep("no X2",5)),
  Month = rep(c(6:10),4), X2 = c(72, 78, 79, 79, 79, 68, 74, 75, 79, 79,72, 78, 79, 75, 75, 72, 78, 79, 82, 82),
  Year = 2024)


Model_output = X22024 %>%
  bind_cols(Predicted_CPUE = predict(zm3, type = "response", newdata = X22024))

ggplot(Model_output, aes(x = Month, y = log(Predicted_CPUE), color = Scenario)) +
  
  geom_point()+ geom_line()



Scenarios2024 = Model_output %>%
  group_by(Scenario) %>%
  summarize(meanCPUE = exp(mean(log(Predicted_CPUE))))

#so it would go from 0.35 in the base case to 0.45 in the fall X2 case, to 0.50 in the summer x2 case. 

ggplot(filter(Wtsum2_ps, Temp == "Tempave"), aes(x = Level, y = diffweight)) + geom_line()+
  ylab("Total summer/fall growth (g)")+ xlab("Average Pseudodiaptomus Biomass (mg/L)")+theme_bw()+
  geom_vline(data = Scenarios2024, aes(xintercept = meanCPUE/1000, color = Scenario))+
  coord_cartesian(xlim = c(0,2), ylim = c(2, 3))+
  scale_color_manual(values = mypal)


#% change in pseudodiaptomus per change in X2 by month.
yrs = data.frame(Year = c(2011:2024))
regions = data.frame(Region = c("Suisun Bay", "Suisun Marsh", "Grizzly Bay", "River"))

newdat = data.frame(Month = rep(c(6:10), 25), X2 = rep(c(61:85), each =5)) %>%
  merge(yrs) %>%
  merge(regions)
newout = newdat %>%
  bind_cols(Predicted_CPUE = predict(zm3, type = "response", newdata = newdat)) %>%
  group_by(Month, Region) %>%
  mutate(lagCPUE = dplyr::lag(Predicted_CPUE), Percent_Change = lagCPUE/Predicted_CPUE) %>%
  ungroup()

#what's the percent change by region?
newout %>% 
  group_by(Region) %>%
  summarize(change = mean(Percent_Change, na.rm =T))



#% change in pseudodiaptomus per change in 4 km X2 by month.


newdat2 = data.frame(Month = rep(c(6:10), 7), X2 = rep(seq(60, 84, by =4), each =5)) %>%
  merge(yrs) %>%
  merge(regions)
newout2= newdat2 %>%
  bind_cols(Predicted_CPUE = predict(zm3, type = "response", newdata = newdat2)) %>%
  group_by(Month, Region) %>%
  mutate(lagCPUE = dplyr::lag(Predicted_CPUE), Percent_Change = lagCPUE/Predicted_CPUE) %>%
  ungroup()

#what's the percent change by region?
newout2 %>% 
  group_by(Region) %>%
  summarize(change = mean(Percent_Change, na.rm =T))



#####seasonal analysis############################



#average biomass of pseudodiaptomus by month in Suisun Bay
zoopsBay = filter(zoopsmwidef, Region %in% c("NE Suisun", "NW Suisun", "SE Suisun", "SW Suisun")) %>%
  group_by(doy) %>%
  summarize(across(limno:pdiapfor, function(x) mean(x, na.rm = TRUE))) %>%
  mutate(level =15,SeasonScenario = "Base") %>%
  mutate(Month = case_match(doy, c(153:181) ~ 6,
                            c(182:213) ~7,
                            c(214:245) ~8,
                            c(246:276) ~9,
                            c(277:305) ~10),
         Season = case_when(Month %in% c(6:7) ~ "Summer",
                            Month %in% c(8) ~ "Summer2",
                            Month %in% c(9:10) ~ "Fall"))
  
zoopsB = zoopsBay


#now replicate this for each biomass level
Levels = c(0.01, 0.5,1, 1.5, 2, 5, 10, 20, 40, 100, 200,500)
seasons = c("Base", "Summer", "Fall")
for(j in 1:3){
  for(i in 1:11){
    Level = Levels[i]
    Pseudo = zoopsBay %>%
      mutate(level = Level, 
             SeasonScenario =seasons[j],
             pdiapfor = case_when(Season == SeasonScenario ~ level*pdiapfor,
                                  TRUE ~ pdiapfor),
             pdiapjuv = case_when(Season == SeasonScenario ~ level*pdiapjuv,
                                  TRUE ~ pdiapjuv))
    zoopsB = bind_rows(zoopsB, Pseudo)
  }}


ggplot(zoopsB, aes(x = doy, y = pdiapfor))+
  facet_grid(SeasonScenario~level)+ geom_line()

zoopsB = filter(zoopsB, level != 15)


pseudosave_monthlyB = ungroup(zoopsB) %>%
  group_by(level, Month, Season, SeasonScenario) %>%
  summarize(pdiapfor = mean(pdiapfor, na.rm =T))

ggplot(pseudosave_monthlyB, aes(x = Month, y = pdiapfor, color = SeasonScenario))+ geom_line()+
  facet_wrap(~level, scales = "free")

pseudosave_monthly2B = ungroup(zoopsB) %>%
  group_by(Season) %>%
  summarize(pdiapfor = mean(pdiapfor, na.rm =T))



#organize correctly
zoopsMMB =  zoopsB  %>%
  select(Month, Season, SeasonScenario, level,doy, 
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) %>%
  arrange(level) %>%
  arrange(SeasonScenario) 

ggplot(zoopsMMB, aes(x = doy, y = pdiapfor, color = SeasonScenario)) +
  facet_wrap(~level, scales = "free") + geom_point()


test_psMB= zoopsMMB%>%
  split(list(zoopsMMB$SeasonScenario, zoopsMMB$level))

zoop_psMB= array(unlist(test_psMB),dim=c(153,17,3, 11), 
                dimnames = list(c(153:305), names(zoopsMMB), unique(zoopsMMB$SeasonScenario),
                                sort(unique(zoopsMMB$level))))


#OK! Now i just need to get rid of the id columns
zoopx_psMB = zoop_psMB[(1:153), c(6:17), c(1:3), c(1:11)]
zoopx2_psMB = apply(zoopx_psMB, c(2,3,4), as.numeric)#days by prey by strata by year!


### 2.

#Mean daily water temperature in Suisun Bay
load("data/WaterQuality20102024.RData")

Tempave = mean(filter(AllWQmean2, Parameter == "watertemperature", DOY %in% c(153:305))$ValueImputed)

SuisunAveTemp = filter(AllWQmean2, Parameter == "watertemperature", DOY %in% c(153:305), 
                       Region %in% c("NE Suisun", "NW Suisun", "SE Suisun", "SW Suisun")) %>%
  group_by(DOY) %>%
  summarize(Value = mean(ValueImputed, na.rm =T))

Tempwide_xB = bind_rows(mutate(SuisunAveTemp, SeasonScenario ="Base"),
                        mutate(SuisunAveTemp, SeasonScenario = "Summer"),
                        mutate(SuisunAveTemp, SeasonScenario ="Fall")) %>%
  pivot_wider(names_from = SeasonScenario, values_from = Value) 


Levels = unique(zoopsB$level)
Tempwidef_xB = filter(Tempwide_xB, DOY %in% c(153:305))

Tempwidef_x2B = bind_rows(mutate(Tempwidef_xB, level = Levels[1]),
                          mutate(Tempwidef_xB, level = Levels[2]),
                          mutate(Tempwidef_xB, level = Levels[3]),
                          mutate(Tempwidef_xB, level = Levels[4]),
                          mutate(Tempwidef_xB, level = Levels[5]),
                          mutate(Tempwidef_xB, level = Levels[6]),
                          mutate(Tempwidef_xB, level = Levels[7]),
                          mutate(Tempwidef_xB, level = Levels[8]),
                          mutate(Tempwidef_xB, level = Levels[9]),
                          mutate(Tempwidef_xB, level = Levels[10]),
                          mutate(Tempwidef_xB, level = Levels[11])) %>%
  arrange(level)

test_xB = Tempwidef_x2B  %>%
  split(list(Tempwidef_x2B$level))

Temp_xB = array(unlist(test_xB),dim=c(153,5, 11), 
                dimnames = list(c(153:305), names(Tempwidef_x2B),
                                unique(Tempwidef_x2B$level)))

#OK! Now i just need to get rid of the id columns
Tempx_xB = Temp_xB[c(1:153), c(2:4), c(1:11)]
Tempx2_xB = apply(Tempx_xB, c(2,3), as.numeric)



#now turbidity (constant)


Turbave = mean(filter(AllWQmean2, Parameter == "turbidity", DOY %in% c(153:305))$ValueImputed)

SuisunAveTurb = filter(AllWQmean2, Parameter == "turbidity", DOY %in% c(153:305), 
                       Region %in% c("NE Suisun", "NW Suisun", "SE Suisun", "SW Suisun")) %>%
  group_by(DOY) %>%
  summarize(Value = mean(ValueImputed))

turbwide_xB = bind_rows(mutate(SuisunAveTurb, SeasonScenario = "Base"),
                        mutate(SuisunAveTurb, SeasonScenario = "Summer"),
                        mutate(SuisunAveTurb,SeasonScenario = "Fall")) %>%
  pivot_wider(names_from = SeasonScenario, values_from = Value) 


Turbwidef_xB = filter(turbwide_xB, DOY %in% c(153:305))
Turbwidef_x2B = bind_rows(mutate(Turbwidef_xB, level = Levels[1]),
                          mutate(Turbwidef_xB, level =  Levels[2]),
                          mutate(Turbwidef_xB, level =  Levels[3]),
                          mutate(Turbwidef_xB, level =  Levels[4]),
                          mutate(Turbwidef_xB, level =  Levels[5]),
                          mutate(Turbwidef_xB, level =  Levels[6]),
                          mutate(Turbwidef_xB, level =  Levels[7]),
                          mutate(Turbwidef_xB, level =  Levels[8]),
                          mutate(Turbwidef_xB, level =  Levels[9]),
                          mutate(Turbwidef_xB, level =  Levels[10]),
                          mutate(Turbwidef_xB, level =  Levels[11])) %>%
  arrange(level)


test_xtB = Turbwidef_x2B %>%
  split(list(Turbwidef_x2B$level))

Turb_xB = array(unlist(test_xtB),dim=c(153,5, 11), 
                dimnames = list(c(153:305), names(Turbwidef_x2B),
                                unique(Turbwidef_x2B$level)))

#OK! Now i just need to get rid of the id columns
turbx_xB = Turb_xB[c(1:153), c(2:4), c(1:11)]
turbx2_xB = apply(turbx_xB, c(2,3), as.numeric)


daylight = daylength(38.15, c(130:310))$Daylength*60
LT.fx <- daylight/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace



zooprunSeasonal = smelt_bioenergetics(PD.mn.array = zoopx2_psMB, obs.temp.dat = Tempx2_xB, obs.turb.dat = turbx2_xB, 
                               start.L = rep(23, 15), ex.strata =c("Base","Summer","Fall"),
                               beta_hat = beta_hat[1:200,], start.year = 1)


#relable the zooplankton levels correctly
unique(zooprunSeasonal$Year)

zooprunSeasonal = mutate(zooprunSeasonal, Level = case_match(Year,1 ~ 0.01,  2 ~ 0.5, 3 ~ 1, 4 ~ 1.5, 5 ~ 2, 6 ~ 5,
                                               7~ 10, 8~20, 9 ~ 40, 10~ 100, 11 ~ 200),
                  Month = case_match(Day, c(0:30) ~ 6,
                                     c(31:62) ~7,
                                     c(63:94) ~8,
                                     c(95:125) ~9,
                                     c(126:153) ~10)) %>%
  rename(SeasonScenario = Stratum)

zooprunsumSeas = zooprunSeasonal  %>%
  group_by(Day,Level, SeasonScenario) %>%
  summarize(MWeight = mean(Weight)) %>%
  filter(MWeight >0)


ggplot(zooprunsumSeas, aes(x = Day, y = MWeight, color = SeasonScenario))+
  facet_wrap(~Level)+ geom_line()

ggplot(zooprunsumSeas, aes(x = Day, y = MWeight, color = as.factor(Level)))+
  facet_wrap(~SeasonScenario)+ geom_line()
#

aveps = group_by(pseudosave_monthlyB, level, SeasonScenario) %>%
  summarize(pdiapfor = mean(pdiapfor)) %>%
  rename(Level = level)
zooprun2maxS = ungroup(zooprunsumSeas) %>% 
  group_by(SeasonScenario, Level) %>%
  summarize(Weight = max(MWeight, na.rm =T)) %>%
  mutate(Levelx = factor(Level, 
                         labels = c("90% less", "50% less","Average", "50% more", "100% more",
                                    "500% more","10x more", "20 x more", "40 x more", "100 x more", "200 x more"))) %>%
  left_join(aveps)

ggplot(zooprun2maxS, aes(x = SeasonScenario, y = Weight)) +
  facet_wrap(~Levelx)+ geom_col()+
  #coord_cartesian(ylim = c(2.2, 3))+
  ylab("Final Weight at end of October")+
  xlab("Season where zooplankton increase")+
  coord_cartesian(ylim = c(1.8, 2.8))


ggplot(filter(zooprun2maxS, SeasonScenario != "Base"), 
       aes(x = Level*100, y = Weight, color = SeasonScenario)) + geom_line()+
  coord_cartesian(xlim = c(0,1500))+ xlab("Zooplankton increase - \nPercent of average")+
  ylab("Final weight at end of October (g)")+ 
  scale_color_manual(values = c("orange2", "green3"), name = "Season when zoops increase",
                     labels = c("Fall: Sep-Oct", "Summer: Jun-Jul"))
  

pseudosave_monthly2B = mutate(pseudosave_monthly2B, SeasonScenario = Season) %>%
  filter(Season != "Summer2")
ggplot(filter(zooprun2maxS, SeasonScenario != "Base"),
       aes(x = pdiapfor, y = Weight, color = SeasonScenario)) + geom_line()+
  coord_cartesian(xlim = c(0,25))+
  #geom_vline(data = pseudosave_monthly2B, aes(xintercept = pdiapfor, color = SeasonScenario))+
  ylab("Final weight at end of October (g)")+ xlab("average zooplankton density (mg/m3)") + 
  scale_color_manual(values = c("orange2", "green3"), name = "Season when zoops increase",
                     labels = c("Fall: Sep-Oct", "Summer: Jun-Jul"))





#########################################################################
## june gates versus september gates #############################


#average biomass of pseudodiaptomus by month in Suisun Marsh
zoopsmarsh = filter(zoopsmwidef, Region %in% c("Suisun Marsh")) %>%
  group_by(doy) %>%
  summarize(across(limno:pdiapfor, function(x) mean(x, na.rm = TRUE))) %>%
  mutate(level =1, MonthScenario = 0) %>%
  mutate(Month = case_match(doy, c(153:181) ~ 6,
                            c(182:213) ~7,
                            c(214:245) ~8,
                            c(246:276) ~9,
                            c(277:305) ~10))
zoopsM = zoopsmarsh

####different months of increase ########################
#Strata will be levels of zooplankton increase, years will be different months of zooplankton increase.
#i'll just use constant temperature

#now replicate this for each biomass level
Levels = c(0.01, 0.5, 1.5, 2, 5, 10, 20, 40, 100, 200,500)
for(j in 6:10){
for(i in 1:10){
  Level = Levels[i]
  Pseudo = zoopsmarsh %>%
    mutate(level = Level, 
           MonthScenario =j,
           pdiapfor = case_when(Month == j ~ level*pdiapfor,
                                TRUE ~ pdiapfor),
           pdiapjuv = case_when(Month == j ~ level*pdiapjuv,
                               TRUE ~ pdiapfor))
  zoopsM = bind_rows(zoopsM, Pseudo)
}}
ggplot(zoopsM, aes(x = doy, y = pdiapfor))+
  facet_grid(MonthScenario~level)+ geom_line()

zoopsM = filter(zoopsM, level != 1)


pseudosave_monthly = ungroup(zoopsM) %>%
  group_by(level, Month, MonthScenario) %>%
  summarize(pdiapfor = mean(pdiapfor, na.rm =T))

pseudosave_monthly2 = ungroup(zoopsmarsh) %>%
  group_by(Month) %>%
  summarize(pdiapfor = mean(pdiapfor, na.rm =T))

#organize correctly
zoopsMM =  zoopsM  %>%
  select(Month, MonthScenario, level,doy, 
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) %>%
  arrange(MonthScenario)

test_psM= zoopsMM%>%
  split(list(zoopsMM$MonthScenario, zoopsMM$level))

zoop_psM= array(unlist(test_psM),dim=c(153,16,5, 10), 
               dimnames = list(c(153:305), names(zoopsMM), unique(zoopsMM$MonthScenario),
                               unique(zoopsMM$level)))


#OK! Now i just need to get rid of the id columns
zoopx_psM = zoop_psM[(1:153), c(5:16), c(1:5), c(1:10)]
zoopx2_psM = apply(zoopx_psM, c(2,3,4), as.numeric)#days by prey by strata by year!


### 2.

#Mean daily water temperature in Suisun Marsh
load("data/WaterQuality20102024.RData")

Tempave = mean(filter(AllWQmean2, Parameter == "watertemperature", DOY %in% c(153:305))$ValueImputed)

SuisunAveTemp = filter(AllWQmean2, Parameter == "watertemperature", DOY %in% c(153:305), Region == "Suisun Marsh") %>%
  group_by(DOY) %>%
  summarize(Value = mean(ValueImputed))

Tempwide_xM = bind_rows(mutate(SuisunAveTemp, MonthScenario =6),
                        mutate(SuisunAveTemp, MonthScenario =7),
                        mutate(SuisunAveTemp, MonthScenario =8),
                        mutate(SuisunAveTemp, MonthScenario =9),
                        mutate(SuisunAveTemp, MonthScenario =10)) %>%
  pivot_wider(names_from = MonthScenario, values_from = Value) 


Levels = unique(zoopsM$level)
Tempwidef_xM = filter(Tempwide_xM, DOY %in% c(153:305))

Tempwidef_x2M = bind_rows(mutate(Tempwidef_xM, level = Levels[1]),
                         mutate(Tempwidef_xM, level = Levels[2]),
                         mutate(Tempwidef_xM, level = Levels[3]),
                         mutate(Tempwidef_xM, level = Levels[4]),
                         mutate(Tempwidef_xM, level = Levels[5]),
                         mutate(Tempwidef_xM, level = Levels[6]),
                         mutate(Tempwidef_xM, level = Levels[7]),
                         mutate(Tempwidef_xM, level = Levels[8]),
                         mutate(Tempwidef_xM, level = Levels[9]),
                         mutate(Tempwidef_xM, level = Levels[10]))

test_xM = Tempwidef_x2M  %>%
  split(list(Tempwidef_x2M$level))

Temp_xM = array(unlist(test_xM),dim=c(153,7, 10), 
               dimnames = list(c(153:305), names(Tempwidef_x2M),
                               unique(Tempwidef_x2M$level)))

#OK! Now i just need to get rid of the id columns
Tempx_xM = Temp_xM[c(1:153), c(2:6), c(1:10)]
Tempx2_xM = apply(Tempx_xM, c(2,3), as.numeric)



#now turbidity (constant)


Turbave = mean(filter(AllWQmean2, Parameter == "turbidity", DOY %in% c(153:305))$ValueImputed)

SuisunAveTurb = filter(AllWQmean2, Parameter == "turbidity", DOY %in% c(153:305), Region == "Suisun Marsh") %>%
  group_by(DOY) %>%
  summarize(Value = mean(ValueImputed))

turbwide_xM = bind_rows(mutate(SuisunAveTurb, MonthScenario =6),
                        mutate(SuisunAveTurb, MonthScenario =7),
                        mutate(SuisunAveTurb, MonthScenario =8),
                        mutate(SuisunAveTurb, MonthScenario =9),
                        mutate(SuisunAveTurb, MonthScenario =10)) %>%
  pivot_wider(names_from = MonthScenario, values_from = Value) 


Turbwidef_xM = filter(turbwide_xM, DOY %in% c(153:305))
Turbwidef_x2M = bind_rows(mutate(Turbwidef_xM, level = Levels[1]),
                         mutate(Turbwidef_xM, level =  Levels[2]),
                         mutate(Turbwidef_xM, level =  Levels[3]),
                         mutate(Turbwidef_xM, level =  Levels[4]),
                         mutate(Turbwidef_xM, level =  Levels[5]),
                         mutate(Turbwidef_xM, level =  Levels[6]),
                         mutate(Turbwidef_xM, level =  Levels[7]),
                         mutate(Turbwidef_xM, level =  Levels[8]),
                         mutate(Turbwidef_xM, level =  Levels[9]),
                         mutate(Turbwidef_xM, level =  Levels[10]))


test_xtM = Turbwidef_x2M  %>%
  split(list(Turbwidef_x2M$level))

Turb_xM = array(unlist(test_xtM),dim=c(153,7, 10), 
               dimnames = list(c(153:305), names(Turbwidef_x2M),
                               unique(Turbwidef_x2M$level)))

#OK! Now i just need to get rid of the id columns
turbx_xM = Turb_xM[c(1:153), c(2:6), c(1:10)]
turbx2_xM = apply(turbx_xM, c(2,3), as.numeric)


daylight = daylength(38.15, c(130:310))$Daylength*60
LT.fx <- daylight/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace



zooprun2 = smelt_bioenergetics(PD.mn.array = zoopx2_psM, obs.temp.dat = Tempx2_xM, obs.turb.dat = turbx2_xM, 
                              start.L = rep(23, 15), ex.strata =c("6","7","8","9","10"),
                              beta_hat = beta_hat[1:200,], start.year = 1)


#relable the zooplankton levels correctly
unique(zooprun2$Year)

zooprun2 = mutate(zooprun2, Level = case_match(Year,1 ~ 0.1,  2 ~ 0.5, 3 ~ 1.5, 4 ~ 2, 5 ~ 5, 6 ~ 10,
                                             7~ 20, 8~40, 9 ~ 100, 10~ 200),
                 Month = case_match(Day, c(0:30) ~ 6,
                                    c(31:62) ~7,
                                    c(63:94) ~8,
                                    c(95:125) ~9,
                                    c(126:153) ~10)) %>%
  rename(MonthScenario = Stratum)

zooprunsum2 = zooprun2 %>%
  group_by(Day, Month, Level, MonthScenario) %>%
  summarize(MWeight = mean(Weight)) %>%
  filter(MWeight >0)


ggplot(zooprunsum2, aes(x = Day, y = MWeight, color = MonthScenario))+
  facet_wrap(~Level)+ geom_line()

ggplot(zooprunsum2, aes(x = Day, y = MWeight, color = as.factor(Level)))+
  facet_wrap(~MonthScenario)+ geom_line()
#

aveps = filter(pseudosave_monthly, MonthScenario == Month) %>%
  mutate(MonthScenario = as.factor(MonthScenario)) %>%
  rename(Level = level)

#why isn't my june scenario working?
zooprun2max = group_by(zooprunsum2, MonthScenario, Level) %>%
  summarize(Weight = max(MWeight, na.rm =T)) %>%
  mutate(Levelx = factor(Level, 
                         labels = c("90% less", "50% less", "50% more", "100% more",
                                    "500% more","10x more", "20 x more", "40 x more", "100 x more", "200 x more"))) %>%
  left_join(aveps)

ggplot(zooprun2max, aes(x = MonthScenario, y = Weight)) +
  facet_wrap(~Levelx)+ geom_col()+
  coord_cartesian(ylim = c(2.2, 3))+
  ylab("Final Weight at end of October")+
  xlab("Month when zooplankton increase")


ggplot(zooprun2max, aes(x = Level, y = Weight, color = MonthScenario)) + geom_line()

pseudosave_monthly2 = mutate(pseudosave_monthly2, MonthScenario = as.factor(Month))
ggplot(zooprun2max, aes(x = pdiapfor, y = Weight, color = MonthScenario)) + geom_line()+
  coord_cartesian(xlim = c(0,50))+
  geom_vline(data = pseudosave_monthly2, aes(xintercept = pdiapfor, color = MonthScenario))

#Now, how much does salinity need to drop in each month to acheive a doubling of zoops?

####salinity per doubling zoops #####
load("outputs/salinityConversions_threeregions.RData")

#first what salinity gives you the average pseudo abundance?
pseudo_convert = filter(sal_conversions, IBMR == "pdiapfor", Region == "Suisun Marsh") %>%
  #mutate(across(starts_with("draw_"), ~exp(.x)-1))%>%
  pivot_longer(starts_with("draw_"), names_prefix="draw_", names_to="draw", values_to="fit")%>%
  mutate(fit=if_else(fit<0, 0, fit)) %>%
  group_by(Month, SalSurf) %>%
  summarize(pdiapfor = round(exp(mean(fit))/1000, 1))
  
  
pseudojuv_convert = filter(sal_conversions, IBMR == "pdiapjuv", Region == "Suisun Marsh")%>%

  pivot_longer(starts_with("draw_"), names_prefix="draw_", names_to="draw", values_to="fit")%>%
  mutate(fit=if_else(fit<0, 0, fit)) %>%
  group_by(Month, SalSurf) %>%
  summarize(pdiapjuv = round(exp(mean(fit))/1000, 1))

#try geometric mean instead
pseudosave_monthly3 = ungroup(zoopsmarsh) %>%
  group_by(Month) %>%
  summarize(pdiapfor = exp(mean(log(pdiapfor), na.rm =T)))

test = c(2,4,5,6,8,100)
mean(test)
exp(mean(log(test)))

#now see what salinity I am predicted to be at

salspredicted = left_join(mutate(pseudosave_monthly2, pdiapfor = round(pdiapfor, 1)), pseudo_convert)
#hm. I'm off the charts. 

#maybe i should use the linera model instead and do the same thing i did for X2. 
load("C:/Users/rhartman/OneDrive - California Department of Water Resources/salinity control gates/SFHA_synthesis/outputs/zoopx2model.RData")


##quantile bioenergetic run ###############################
pseudosaveQ = ungroup(zoopsQS) %>%
  group_by(Quantile) %>%
  summarize(pdiapfor = exp(mean(log(pdiapfor), na.rm =T)),
            pdiapjuv = exp(mean(log(pdiapjuv), na.rm =T))) %>%
            mutate(pdiaptot = pdiapfor+pdiapjuv)

zooprunq = smelt_bioenergetics(PD.mn.array = zoopx2_psq, obs.temp.dat = Tempx2_xq, obs.turb.dat = turbx2_xq, 
                              start.L = rep(23, 15), ex.strata =c("Temp_15","Temp_1","Temp_05","Tempave","Temp05","Temp1","Temp15"),
                              beta_hat = beta_hat[1:200,], start.year = 2010)


#relable the zooplankton levels correctly
unique(zooprunq$Year)

zooprunq = mutate(zooprunq, Quantile = case_match(Year, 2010 ~ 0.05,2011 ~ 0.25, 2012 ~0.5, 2013 ~ 0.75, 2014 ~ 0.95),
                 Month = case_match(Day, c(0:30) ~ 6,
                                    c(31:62) ~7,
                                    c(63:94) ~8,
                                    c(95:125) ~9,
                                    c(126:153) ~10)) %>%
  rename(Temp = Stratum)

zooprunsumq = zooprunq %>%
  group_by(Day, Month, Quantile, Temp) %>%
  summarize(MWeight = mean(Weight))

ggplot(zooprunsumq, aes(x = Day, y = MWeight, color = as.factor(Temp))) +
  geom_smooth()+
  facet_wrap(~Quantile)


#OK, how much did they grow over the whole summer (or fall)?

#now do it keeping hte variability
Wtsum2_ps2q = zooprunq %>%
  mutate(doy = Day+153) %>%
  left_join(select(Pseudo, doy, pdiapfor, pdiapjuv)) %>%
  group_by(Quantile, Temp, s) %>%
  summarize(startweight = first(Weight), endwieght = last(Weight), diffweight = endwieght-startweight) %>%
  group_by(Quantile, Temp) %>%
  summarize(meadidff = mean(diffweight), sdweight = sd(diffweight), minweight = mean(diffweight) - sdweight, max = mean(diffweight)+ sdweight) %>%
  left_join(mutate(pseudosaveQ, Quantile = as.numeric(Quantile))) %>%
  mutate(Temp = factor(Temp, levels = c("Temp_15", "Temp_1", "Temp_05", "Tempave", 
                                        "Temp05", "Temp1", 
                                        "Temp15")),
         Tempnum = as.numeric(Temp)) 


#saturation curve with variance
ggplot(Wtsum2_ps2q, aes(x = pdiaptot, y = meadidff, color = as.factor(Temp))) + 
  geom_ribbon(data = Wtsum2_ps2q, aes(x = pdiaptot, ymin = minweight, ymax = max, fill = as.factor(Temp), color = as.factor(Temp)), 
              alpha = 0.2, inherit.aes = F, linetype =3, size = 0.25)+
  ylab("Total summer/fall growth (g)")+ xlab("Average Pseudodiaptomus biomass (mg/m3)")+ geom_line(size =1)+
  scale_color_brewer(palette = "Dark2", name = "Temperature\nscenario",
                     labels = c("Mean Temp - 1.5 C", "Mean Temp - 1 C", "Mean Temp - 0.5 C", 
                                "Mean Temp", "Mean Temp + 0.5 C", "Mean Temp + 1 C", "Mean Temp + 1.5 C"))+ 
  scale_fill_brewer(palette = "Dark2", name = "Temperature\nscenario",
                    labels = c("Mean Temp - 1.5 C", "Mean Temp - 1 C", "Mean Temp - 0.5 C", 
                               "Mean Temp", "Mean Temp + 0.5 C", "Mean Temp + 1 C", "Mean Temp + 1.5 C"))+theme_bw()


#put them both together maybe? ######################

AllWtsum = bind_rows(Wtsum2_ps2, Wtsum2_ps2q) %>%
  filter(pdiaptot <60)

#saturation curve with variance
ggplot(AllWtsum, aes(x = pdiaptot, y = meadidff, color = as.factor(Temp), fill = as.factor(Temp))) + 
 # geom_ribbon(aes(x = pdiaptot, ymin = minweight, ymax = max, fill = as.factor(Temp), color = as.factor(Temp)), 
  #            alpha = 0.2, inherit.aes = F, linetype =3, size = 0.25)+
  ylab("Total summer/fall growth (g)")+ xlab("Average Pseudodiaptomus biomass (mg/m3)")+ geom_line(size =1)+
  scale_color_brewer(palette = "Dark2", name = "Temperature\nscenario",
                     labels = c("Mean Temp - 1.5 C", "Mean Temp - 1 C", "Mean Temp - 0.5 C", 
                                "Mean Temp", "Mean Temp + 0.5 C", "Mean Temp + 1 C", "Mean Temp + 1.5 C"))+ 
  scale_fill_brewer(palette = "Dark2", name = "Temperature\nscenario",
                    labels = c("Mean Temp - 1.5 C", "Mean Temp - 1 C", "Mean Temp - 0.5 C", 
                               "Mean Temp", "Mean Temp + 0.5 C", "Mean Temp + 1 C", "Mean Temp + 1.5 C"))+theme_bw()+
  coord_cartesian(xlim = c(0, 30))

#This is all jagged because the quantiles for pdiapjuv and ptiaptot may or may not hav elined up with ech other. 