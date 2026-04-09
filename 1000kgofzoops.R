#extra zooplankton growth
#how much growth would we get if we dumped an extra 1000kg of zoops into the marsh?

#how much for 6000 kg of zoops?

#do I want to get into the fact that the different actions will give you different types of zoops? Probably. 

#copepods have about a 2-week life span, so I'll divide the time into two week intervals and add some in each interval

#5 months, so 10, two-week intervals. 100kg per fortnight.

#I'll just do Suisun Marsh to start.

# 1. Read model parameters and data
# 3. Model bioenergetics and reference points
# 4. Get summaries
# 5. Graphs
library(tidyverse)
library(abind)
library(zoo)
library(lubridate)

library(ggplot2)

library(gsw)
library(chillR)
library(arrayhelpers)
source("BEMfunction.R")

#

load("data/zoopsmwide.RData")

zoopsmwidef = filter(zoopsmAve, doy %in% c(153:305))%>%
  arrange(Region)

#average biomass of pseudodiaptomus by month in Suisun Marsh - all years - 2010-2014
zoopsMarsh = filter(zoopsmwidef, Region %in% c("Suisun Marsh"))  %>%
  mutate(Taxa = "Baseline")
zoopsS = zoopsMarsh

#how do I add 1000Kg?
#first calculate total biomass in the whole marsh, add 1000kg, then divide by volume again.

#volume of Suisun Marsh is 35442275 m3
#1000kg/35442275

1000/35442275

#convert to mg
1000/35442275*1000000

# increase of 28.2 mg over the entire summer. Let's assume we add some ever 
#two weeks, and they all live about 2 weeks. So 2.82 mg increase per day.
#I'll also go through and see how this would be different for each taxon

Taxa = names(zoopsMarsh)[6:17]

for(i in 1:12){
  Pseudo = zoopsMarsh
  Pseudo[,Taxa[i]]  = Pseudo[,Taxa[i]]+2.82
  Pseudo$Taxa = Taxa[i]
  zoopsS = bind_rows(zoopsS, Pseudo)
}


#now the different years

zoops1000 = zoopsS %>%
  mutate(Taxa = factor(Taxa, levels = c("Baseline", "limno","othcaljuv","pdiapjuv","othcalad" ,"acartela",
                                        "othclad",    "allcopnaup",
                                        "daphnia",    "othcyc",     "other",      "eurytem",    "pdiapfor"))) %>%
  select(Year, Taxa,doy, 
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) %>%
  arrange(Taxa)

test_1000= zoops1000%>%
  split(list(zoops1000$Taxa, zoops1000$Year))

zoop_1000= array(unlist(test_1000),dim=c(153,15, 13,15), 
                dimnames = list(c(153:305), names(zoops1000),
                                unique(zoops1000$Taxa), unique(zoops1000$Year)))


#OK! Now i just need to get rid of the id columns
zoopx_1000 = zoop_1000[(1:153), c(4:15), c(1:13), c(1:15)] #Days by taxa by taxon increase by year
zoopx2_1000 = apply(zoopx_1000, c(2,3,4), as.numeric) 



######################################

#OK, now temperature and turbidity
load("data/WaterQuality20102024.RData")


TempsS = AllWQmean2 %>%
  filter(DOY %in%c(153:305), , Parameter == "watertemperature",
         Region == "Suisun Marsh") 

Tempwide_x = bind_rows(mutate(TempsS,Taxon = "Baseline"),
  mutate(TempsS,Taxon = Taxa[1]),
                       mutate(TempsS,Taxon = Taxa[2]),
                       mutate(TempsS,Taxon = Taxa[3]),
                       mutate(TempsS,Taxon = Taxa[4]),
                       mutate(TempsS,Taxon = Taxa[5]),
                       mutate(TempsS,Taxon = Taxa[6]),
                       mutate(TempsS,Taxon = Taxa[7]),
                       mutate(TempsS,Taxon = Taxa[8]),
                       mutate(TempsS,Taxon = Taxa[9]),
                       mutate(TempsS,Taxon = Taxa[10]),
                       mutate(TempsS,Taxon = Taxa[11]),
                       mutate(TempsS,Taxon = Taxa[12])) %>%
  ungroup() %>%
pivot_wider(names_from = Taxon, values_from =ValueImputed) %>%
  select(Year,DOY, Baseline,
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor)


test_x = Tempwide_x  %>%
  split(list(Tempwide_x$Year))

Temp_x = array(unlist(test_x),dim=c(153,15, 15), 
               dimnames = list(c(153:305), names(Tempwide_x),
                               unique(Tempwide_x$Year)))

#OK! Now i just need to get rid of the id columns
Tempx_x = Temp_x[c(1:153), c(3:15), c(1:15)]
Tempx2_x = apply(Tempx_x, c(2,3), as.numeric)

##now turbidity #################################


TurbS = AllWQmean2 %>%
  filter(DOY %in%c(153:305), , Parameter == "turbidity",
         Region == "Suisun Marsh") 

TurbS_x = bind_rows(mutate(TurbS,Taxon = "Baseline"),mutate(TurbS,Taxon = Taxa[1]),
                       mutate(TurbS,Taxon = Taxa[2]),
                       mutate(TurbS,Taxon = Taxa[3]),
                       mutate(TurbS,Taxon = Taxa[4]),
                       mutate(TurbS,Taxon = Taxa[5]),
                       mutate(TurbS,Taxon = Taxa[6]),
                       mutate(TurbS,Taxon = Taxa[7]),
                       mutate(TurbS,Taxon = Taxa[8]),
                       mutate(TurbS,Taxon = Taxa[9]),
                       mutate(TurbS,Taxon = Taxa[10]),
                       mutate(TurbS,Taxon = Taxa[11]),
                       mutate(TurbS,Taxon = Taxa[12])) %>%
  ungroup() %>%
  pivot_wider(names_from = Taxon, values_from =ValueImputed) %>%
  select(Year,DOY, Baseline,
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor)


test_xt = TurbS_x  %>%
  split(list(TurbS_x$Year))

Turb_x = array(unlist(test_xt),dim=c(153,15, 15), 
               dimnames = list(c(153:305), names(TurbS_x),
                               unique(TurbS_x$Year)))

#OK! Now i just need to get rid of the id columns
Turbx_x = Turb_x[c(1:153), c(3:15), c(1:15)]
Turbx2_x = apply(Turbx_x, c(2,3), as.numeric)

### 3. Model bioenergetics  #############################

zooprun1000 = smelt_bioenergetics(PD.mn.array = zoopx2_1000 , obs.temp.dat = Tempx2_x, 
                                  obs.turb.dat = Turbx2_x, 
                              start.L = rep(23, 15), 
                              ex.strata = names(TurbS_x)[3:15],
                              beta_hat = beta_hat[1:200,], start.year = 2010)


#relable the zooplankton levels correctly
unique(zooprun1000$Year)
unique(zooprun1000$Stratum)

zooprun1000 = zooprun1000%>%
  rename(Taxa = Stratum)


#ggplot(zooprun1000, aes(x = Day, y = Weight, color = as.factor(Year)))+
#  facet_wrap(~Taxa)+
#  geom_smooth()



Wt.summary1000 = zooprun1000 %>%
  group_by(Taxa, Year, Day) %>%
  summarize(MWeight = mean(Weight, na.rm =T), sdWeight = sd(Weight, na.rm =T),
            MLength = mean(Length, na.rm =T), sdLength = sd(Length, na.rm =T)) 


ggplot(Wt.summary1000 , aes(x = Day, y = MWeight,  group = as.factor(Year), color = as.factor(Year))) +
  geom_smooth()+
  facet_wrap(~Taxa)+
  ylab("Predicted Weight in Grams")

ggplot(Wt.summary1000 , aes(x = Day, y = MWeight,  color = Taxa)) +
  geom_smooth()+
  facet_wrap(~Year)+
  ylab("Predicted Weight in Grams")+
  scale_color_manual(values = c("black", "slategrey", "slategray1", "slategray2", "skyblue", "skyblue2", "lightgreen",
                                "green3", "green4", "cyan2", "cyan3", "pink", "pink3"))





ggplot(filter(Wt.summary1000, Taxa %in% c("Baseline", "eurytem", "pdiapfor", "limno")),
              aes(x = Day, y = MWeight,  color = Taxa)) +
  geom_smooth()+
  facet_wrap(~Year)+
  ylab("Predicted Weight in Grams")+
  scale_color_manual(values = c("black", "slategrey", "slategray1", "slategray2", "skyblue", "skyblue2", "lightgreen",
                                "green3", "green4", "cyan2", "cyan3", "pink", "pink3"))

#OK, this makes sense. Taxa with the same consumption rate raise the growth rate the same amount

#difference to baseline ##############################

baseline = filter(Wt.summary1000, Taxa == "Baseline")  %>%
  rename(BaseWeight = MWeight) %>%
  ungroup() %>%
  select(Day, Year, BaseWeight)

Wt.summary1000 = left_join(Wt.summary1000, baseline) %>%
  mutate(Diff = MWeight-BaseWeight)

ggplot(Wt.summary1000, aes(x = Day, y = Diff,  color = Taxa)) +
  geom_smooth()+
  facet_wrap(~Year)+
  ylab("difference in weight from baseline in Grams")+
  scale_color_manual(values = c("black", "slategrey", "slategray1", "slategray2", "skyblue", "skyblue2", "lightgreen",
                                "green3", "green4", "cyan2", "cyan3", "pink", "pink3"))

#values at the end

maxweights =Wt.summary1000 %>%
  group_by(Year, Taxa) %>%
  summarise(Diff = max(Diff))

ggplot(maxweights, aes(x = as.factor(Year), y = Diff)) + geom_boxplot()+
  ylab("Increase in weight with food subsidy")

#average growth rate for all the critters smelt eat
filter(ungroup(maxweights), Diff !=0) %>%
  summarise(mean(Diff))

#maximum possible - use for zoop tanks
filter(ungroup(maxweights), Diff !=0, Taxa == "pdiapfor") %>%
  summarise(mean(Diff))
