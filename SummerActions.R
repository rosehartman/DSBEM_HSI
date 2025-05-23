#OK, not the results for the summer flow scenarios
#this is based on aanalysis from Will Smith, USFWS, as used in Smith and Nobriga, 2023
#It was adapted by Rosemary Hartman for use in the summer flow evaluation, spring 2025

library(tidyverse)
library(tidyverse)
library(abind)
library(zoo)
library(lubridate)
library(lattice)
library(ggplot2)
library(gridExtra)

library(arrayhelpers)
library(chillR)

#this is the function that actually runs the bioenergetics.
source("BEMfunction_summer.R")
### 1. Read model parameters and data6
n.strata <- 6 #"Confluence"   "NE Suisun"    "NW Suisun"    "SE Suisun"    "Suisun Marsh" "SW Suisun" 
ex.strata = sort(c( "NE Suisun", "SE Suisun", "NW Suisun",  "Confluence", "SW Suisun",  "Suisun Marsh"))
n.prey <- 12 #
n.days<-153 # summer and fall
n.scenarios = 9
ex.scenarios = sort(c("10800-3-m", "10800-60d", "10800-b",   "15000-3-m", "15000-60d", 
                      "15000-b",   "21000",     "Noaction",  "fallx2"))

#THe model was set up to run over years, not scenarios, so I recoded the years as scenarios.
yr.seq<- ex.scenarios
n.yrs<-9


#get the regional temperature quantiles in the right format
#Will set up the model to use multi-dimentional arrays as input, with days by region by year. 
#It is quite annoying, and it is quite easy to mess up.

tempRegional = read_csv("data/temperature_quantile_50_fitsubregion.csv") %>%
  select(day_of_year, `Confluence`,`NE Suisun`, `NW Suisun`,`SE Suisun`,    `Suisun Marsh`,`SW Suisun`)

#Quickly put it in a different format to make a quick pot
tempslong = pivot_longer(tempRegional, cols = c(`NE Suisun`, `SE Suisun`,`NW Suisun`, `Confluence`, `SW Suisun`,  `Suisun Marsh`),
                         names_to = "Region", values_to = "Temperature")

#check that it makes sense (it does)
ggplot(tempslong, aes(x = day_of_year, y = Temperature, color = Region)) + geom_line()

#temperature is the same for each year and scenario,
#so I set up a dtataframe that repeats the temperature data for each scenario.
#There's probably a better way to do this, but that would invovle rewriting a buch of code.

tempsall = bind_rows(mutate(tempRegional, Scenario = "10800-3-m"),
                     mutate(tempRegional, Scenario = "10800-60d"),
                     mutate(tempRegional, Scenario = "10800-b"),
                            mutate(tempRegional, Scenario = "15000-3-m"),
                     mutate(tempRegional, Scenario = "15000-60d"),
                     mutate(tempRegional, Scenario = "15000-b"),
                     mutate(tempRegional, Scenario = "21000"),
                     mutate(tempRegional, Scenario = "fallx2"),
                     mutate(tempRegional, Scenario = "Noaction"))%>%
  rename(DOY = day_of_year) %>%
  filter(DOY %in% c(153:305)) # just June-October

#Now I have to turn it into an array, that is set up as rows x regions x scnarios 
test = tempsall  %>%
  split(list(tempsall$Scenario))

TempsR = array(unlist(test),dim=c(153,8, 9), 
             dimnames = list(c(153:305), names(tempsall),
                             unique(tempsall$Scenario)))

#OK! Now i just need to get rid of the id columns
TempxR = TempsR[c(1:153), c(2:7), c(1:9)]

#and make it numeric
Tempx2R = apply(TempxR, c(2,3), as.numeric)


#Repeat the process for turbididty

#They gave me log-transformed turbidity, so also have to back-tranform it.
turbRegional = read_csv("data/ln_turbidity_quantile_50_fitsubregion.csv") %>%
  select(day_of_year, `Confluence`,`NE Suisun`, `NW Suisun`,`SE Suisun`,  `Suisun Marsh`, `SW Suisun`) %>%
  mutate(across(c(`Confluence`,`NE Suisun`, `NW Suisun`,`SE Suisun`,  `SW Suisun`,  `Suisun Marsh`), exp))

#reformat so I can graph it
turblong = pivot_longer(turbRegional, cols = c(`NE Suisun`, `SE Suisun`,`NW Suisun`, `Confluence`, `SW Suisun`,  `Suisun Marsh`),
                         names_to = "Region", values_to = "Turbidity")
#check that it makes sense
ggplot(turblong, aes(x = day_of_year, y = Turbidity, color = Region)) + geom_line()

#create a data frame that repeats the data for each scenario.
turball = bind_rows(mutate(turbRegional, Scenario = "10800-3-m"),
                     mutate(turbRegional, Scenario = "10800-60d"),
                     mutate(turbRegional, Scenario = "10800-b"),
                     mutate(turbRegional, Scenario = "15000-3-m"),
                     mutate(turbRegional, Scenario = "15000-b"),
                     mutate(turbRegional, Scenario = "15000-60d"),
                     mutate(turbRegional, Scenario = "21000"),
                     mutate(turbRegional, Scenario = "fallx2"),
                     mutate(turbRegional, Scenario = "Noaction"),)%>%
  rename(DOY = day_of_year) %>%
  filter(DOY %in% c(153:305))

test = turball %>%
  split(list(turball$Scenario))

TurbR = array(unlist(test),dim=c(153,8, 9), 
               dimnames = list(c(153:305), names(turball),
                               unique(turball$Scenario)))

#OK! Now i just need to get rid of the id columns
TurbxR = TurbR[c(1:153), c(2:7), c(1:9)]
Turbx2R = apply(TurbxR, c(2,3), as.numeric)

#################################################################
#now the zooplankton.
#I use the zooplankton that I adjusted with the salinity GAM resutls.  

load("outputs/SalinityStuff.RData")

#pull out just the median values from the results
#Ideally I"d like to include the uncertainties inthe zooplankton carried through to the 
#smelt model, but I don't have time to figure that out now.


zoopsmedian = pivot_longer(step4,  cols = c(`10800-3-m_median`:last_col()),
                           names_to = "Scenario", values_to = "Biomass") %>%
  filter(!str_detect(Scenario, "95"), Region %in% c("Confluence", "NW Suisun", "SW Suisun",
                                                    "NE Suisun", "East Suisun", "Suisun Marsh",
                                                    "SE Suisun")) 

zoopsmAveSummer = zoopsmedian%>%
  filter(DOY %in% c(152:304))%>%
  group_by(Region, Scenario, IBMR) %>%
  mutate(Biomass = na.locf(Biomass)) %>% #if we are missing data, carry last observation forward
  ungroup() %>%
  pivot_wider(id_cols = c(Region, DOY,  Scenario),
              names_from = IBMR, values_from = Biomass) %>%
  select(Region,  DOY,  Scenario, #put it in the right order for the model
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) %>%
  arrange(Scenario)

#create an array that is days x taxa x region x scenario
test = zoopsmAveSummer %>%
  split(list(zoopsmAveSummer$Region, zoopsmAveSummer$Scenario))

zoop = array(unlist(test),dim=c(153,15,6, 9), 
             dimnames = list(c(152:304), names(zoopsmAveSummer), unique(zoopsmAveSummer$Region),
                             sort(unique(zoopsmAveSummer$Scenario))))



#OK! Now i just need to get rid of the id columns
zoopx = zoop[c(1:153), c(4:15), c(1:6), c(1:9)]
zoopx2 = replace_na(apply(zoopx, c(2,3,4), as.numeric),0)

######################2024 version##############


zoopsmedian2024 = pivot_longer(step4_2024,  cols = c(`10800-3-m_median`:last_col()),
                           names_to = "Scenario", values_to = "Biomass") %>%
  filter(!str_detect(Scenario, "95"), Region %in% c("Confluence", "NW Suisun", "SW Suisun",
                                                    "NE Suisun", "East Suisun", "Suisun Marsh",
                                                    "SE Suisun")) 

zoopsmAveSummer2024 = zoopsmedian2024%>%
  filter(DOY %in% c(152:304))%>%
  group_by(Region, Scenario, IBMR) %>%
  mutate(Biomass = na.locf(Biomass)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Region, DOY,  Scenario),
              names_from = IBMR, values_from = Biomass) %>%
  select(Region,  DOY,  Scenario,
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) %>%
  arrange(Scenario)



test = zoopsmAveSummer2024 %>%
  split(list(zoopsmAveSummer2024$Region, zoopsmAveSummer2024$Scenario))

zoop2024 = array(unlist(test),dim=c(153,15,6, 9), 
             dimnames = list(c(152:304), names(zoopsmAveSummer2024), unique(zoopsmAveSummer2024$Region),
                             sort(unique(zoopsmAveSummer2024$Scenario))))



#OK! Now i just need to get rid of the id columns
zoopx2024 = zoop2024[c(1:153), c(4:15), c(1:6), c(1:9)]
zoopx22024 = replace_na(apply(zoopx2024, c(2,3,4), as.numeric),0)
PD.mn.array2024 = zoopx22024 #days by prey by strata by year!



save(zoopx2, zoopsmAveSummer,zoopx22024, zoopsmAveSummer2024, file = "zoopsSummerAction.RData")



################################################################
#day lengths
daylight = daylength(38.15, c(130:310))$Daylength*60
LT.fx <- daylight/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace

beta_hat <- read.table(file='data/beta_hat.txt') # MC filtering coefficients
#this is model uncertainty


### 3. Model bioenergetics 
#load("zoopsSummerAction.RData")

#This is with teh 2016 data
sfharun = smelt_bioenergetics_summer(PD.mn.array = zoopx2, obs.temp.dat = Tempx2R, obs.turb.dat = Turbx2R, 
                              start.L = rep(23, 15), ex.strata = ex.strata, beta_hat = beta_hat[1:200,], 
                              ex.scenarios = ex.scenarios)

#save(sfharun, file = "outputs/sfharun.RData")
#load("outputs/sfharun.RData")

#here is the 2024 data
sfharun2024 = smelt_bioenergetics_summer(PD.mn.array = zoopx22024, obs.temp.dat = Tempx2R, obs.turb.dat = Turbx2R, 
                                     start.L = rep(23, 9), ex.strata = ex.strata, beta_hat = beta_hat[1:200,], 
                                     ex.scenarios = ex.scenarios)

#save(sfharun2024, file = "outputs/sfharun2024.RData")

#put both datasets together and calculate the mean and standard deviation of the outputs
sfhasummary = bind_rows(mutate(sfharun, Year = "2016"),
                        mutate(sfharun2024, Year = "2024"))%>%
  group_by(Stratum, Scenario, Day, Year) %>%
  summarize(sdLength = sd(Length), sdWeight = sd(Weight),Length = mean(Length), Weight = mean(Weight)) %>%
  group_by(Stratum, Scenario, Year) %>%
  mutate(Growth = (Weight -dplyr::lag(Weight))/Weight) %>%
  ungroup()


ggplot(sfhasummary, aes(x = Day, y = Growth, color = Scenario)) + geom_line()+
  facet_wrap(Year~Stratum)

ggplot(filter(sfhasummary, Stratum == "SE Suisun", Year == "2016"), aes(x = Day, y = Growth, color = Scenario)) + 
  geom_line()+
  coord_cartesian(ylim = c(0.008, 0.015), xlim = c(75, 155))
#  facet_wrap(Year~Stratum)


ggplot(sfhasummary, aes(x = Day, y = Length, color = Scenario)) + geom_line()+
  facet_wrap(Year~Stratum)

ggplot(sfhasummary, aes(x = Day, y = Weight, color = Scenario)) + geom_line()+
  facet_wrap(~Stratum)

ggplot(filter(sfhasummary, Stratum == "Suisun Marsh", Day >100), aes(x = Day, y = Weight, color = Scenario)) + 
  geom_line()+
  facet_wrap(Year~Stratum)

#Now I want to see what the maximum length and weight they got was.
#since they all started out the same size, this tells me how much they grew.
totagrowth = sfhasummary %>%
  group_by(Stratum, Scenario, Year) %>%
  summarize(Length = max(Length, na.rm =T), Weight = max(Weight, na.rm =T), 
            sdLength = max(sdLength), sdWeight = max(sdWeight))%>%
  mutate(Scenario = factor(Scenario, levels = c("Noaction", "fallx2", "10800-3-m", "10800-60d", "10800-b",
                                                "15000-3-m", "15000-60d", "15000-b", "21000")))

##plot of legnth at the end of the summer ##############
ggplot(totagrowth, aes(x = Scenario, y = Length, group = Year, fill = Scenario)) + 
  geom_col(aes(alpha = Year), position = "dodge", color = "grey")+
  scale_alpha_manual(values = c(0.6, 1))+
  facet_wrap(~Stratum)+ coord_cartesian(ylim = c(52, 70))+
  geom_errorbar(aes(ymin = Length - sdLength, ymax = Length+sdLength), position = "dodge")+
  theme(axis.text.x = element_text(angle = 90))

##plot of weight at the end of the summer ##############
ggplot(totagrowth, aes(x = Scenario, y = Weight, group = Year, fill = Scenario)) + 
  geom_col(aes(alpha = Year), position = "dodge", color = "grey")+
  scale_alpha_manual(values = c(0.6, 1))+
  facet_wrap(~Stratum)+ coord_cartesian(ylim = c(1, 3))+
  geom_errorbar(aes(ymin = Weight - sdWeight, ymax = Weight+sdWeight), position = "dodge")+
  theme(axis.text.x = element_text(angle = 90))


save(sfharun, totagrowth,sfharun2024, sfhasummary, file = "outputs/sfharunbothyears.RData")
load("outputs/sfharunbothyears.RData")

###########################################################
#look at growth by month
totagrowth_month = sfhasummary %>%
  mutate(Month = case_match(Day, c(0:30) ~ 6,
                            c(31:62) ~ 7,
                            c(63:94) ~ 8,
                            c(95:125) ~ 9,
                            c(126:157) ~ 10)) %>%
  group_by(Stratum, Scenario, Year, Month) %>%
  summarize(Length = max(Length, na.rm =T), Weight = max(Weight, na.rm =T), 
            sdLength = max(sdLength), sdWeight = max(sdWeight))%>%
  mutate(Scenario = factor(Scenario, levels = c("Noaction", "fallx2", "10800-3-m", "10800-60d", "10800-b",
                                                "15000-3-m", "15000-60d", "15000-b", "21000")),
         Laggrowth = dplyr::lag(Weight), NewGrowth = Weight-Laggrowth,
         NewGrowth = case_when(Month ==6 ~ Weight - 0.07,
                               TRUE ~ NewGrowth))
  

ggplot(totagrowth_month, aes(x = Scenario, y = NewGrowth, color = Year)) +
  facet_grid(Month ~ Stratum, scales = "free_y")+
  geom_point(position = "dodge")+ ylab("Monthly weigth increase")+
  theme(axis.text.x = element_text(angle = 90))

###############################################################################################
#OK, now I need to add the salinities and calculate the growth rate in the 
#areas where salinity is <6

#combine both salinity files
salsummary = bind_rows(sal2016 %>%
  mutate(DOY = yday(Date), Year = "2016") %>%
  select(-Scenario2),
  sal2024 %>%
    mutate(DOY = yday(Date), Year = "2024") %>%
    select(-Scenario2, -Scenariox)) %>%
  rename(Stratum = Region)

#mass-specific growth rate in each region
growth_6 = sfhasummary %>%
  group_by(Year, Scenario) %>%
  mutate(lagweight = dplyr::lag(Weight), #lag of weigth to make it easeir to calculate growth
         DOY = Day +153, #convert "day" to "Day of year"
         Growth = (Weight-lagweight)/Weight,
         Growth2 = (Weight-lagweight)) %>% #grwoth rate g/g/day
  ungroup() %>%
  left_join(salsummary) %>%
  mutate(Good = case_when(Salinity <= 6 ~ TRUE,
                          TRUE ~ FALSE)) %>% #new variable to say whether salinity was <6
  filter(Good, Day !=0) %>% 
  group_by(Scenario, Year, Day) %>% 
  summarize(growth2 = mean(Growth2, na.rm =T),
            growth = mean(Growth, na.rm =T)) %>% #mean mass-specific grwoth rate in areas where salinity was <6
  group_by(Scenario, Year) %>%
  mutate(Weight = case_when(Day ==1 ~0.07329784)) %>% #starting weight was 0.07 (used in next step)
  ungroup()

#function to calculate growth from mass-specific growth rates and starting weight
growing = function(x) {
  x$Weight[1] = 0.07329784
  for(i in 2:max(x$Day)){
x$Weight[i] = x$Weight[i-1]*(1+x$growth[i])
  }
  return(x)
  }

growing2 = function(x) {
  x$Weight[1] = 0.07329784
  for(i in 2:max(x$Day)){
    x$Weight[i] = x$Weight[i-1]+x$growth2[i]
  }
  return(x)
}


#now apply the function to the growth rate data
test = growth_6 %>%
  group_by(Scenario, Year) %>%
  do(Weight = growing(.))

testx = growth_6 %>%
  group_by(Scenario, Year) %>%
  do(Weight = growing2(.))

#that didn't come out like i intended, this hould fix it.
test2 = bind_rows(test$Weight)
test2x= bind_rows(testx$Weight) %>%
  rename(Weight2 = Weight) %>%
  left_join(test2)
#plot the results
ggplot(test2x, aes(x = Day, y = growth, color = Scenario))+
  geom_line() + facet_wrap(~Year)+ylab("Growth Rate g/g/day")
#OK, well that's confusing. 

ggplot(test2x, aes(x = Day, y = growth2, color = Scenario))+
  geom_line() + facet_wrap(~Year)+ylab("Growth Rate g/day")+
  coord_cartesian(xlim = c(120, 155), ylim = c(0.018, 0.021))

ggplot(test2x, aes(x = Day, y = Weight, color = Scenario))+
  geom_line() + facet_wrap(~Year)+ylab("Weight (g)")

ggplot(test2x, aes(x = Day, y = Weight2, color = Scenario))+
  geom_line() + facet_wrap(~Year)+ylab("Weight (g)")


#calculate mean growth rate and maximum weigth if fish are restricted to
#salinities of <6

#this is what does into the consequence table.
meangrowth = group_by(test2x, Year, Scenario) %>%
  summarise(Growth = mean(growth), Growth2 = mean(growth2),total = max(Weight, na.rm =T),
            total2 = max(Weight2, na.rm =T))%>%
  mutate(Scenario = factor(Scenario, levels = c("Noaction", "fallx2", "10800-3-m",
                                                "10800-60d", "10800-b", "15000-3-m", "15000-60d",
                                                "15000-b", "21000")))

write.csv(meangrowth, "outputs/meangrowth_summeractions.csv")

ggplot(meangrowth, aes(x = Scenario, y = total, fill = Scenario, group = Year)) + 
  geom_col(aes(alpha = Year), position = "dodge")+
  scale_alpha_manual(values = c(0.7, 1))+
  coord_cartesian(ylim = c(1.5,2.1))+
  ylab("total summer growth (g)")

ggplot(meangrowth, aes(x= Scenario, y = total2, fill = Scenario, group = Year)) + 
  geom_col(aes(alpha = Year), position = "dodge")+
  scale_alpha_manual(values = c(0.7, 1))+
  coord_cartesian(ylim = c(1.5,2.3))+
  ylab("total summer growth (g)")


#####################################################################
#why is X2 lower than no action?
#Maybe growth is higher in the confluence and they are restricted tot the confluence in no-action,
#but are in the Bay with lower food in the Action scenario?

growth_6a = sfhasummary %>%
  group_by(Year, Scenario) %>%
  mutate(lagweight = dplyr::lag(Weight), DOY = Day +153,
         Growth = (Weight-lagweight)) %>%
  ungroup() %>%
  left_join(salsummary) %>%
  mutate(Good = case_when(Salinity <= 6 ~ TRUE,
                          TRUE ~ FALSE)) %>%
  filter(DOY>153)

ggplot(filter(growth_6a, Scenario %in% c("Noaction", "fallx2")),
       aes(x = DOY, y = Stratum))+ geom_tile(aes(color = Good, fill = Growth))+
  scale_alpha_manual(values = c(0.5, 1))+
  facet_wrap(Scenario~Year, nrow =2)
 

test = growth_6a %>%
  filter(!Stratum  %in% c("SE Suisun", "NE Suisun", "NW Suisun")) %>%
  group_by(Year, Scenario) %>%
  mutate(lagweight = dplyr::lag(Weight), DOY = Day +153,
         Growth = (Weight-lagweight)) %>%
  ungroup() %>%
  left_join(salsummary) %>%
  mutate(Good = case_when(Salinity <= 6 ~ TRUE,
                          TRUE ~ FALSE)) %>%
  filter(DOY>153) %>% 
  group_by(Scenario, Year, Day) %>% 
  summarize(growth = mean(Growth, na.rm =T)) %>% #mean mass-specific grwoth rate in areas where salinity was <6
  group_by(Scenario, Year) %>%
  mutate(Weight = case_when(Day ==1 ~0.07329784)) %>% #starting weight was 0.07 (used in next step)
  ungroup()



#now apply the function to the growth rate data
test2 = test %>%
  group_by(Scenario, Year) %>%
  do(Weight = growing(.))

#that didn't come out like i intended, this hould fix it.
test2a = bind_rows(test2$Weight)

#plot the results
ggplot(filter(test2a,Day >1), aes(x = Day, y = growth, color = Scenario))+
  geom_line() + facet_wrap(~Year)+ylab("Growth Rate g/g/day")
#OK, well that's confusing. 

ggplot(filter(test2a,Day >1), aes(x = Day, y = Weight, color = Scenario))+
  geom_line() + facet_wrap(~Year)+ylab("Weight (g)")


#calculate mean growth rate and maximum weigth if fish are restricted to
#salinities of <6

#this is what does into the consequence table.
meangrowth2 = group_by(test2a, Year, Scenario) %>%
  summarise(Growth = mean(growth), total = max(Weight, na.rm =T)) %>%
  mutate(Scenario = factor(Scenario, levels = c("Noaction", "fallx2", "10800-3-m",
                                                "10800-60d", "10800-b", "15000-3-m", "15000-60d",
                                                "15000-b", "21000")))

ggplot(meangrowth2, aes(x = Scenario, y = total, fill = Scenario, group = Year)) + 
  geom_col(aes(alpha = Year), position = "dodge")+
  scale_alpha_manual(values = c(0.7, 1))+
  coord_cartesian(ylim = c(.5,.7))+
  ylab("total summer growth (g)")

#OK! Fall X2 is higher in this situation. So i'd probably a matter of Suisun declining in suitablity

##################################################################

growth_6_A = sfhasummary %>%
  mutate(DOY = Day+154) %>%
  left_join(salsummary) %>%
  mutate(Good = case_when(Salinity <= 6 ~ TRUE,
                          TRUE ~ FALSE)) %>%
  filter(DOY>154, Good) %>%
  group_by(Year, Scenario, DOY)  %>%
  summarize(Weight = mean(Weight))


ggplot(growth_6_A, aes(x = DOY, y = Weight, color = Scenario)) + geom_line()+
  facet_wrap(~Year)


#########################################################################
#what's going on with Suisun marsh?
#it looks like actions do worse than no actions, which doesn't make sesne because 
#zooplankton increse with actions.
ggplot(filter(zoopsmedian, Region == "Suisun Marsh"), aes(x = DOY, y = Biomass, fill = (IBMR)))+ geom_col() +
  facet_wrap(~Scenario)

ggplot(filter(zoopsmedian, Region == "NW Suisun"), aes(x = DOY, y = Biomass, fill = (IBMR)))+ geom_col() +
  facet_wrap(~Scenario)

ggplot(filter(zoopsmedian2024, Region == "Suisun Marsh"), aes(x = DOY, y = Biomass, fill = (IBMR)))+ geom_col() +
  facet_wrap(~Scenario)

#set it up with just Suisun Marsh

zoopex = abind(zoopx22024[,,,], zoopx2[,,,], along =3)
zoopex = zoopex[,,c(6,12),]
dimnames(zoopex)[[3]] = c("2016", "2024")

tempex = abind(Tempx2R, Tempx2R, along =2)
tempex = tempex[,c(6,12),]
dimnames(tempex)[[2]] = c("2016", "2024")

turbex = abind(Turbx2R, Turbx2R, along =2)
turbex = turbex[,c(6,12),]

dimnames(turbex)[[2]] = c("2016", "2024")

suisunrun = smelt_bioenergetics_summer(PD.mn.array = zoopex, 
                                                                 obs.temp.dat = tempex, 
                                                                  obs.turb.dat = turbex, 
                                                                  start.L = rep(23, 15), 
                                                                  beta_hat = beta_hat[1,], 
                                       ex.strata = c("2016","2024"),
                                                                  ex.scenarios = ex.scenarios)

suisun =  suisunrun%>%
  group_by(Stratum, Scenario, Day) %>%
  summarize(sdLength = sd(Length), sdWeight = sd(Weight),Length = mean(Length), Weight = mean(Weight))

ggplot(suisun, aes(x = Day, y = Length, color = Scenario)) + geom_line()+
  facet_wrap(~Stratum)+
  coord_cartesian(xlim = c(140,155), ylim = c(53,55))+
  scale_color_brewer(palette = "Set3")


