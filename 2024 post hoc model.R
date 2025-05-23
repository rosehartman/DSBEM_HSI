#2024 post-hoc analysis for ITP ammendment


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
n.scenarios = 4
ex.scenarios = sort(c("NoAction",  "X2Only",    "GatesOnly", "Gates+X2"))

#THe model was set up to run over years, not scenarios, so I recoded the years as scenarios.
yr.seq<- ex.scenarios
n.yrs<-4


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

tempsall24 = bind_rows(mutate(tempRegional, Scenario = "Gates+X2"),
                     mutate(tempRegional, Scenario = "GatesOnly"),
                     mutate(tempRegional, Scenario = "NoAction"),
                     mutate(tempRegional, Scenario = "X2Only"))%>%
  rename(DOY = day_of_year) %>%
  filter(DOY %in% c(153:305)) # just June-October

#Now I have to turn it into an array, that is set up as rows x regions x scnarios 
test = tempsall24  %>%
  split(list(tempsall24$Scenario))

TempsR24 = array(unlist(test),dim=c(153,8, 4), 
               dimnames = list(c(153:305), names(tempsall24),
                               unique(tempsall24$Scenario)))

#OK! Now i just need to get rid of the id columns
TempxR24 = TempsR24[c(1:153), c(2:7), c(1:4)]

#and make it numeric
Tempx2R24 = apply(TempxR24, c(2,3), as.numeric)


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
turball24 = bind_rows(mutate(turbRegional, Scenario = "Gates+X2"),
                    mutate(turbRegional, Scenario = "GatesOnly"),
                    mutate(turbRegional, Scenario = "NoAction"),
                    mutate(turbRegional, Scenario = "X2Only"))%>%
  rename(DOY = day_of_year) %>%
  filter(DOY %in% c(153:305))

test = turball24 %>%
  split(list(turball24$Scenario))

TurbR24 = array(unlist(test),dim=c(153,8, 4), 
              dimnames = list(c(153:305), names(turball24),
                              unique(turball24$Scenario)))

#OK! Now i just need to get rid of the id columns
TurbxR24 = TurbR24[c(1:153), c(2:7), c(1:4)]
Turbx2R24 = apply(TurbxR24, c(2,3), as.numeric)

#################################################################
#now the zooplankton.
#I use the zooplankton that I adjusted with the salinity GAM resutls.  

step424 = read_csv("Outputs/zoop sal adjustments_2024gates.csv")

#pull out just the median values from the results
#Ideally I"d like to include the uncertainties inthe zooplankton carried through to the 
#smelt model, but I don't have time to figure that out now.


zoopsmedian24 = pivot_longer(step424,  cols = c(`Gates+X2_median`:last_col()),
                           names_to = "Scenario", values_to = "Biomass") %>%
  filter(!str_detect(Scenario, "95"), region %in% c("Confluence", "NW Suisun", "SW Suisun",
                                                    "NE Suisun", "East Suisun", "Suisun Marsh",
                                                    "SE Suisun")) 

zoopsmAveSummer24 = zoopsmedian24%>%
  filter(DOY %in% c(152:304))%>%
  group_by(region, Scenario, IBMR) %>%
  mutate(Biomass = na.locf(Biomass)) %>% #if we are missing data, carry last observation forward
  ungroup() %>%
  pivot_wider(id_cols = c(region, DOY,  Scenario),
              names_from = IBMR, values_from = Biomass) %>%
  select(region,  DOY,  Scenario, #put it in the right order for the model
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) %>%
  arrange(Scenario)

#create an array that is days x taxa x region x scenario
test24 = zoopsmAveSummer24 %>%
  split(list(zoopsmAveSummer24$region, zoopsmAveSummer24$Scenario))

zoop24 = array(unlist(test24),dim=c(153,15,6, 4), 
             dimnames = list(c(152:304), names(zoopsmAveSummer24), sort(unique(zoopsmAveSummer24$region)),
                             sort(unique(zoopsmAveSummer24$Scenario))))



#OK! Now i just need to get rid of the id columns
zoopx24 = zoop24[c(1:153), c(4:15), c(1:6), c(1:4)]
zoopx224 = replace_na(apply(zoopx24, c(2,3,4), as.numeric),0)




################################################################
#day lengths
daylight = daylength(38.15, c(130:310))$Daylength*60
LT.fx <- daylight/(daylength(38.15, 173)$Daylength*60) #daylenght divided by daylength at summer solstace

beta_hat <- read.table(file='data/beta_hat.txt') # MC filtering coefficients
#this is model uncertainty


### 3. Model bioenergetics 

run2024 = smelt_bioenergetics_summer(PD.mn.array = zoopx224, obs.temp.dat = Tempx2R24, obs.turb.dat = Turbx2R24, 
                                     start.L = rep(23, 15), ex.strata = ex.strata, beta_hat = beta_hat[1:200,], 
                                     ex.scenarios = ex.scenarios)

save(run2024, file = "outputs/sfharun_2024posthoc.RData")
#load("outputs/sfharun_2024posthoc.RData")

#put both datasets together and calculate the mean and standard deviation of the outputs
summary_2024 = run2024%>%
  mutate(Year = "2024") %>%
  group_by(Stratum, Scenario, Day, Year) %>%
  summarize(sdLength = sd(Length), sdWeight = sd(Weight),Length = mean(Length), Weight = mean(Weight), 
            energy = mean(energy), food = mean(Food), Limno = mean(Limno)) %>%
  group_by(Stratum, Scenario, Year) %>%
  mutate(Growth = (Weight -dplyr::lag(Weight))/Weight) %>%
  ungroup()


ggplot(summary_2024, aes(x = Day, y = Growth, color = Scenario)) + geom_line()+
  facet_wrap(Year~Stratum)

#OOOOHHHHH, they're eating more limnoithona...

ggplot(summary_2024, aes(x = Day, y = Limno, color = Scenario)) + geom_line()+
  facet_wrap(Year~Stratum)

ggplot(summary_2024, aes(x = Day, y = energy, color = Scenario)) + geom_line()+
  facet_wrap(Year~Stratum)

ggplot(filter(summary_2024, Stratum == "SE Suisun"), aes(x = Day, y = Growth, color = Scenario)) + 
  geom_line()+
  coord_cartesian(ylim = c(0.008, 0.015), xlim = c(75, 155))
#  facet_wrap(Year~Stratum)


ggplot(summary_2024, aes(x = Day, y = Length, color = Scenario)) + geom_line()+
  facet_wrap(~Stratum)

ggplot(summary_2024, aes(x = Day, y = Weight, color = Scenario)) + geom_line()+
  facet_wrap(~Stratum)

#Now I want to see what the maximum length and weight they got was.
#since they all started out the same size, this tells me how much they grew.
totagrowth24 = summary_2024 %>%
  group_by(Stratum, Scenario, Year) %>%
  summarize(Length = max(Length, na.rm =T), Weight = max(Weight, na.rm =T), 
            sdLength = max(sdLength), sdWeight = max(sdWeight))%>%
  mutate(Scenario = factor(Scenario, levels = c("NoAction", "X2Only", "GatesOnly", "Gates+X2")))

##plot of legnth at the end of the summer ##############
ggplot(totagrowth24, aes(x = Scenario, y = Length, group = Year, fill = Scenario)) + 
  geom_col( color = "grey")+
  scale_alpha_manual(values = c(0.6, 1))+
  facet_wrap(~Stratum)+ coord_cartesian(ylim = c(52, 70))+
  geom_errorbar(aes(ymin = Length - sdLength, ymax = Length+sdLength), position = "dodge")+
  theme(axis.text.x = element_text(angle = 90))

##plot of weight at the end of the summer ##############
ggplot(totagrowth24, aes(x = Scenario, y = Weight, group = Year, fill = Scenario)) + 
  geom_col(color = "grey")+
  scale_alpha_manual(values = c(0.6, 1))+
  facet_wrap(~Stratum)+ coord_cartesian(ylim = c(1, 3))+
  geom_errorbar(aes(ymin = Weight - sdWeight, ymax = Weight+sdWeight), position = "dodge")+
  theme(axis.text.x = element_text(angle = 90))


###########################################################
#look at growth by month
totagrowth_month24 = summary_2024 %>%
  mutate(Month = case_match(Day, c(0:30) ~ 6,
                            c(31:62) ~ 7,
                            c(63:94) ~ 8,
                            c(95:125) ~ 9,
                            c(126:157) ~ 10)) %>%
  group_by(Stratum, Scenario, Year, Month) %>%
  summarize(Length = max(Length, na.rm =T), Weight = max(Weight, na.rm =T), 
            sdLength = max(sdLength), sdWeight = max(sdWeight))%>%
  mutate(Scenario = factor(Scenario, levels = c("NoAction", "X2Only", "GatesOnly", "Gates+X2")),
         Laggrowth = dplyr::lag(Weight), NewGrowth = Weight-Laggrowth,
         NewGrowth = case_when(Month ==6 ~ Weight - 0.07,
                               TRUE ~ NewGrowth))


ggplot(totagrowth_month24, aes(x = Scenario, y = NewGrowth, color = Year)) +
  facet_grid(Month ~ Stratum, scales = "free_y")+
  geom_point(position = "dodge")+ ylab("Monthly weigth increase")+
  theme(axis.text.x = element_text(angle = 90))

ggplot(totagrowth_month24, aes(x = Scenario, y = Weight, color = Year)) +
  facet_grid(Month ~ Stratum, scales = "free_y")+
  geom_point(position = "dodge")+ ylab("Weight")+
  theme(axis.text.x = element_text(angle = 90))

###############################################################################################
#OK, now I need to add the salinities and calculate the growth rate in the 
#areas where salinity is <6

#combine both salinity files
salsummary =  dat2020%>%
  rename(Stratum = region, Salinity = mean) %>%
  mutate(Date = date(time), DOY = yday(Date))

#mass-specific growth rate in each region
growth_6_24 = summary_2024 %>%
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
test = growth_6_24 %>%
  group_by(Scenario, Year) %>%
  do(Weight = growing(.))

testx = growth_6_24 %>%
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
  coord_cartesian(xlim = c(120, 155), ylim = c(0.015, 0.021))

test2x = mutate(test2x, DOY )

ggplot(test2x, aes(x = Day, y = Weight, color = Scenario))+
  geom_line() + facet_wrap(~Year)+ylab("Weight (g)")+
  scale_x_continuous(breaks = c(0,31,62,90, 123), labels = c("jun", "jul", "aug", "sep", "oct"))


#calculate mean growth rate and maximum weigth if fish are restricted to
#salinities of <6

#this is what does into the consequence table.
meangrowth24 = group_by(test2x, Year, Scenario) %>%
  summarise(Growth = mean(growth), Growth2 = mean(growth2),total = max(Weight, na.rm =T),
            total2 = max(Weight2, na.rm =T))%>%
  mutate(Scenario = factor(Scenario, levels = c("NoAction", "X2Only", "GatesOnly", "Gates+X2")))

write.csv(meangrowth24, "outputs/meangrowth_2024actions.csv", row.names = F)

ggplot(meangrowth24, aes(x = Scenario, y = total, fill = Scenario, group = Year)) + 
  geom_col()+
  scale_alpha_manual(values = c(0.7, 1))+
  coord_cartesian(ylim = c(1.5,2.1))+
  ylab("total summer growth (g)")

ggplot(meangrowth24, aes(x= Scenario, y = total2, fill = Scenario, group = Year)) + 
  geom_col()+
  scale_alpha_manual(values = c(0.7, 1))+
  coord_cartesian(ylim = c(1.5,2.3))+
  ylab("total summer growth (g)")


##########################################################
#map of regions for ITP ammendment
library(sf)
library(ggspatial)
library(deltamapr)

Regions = mutate(R_DSIBM,
                  Region = case_match(SUBREGION, "NW Suisun" ~ "Grizzly Bay",
                                      c("NE Suisun", "SE Suisun", "SW Suisun") ~ "Suisun bay",
                                      c("Yolo Bypass", "Sacramento River") ~ "North Delta",
                                      c("Lower Sacramento River", "Confluence") ~ "Lower Sacramento",
                                      .default = SUBREGION))
Reglable = group_by(Regions, Region) %>%
  summarize(across(everything(), first))


ggplot()+
  geom_sf(data = WW_Delta )+
  geom_sf(data = Regions, aes(fill = Region), alpha = 0.5)+
  geom_sf_label(data = Reglable, aes(label = Region))+
  coord_sf(ylim = c(38, 38.3), xlim = c(-122.2, -121.7))+
  ylab(NULL)+xlab(NULL)+
  scale_fill_discrete(guide = "none") +
  annotation_north_arrow(location = "tl")+
  annotation_scale()+
  theme_bw()
  
