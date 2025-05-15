#grab cdec data for input

library(cder)
library(chillR)
library(tidyverse)
library(readxl)
library(sf)
library(abind)
library(zooper)
library(deltamapr)
library(discretewq)
library(wql)
library(here)

#27 is turbidity in NTU, 221 is turbidity in FNU, 25 is water temperature in farenhieght
# CDEC_cages = cdec_query(c("BDL", "RVB"), sensors = c(27, 221, 25), 
#                         start.date = ymd("2023-08-30"), end.date = ymd("2023-10-14"))
# save(CDEC_cages, file = "CDEC_cages.RData")
#load("CDEC_cages.RData")

WQ20172024 = read_csv("data/SMSCG_wq_data_2017-2024.csv") %>%
  rename(StationID = station) %>%
 mutate(Date = date(date_time_pst), Year = year(Date)) %>%
  select(Year, StationID, Date, watertemperature, turbidity, salinity) %>%
  pivot_longer(cols = c(watertemperature, turbidity, salinity), names_to = "Parameter", values_to = "Value")

cagedates = data.frame(Year = c(2019, 2023, 2024), 
                       StartDate = c(ymd("2019-10-9"), ymd("2023-08-30"), ymd("2024-08-27")),
                       EndDate = c(ymd("2019-11-6"), ymd("2023-10-12"), ymd("2024-09-26")))

CDEC_cages = left_join(WQ20172024, cagedates) %>%
  filter(Date>= StartDate & Date <= EndDate, StationID %in% c("BDL", "RVB"))

CDEC_cagesdaily = CDEC_cages %>%
  group_by(Date, StationID, Parameter) %>%
  summarize(Value = mean(Value, na.rm =T))

CDEC_wide = pivot_wider(CDEC_cagesdaily, id_cols = Date, names_from = c(Parameter, StationID), values_from = Value,
                        values_fn = mean)

#OK, now I need day lenght in minutes

#input data is julian day and latitude. Latitude 38.16 should work for both places
CDEC_wide = mutate(CDEC_wide, jday = yday(Date), daylight = daylength(38.15, jday)$Daylength*60) %>%
  rename(NTU.RVB = `turbidity_RVB`, NTU.BDL = `turbidity_BDL`,
         Temp.BDL = watertemperature_BDL,Temp.RVB = watertemperature_RVB)

save(CDEC_wide, file = c("data/CDEC_wide_cages.RData"))
##########################################################################
#now organize the zooplankton data


# zoops = read_excel("data/DWR_CageZoop2023_Complete_TEC_1.22.2024.xlsx",
#                    sheet = "Zoop Data") 

load("C:/Users/rhartman/OneDrive - California Department of Water Resources/smelt cages/Smelt-cages/data/allcagezoops.RData")
# I think I want to average the zooplankton from all samples over the previous two weeks. 
names(allzoops)
zoops = mutate(allzoops, uniqueID = paste(Date, Location, station),
               InOut = case_when(str_detect(station, "Outside") ~ "Outside",
                              str_detect(Location, "Inside") ~ "Inside"),
                                                CageNum = str_sub(station, 1, 6)) %>%
  mutate(Taxon = str_remove(Taxon, "\\."),
         Taxon = str_remove(Taxon, " spp"),
         Taxon = str_remove(Taxon, " sp"),
                            Taxon = str_remove(Taxon, " UNID"))

write.csv(sort(unique(zoops$Taxon)), file = "Cagezooptaxa.csv")

#data from 2019 was only outside, whereas 2023 and 2024 was insided and outside. Hm.

#I made up the biomass for tsome of the 'other' taxa
lookup = read_csv("Data/Zooptaxa.csv")

zoopIBMR = left_join(zoops, lookup, by = c("Taxon" = "Taxa")) %>%
  mutate(BPUE = CPUE*Biomass/1000) %>% #convert from #=ug per ten gallons to mg/m3
  group_by(IBMR, Location, Date, uniqueID) %>%
    summarize(BPUE = sum(BPUE, na.rm =T)) 

#add in zeros
zoopIBMRzeros = zoopIBMR %>%
  filter(!is.na(IBMR)) %>%
  select(BPUE, uniqueID, Location, Date, IBMR) %>%
  pivot_wider(names_from = IBMR, values_from=BPUE, values_fill = 0) %>%
  mutate(Year = year(Date)) %>%
  left_join(cagedates) %>%
  mutate(Location = case_when(Location %in% c("RV", "Rio Vista") ~ "RVB",
                              Location %in% c("SM", "Montezuma") ~ "BDL",
                              TRUE ~ Location)) %>%
  filter(Date>= StartDate & Date <= EndDate, Location %in% c("BDL", "RVB"))


#now put it in the right order and expand to all dates
library(zoo)
zoopIBMR2 = zoopIBMRzeros %>%
  select(Location, Date, limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, 
         allcopnaup, daphnia, othcyc, other, eurytem, pdiapfor)%>%
  mutate(Period = case_when(Date  >= ymd("2019-10-08") & Date <= ymd("2019-10-16") ~ "A",
                            Date  > ymd("2019-10-16") & Date <= ymd("2019-10-21") ~ "B",
                            Date  > ymd("2019-10-21") & Date <= ymd("2019-10-28") ~ "C",
                            Date  > ymd("2019-10-28") & Date <= ymd("2019-11-06") ~ "D",
 Date  >= ymd("2023-8-30") & Date <= ymd("2023-9-15") ~ "E",
                            Date  > ymd("2023-9-15") & Date <= ymd("2023-9-29") ~ "F",
 Date  > ymd("2023-9-29") & Date <= ymd("2023-10-11") ~ "G",
 Date >= ymd("2024-8-27") & Date <= ymd("2024-9-14") ~ "H",
 Date > ymd("2024-9-14") & Date <= ymd("2024-9-28") ~ "I",
                            TRUE ~ "J"))  %>%
  group_by(Location, Period) %>%
  summarize(across(limno:pdiapfor, mean))

dates = data.frame(Date = c(seq.Date(ymd("2019-10-08"), ymd("2019-11-06"), by =1),
                            seq.Date(as.Date("2023-08-30"), as.Date("2023-10-15"), by =1),
                            seq.Date(as.Date("2024-08-27"), as.Date("2024-09-27"), by =1)))
sites = data.frame(Location = c("RVB", "BDL"))



alldates = merge(dates, sites) %>%
                               mutate(Period = case_when(Date  >= ymd("2019-10-08") & Date <= ymd("2019-10-16") ~ "A",
                                                         Date  > ymd("2019-10-16") & Date <= ymd("2019-10-21") ~ "B",
                                                         Date  > ymd("2019-10-21") & Date <= ymd("2019-10-28") ~ "C",
                                                         Date  > ymd("2019-10-28") & Date <= ymd("2019-11-06") ~ "D",
                                                         Date  >= ymd("2023-8-30") & Date <= ymd("2023-9-15") ~ "E",
                                                         Date  > ymd("2023-9-15") & Date <= ymd("2023-9-29") ~ "F",
                                                         Date  > ymd("2023-9-29") & Date <= ymd("2023-10-11") ~ "G",
                                                         Date >= ymd("2024-8-27") & Date <= ymd("2024-9-14") ~ "H",
                                                         Date > ymd("2024-9-14") & Date <= ymd("2024-9-28") ~ "I",
                                                         TRUE ~ "J"))

zoopIBMR3 = full_join(zoopIBMR2, alldates, by = c("Period", "Location"))

zooptest = pivot_longer(zoopIBMR3, cols = c(limno:pdiapfor), names_to = "IBMR", values_to = "BPUE") %>%
  mutate(Year = year(Date))

#they don't eat nauplii or 'other'
ggplot(filter(zooptest, !IBMR %in% c("allcopnaup", "other")), aes(x = Date, y = BPUE, fill = IBMR)) +
  geom_area() +
  facet_wrap(Location ~ Year, scales = "free_x")+
  coord_cartesian(ylim = c(0,5))+
  ylab("Biomass per cubic meter")

ggplot(filter(zooptest, IBMR %in% c("pdiapfor", "pdiapjuv", "eurytem")), aes(x = Date, y = BPUE, fill = IBMR)) +
  geom_area() +
  facet_wrap(Location ~ Year, scales = "free_x")

#put other back in there
ggplot(filter(zooptest, !IBMR %in% c("allcopnaup")), aes(x = Date, y = BPUE, fill = IBMR)) +
  geom_area() +
  facet_wrap(Location ~ Year, scales = "free_x")+
  coord_cartesian(ylim = c(0,10))+
  ylab("Biomass per cubic meter")


#then I need this as an array, one for each strata, one row for each day
#and I'll do a different array for each year
montzoops19 = filter(zoopIBMR3, Location == "BDL", year(Date) == 2019) %>%
  ungroup() %>%
  select(-Location, -Date, -Period) %>%
  as.matrix()
  
Riozoops19 = filter(zoopIBMR3,Location == "RVB", year(Date) == 2019)%>%
                         ungroup() %>%
                         select(-Location, -Date, -Period) %>%
                         as.matrix()
#
zoop19 <- abind(montzoops19, Riozoops19, along=3)


montzoops23 = filter(zoopIBMR3, Location == "BDL", year(Date) == 2023) %>%
  ungroup() %>%
  select(-Location, -Date, -Period) %>%
  as.matrix()

Riozoops23 = filter(zoopIBMR3, Location == "RVB", year(Date) == 2023)%>%
                         ungroup() %>%
                         select(-Location, -Date, -Period) %>%
                         as.matrix()
#
zoop23 <- abind(montzoops23, Riozoops23, along=3)

montzoops24 = filter(zoopIBMR3, Location == "BDL", year(Date) == 2024) %>%
                          ungroup() %>%
                          select(-Location, -Date, -Period) %>%
                          as.matrix()
Riozoops24 = filter(zoopIBMR3, Location == "RVB", year(Date) == 2024)%>%
                         ungroup() %>%
                         select(-Location, -Date, -Period) %>%
                         as.matrix()
#
zoop24 <- abind(montzoops24, Riozoops24, along=3)

save(zoop19, zoop23, zoop24, file = "data/cagezoops_byyear.Rdata")


##############################################################################
#zoops from the long-term surveys instead.

fmwt2019 = Zoopsynther(Data_type = "Community", Sources = c("EMP", "FMWT", "STN"), Size_class = "Meso",
                       Years = c(2019:2024))
crosswalk = read_csv("data/zoopstaxa.csv") %>%
  select(Taxlifestage, IBMR, CarbonWeight_ug) %>%
  distinct()

fmwt24 = read_csv("C:/Users/rhartman/OneDrive - California Department of Water Resources/salinity control gates/SMSCG/Data/smscgto2024_zooplankton_long.csv")

fmwt2019_24 = bind_rows(fmwt2019, fmwt24) %>%
  mutate(Taxlifestage = str_replace_all(Taxlifestage, "Eurytemora affinis", "Eurytemora carolleeae"))

fmwt2019a = left_join(fmwt2019_24, cagedates) %>%
  filter(Date > StartDate-10 & Date < EndDate +10) %>%
  mutate(Week = week(Date)) %>%
  left_join(crosswalk) %>%
  mutate(Location = case_when(Station %in% c("606", "609", "611", "NZ032") ~ "BDL",
                              Station %in% c("706", "707", "708", "709", "711", "NZ068", "D22") ~ "RVB")) %>%
  filter(!is.na(Location))

#weekly regional averages

fmwtave = fmwt2019a %>%
  mutate(BPUE = case_when(is.na(BPUE) ~CPUE*CarbonWeight_ug,
                          TRUE ~ BPUE)) %>%
  group_by(SampleID, Location, Date, Week, Year, IBMR)  %>%
  summarize(BPUE = sum(BPUE, na.rm = T)) %>%
  group_by(Location, Week, Year, IBMR) %>%
  summarize(BPUE = mean(BPUE))

#sigh. Missing pseudo juv from 2024, lets see if we can fix tat.

fmwtaveb = fmwtave %>%
  pivot_wider(id_cols = c(Location, Week, Year), names_from = IBMR, values_from = BPUE, values_fill = 0) %>%
  pivot_longer(cols = c(acartela:`NA`), values_to = "BPUE", names_to = "IBMR")

unique(fmwtave$Week)
#we're missing a few weeks
library(imputeTS)
allweeks = data.frame(Year = c(2019,2019, 2019, 2019,2019,
                               rep(2023,8),
                               rep(2024, 5)),
                      Week = c(41,42,43,44,45,
                               35,36,37,38,39,40,41,42,
                               35,36,37,38,39)) %>%
  merge(data.frame(IBMR = unique(fmwtave$IBMR))) %>%
 merge(data.frame(Location = c("RVB", "BDL")))

fmwtavea = full_join(fmwtaveb, allweeks) %>%
  filter(!is.na(IBMR)) %>%
  arrange(Year, Location, IBMR, Week) %>%
  group_by(Year, Location, IBMR) %>% 
  mutate(BPUE = na_locf(BPUE, na_remaining = "mean"))
#expand to all dates

alldates2 = merge(dates, sites) %>%
  mutate(Week = week(Date), Year = year(Date))


fmwtweekly = left_join(fmwtavea, alldates2) %>%
  pivot_wider(names_from = "IBMR", values_from = "BPUE") %>%
  select(Location, Date, limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, 
         allcopnaup, daphnia, othcyc, other, eurytem, pdiapfor) %>%
  filter(!is.na(Date))


fmwtweekly_long = pivot_longer(fmwtweekly, cols = c(limno:pdiapfor), names_to = "IBMR", values_to = "BPUE") %>%
  mutate(Year = year(Date))

ggplot(fmwtweekly_long, aes(x = Date, y = BPUE, fill = IBMR))+
  geom_area()+
  facet_grid(Location~Year, scales= "free_x")


#then I need this as an array, one for each strata, one row for each day
#and I'll do a different array for each year
montzoops19fmwt = filter(fmwtweekly , Location == "BDL", year(Date) == 2019) %>%
  ungroup() %>%
  select(-Location, -Date, -Year) %>%
  as.matrix()

Riozoops19fmwt = filter(fmwtweekly,Location == "RVB", year(Date) == 2019)%>%
  ungroup() %>%
  select(-Location, -Date, -Year) %>%
  as.matrix()
#
zoop19fmwt <- abind(montzoops19fmwt, Riozoops19fmwt, along=3)


montzoops23fmwt = filter(fmwtweekly, Location == "BDL", year(Date) == 2023) %>%
  ungroup() %>%
  select(-Location, -Date, -Year) %>%
  as.matrix()

Riozoops23fmwt = filter(fmwtweekly, Location == "RVB", year(Date) == 2023)%>%
  ungroup() %>%
  select(-Location, -Date, -Year) %>%
  as.matrix()
#
zoop23fmwt <- abind(montzoops23fmwt, Riozoops23fmwt, along=3)

montzoops24fmwt = filter(fmwtweekly, Location == "BDL", year(Date) == 2024) %>%
  ungroup() %>%
  select(-Location, -Date, -Year) %>%
  as.matrix()
#
Riozoops24fmwt = filter(fmwtweekly, Location == "RVB", year(Date) == 2024)%>%
  ungroup() %>%
  select(-Location, -Date, -Year) %>%
  as.matrix()
#
#
zoop24fmwt <- abind(montzoops24fmwt, Riozoops24fmwt, along=3)

save(zoop19fmwt, zoop23fmwt, zoop24fmwt, file = "data/cagezoops_byyearfmwt.Rdata")


#####################dates##############################################################################
#Now let's try adjusting the zooplankton data so we have the same average across the entire 6 weeks per site
#this will allow us to see the impact of temperature better


zoopIBMRx = left_join(zoopavetaxa2, lookup) %>%
  mutate(BM = Count*Biomass, BPUE = (BM/1000)*(264.17/10)) %>% #convert from #=ug per ten gallons to mg/m3
  group_by(IBMR, Site) %>%
  summarize(BPUE = mean(BPUE, na.rm =T))

zoopIBMR2x = pivot_wider(zoopIBMRx, id_cols = c(Site), names_from = IBMR, values_from = BPUE,
                        values_fill = 0)
names(zoopIBMR2x)

#missing daphnia, eurytemora, 
#now put it in the right order and expand to all dates

zoopIBMR2x = mutate(zoopIBMR2x, daphnia = 0, eurytem =0) %>%
  select(Site, limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, 
         allcopnaup, daphnia, othcyc, other, eurytem, pdiapfor, amphipod)

alldatesx = data.frame(Date = rep(seq.Date(as.Date("2023-08-30"), as.Date("2023-10-15"), by =1), 2),
                      Site = c(rep("Rio Vista", 47), rep("Montezuma", 47))) 

zoopIBMR3x = full_join(zoopIBMR2x, alldatesx, by = c("Site"))

#then I need this as an array, one for each strata, one row for each day
montzoopsx = as.matrix(filter(zoopIBMR3x, Site == "Montezuma") %>%
                        ungroup() %>%
                        select(-Site, -Date))
Riozoopsx = as.matrix(filter(zoopIBMR3x, Site == "Rio Vista")%>%
                       ungroup() %>%
                       select(-Site, -Date))
#
zoopx <- abind(montzoopsx, Riozoopsx, along=3)

#now let's see what happens when zooplantkton at both sites are the same
zoopsx2 = abind(Riozoopsx, Riozoopsx, along=3)


########################################################
#data for water quality, 2010-2024

alldays = seq(ymd("2010-01-01"), ymd("2024-12-31"), by = "day")
Params = c("salinity", "turbidity", "watertemperature")
Regs = c("Lower Sacramento River",  "Confluence", "Suisun Marsh",
         "SW Suisun",  "SE Suisun", "NE Suisun", "NW Suisun")

allvals = data.frame(Date = alldays) %>%
  merge(data.frame(Parameter = Params)) %>%
  merge(data.frame(Region = Regs)) %>%
  mutate(DOY = yday(Date), Year = year(Date))


#IBMR regions
Regions = deltamapr::R_DSIBM %>%
  filter(SUBREGION %in% c("Lower Sacramento River",  "Confluence", "Suisun Marsh",
                          "SW Suisun",  "SE Suisun", "NE Suisun", "NW Suisun")) %>%
  rename(Region = SUBREGION)

write_sf(deltamapr::R_DSIBM, dsn= "data/IBMR_Regions.shp")


#i've got hte qc'd data for 2017-2024 from Morgan
WQ20172024 = read_csv("data/SMSCG_wq_data_2017-2024.csv") %>%
  rename(StationID = station, ObsDate = date_time_pst) %>%
  select(StationID, ObsDate, watertemperature, turbidity, salinity) %>%
  pivot_longer(cols = c(watertemperature, turbidity, salinity), names_to = "Parameter", values_to = "Value")

#grab the 2010-2016 data from CDEC
library(cder)
WQ20102016 = cdec_query(stations = unique(WQ20172024$StationID), sensors = c(100, #spc in microsemes
                                                                           25 #temperature in degrees F
                                                                           ),
                        start.date = ymd("2010-01-01"), end.date = ymd("2016-12-31"))
#Hm. It looks like we have no turbidity data for that entire time period, unless something went wrong. 

WQturb = cdec_query(stations = unique(WQ20172024$StationID), sensors = c(27, #Turbidity NTU
                                                                           221 #Turbidity FNU
),
start.date = ymd("2010-01-01"), end.date = ymd("2016-12-31"))
#OK, never mind, we just didn't get all the data!



#We don't have any data from the SW suisun region, so I"m going to pull MRZ
#we might as well grab PCT and ANH to improve our data in those regions too. 
extradata = cdec_query(c("PCT", "ANH", "MRZ"), sensors = c(100, #spc in microsemes
                                                           25, #temperature in degrees F
                                                           27, #Turbidity NTU
                                                           221), #Turbidity FNU),
                                                           start.date = ymd("2010-01-01"), end.date = ymd("2024-12-31"))

#something weird happewdn with this data
extradata1 = filter(extradata, !(SensorType == "EL COND"& ObsDate > ymd("2022-07-06") & 
                                   ObsDate < ymd("2022-07-21") & StationID == "PCT"),
                    !(SensorType == "EL COND"& ObsDate > ymd("2023-01-01") &
                        ObsDate < ymd("2023-01-6") & StationID == "PCT"))


EMM = cdec_query(c("EMM", "SDI"), sensors = c(100, #spc in microsemes
                                                           25, #temperature in degrees F
                                                           27, #Turbidity NTU
                                                           221), #Turbidity FNU),
                       start.date = ymd("2010-01-01"), end.date = ymd("2024-12-31"))
WQ20102016a = bind_rows(WQ20102016, WQturb, extradata1)
#filter out bad values, summarize by day
summary(WQ20102016a)


#yup, th'ats wrong

WQ20102016a = filter(WQ20102016a, Value >0, !is.na(Value),Value <900000, !(SensorUnits == "DEG F" & Value >82),
                     !(SensorUnits == "DEG F" & Value < 35),
                     !(SensorType == "EL COND" & Value > 55000),
                     !(SensorUnits == "NTU" & Value >500)) %>%
  mutate(Value2 = case_when(SensorUnits == "DEG F" ~ (Value-32)*5/9,
                            SensorUnits == "uS/cm" ~ ec2pss(Value/1000, 25),
                            TRUE ~ Value),
         Parameter = recode(SensorType, "EL COND" = "salinity", 
                            "TEMP W" = "watertemperature", "TURB W" = "turbidity")) %>%
  select(StationID, ObsDate, Parameter, Value2) %>%
  rename(Value = Value2)

ggplot(WQ20102016a, aes(x = ObsDate, y = Value, color = StationID))+geom_line()+
  facet_wrap(~Parameter, scales = "free_y", nrow =3)

#now I just need to add the regions. Sigh. 
cdecstations = read_csv("data/station_data.csv") %>%
  filter(!is.na(Latitude) )%>%
  select(station, Latitude, Longitude) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = F) %>%
  st_transform(crs = st_crs(Regions)) %>%
  st_join(Regions) %>%
  st_drop_geometry() %>%
  select(station, Latitude, Longitude, Region) %>%
  mutate(Region = case_when(station == "RVB" ~ "Lower Sacramento River",
                            TRUE ~ Region), StationID = station) %>%
  select(-station)

########map##########
#map of stations for synthesis report
stations_all <- read_csv("data/station_data.csv")%>%
  filter(!is.na(Latitude) )%>%
  select(station, Latitude, Longitude) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = F)
smeltregions = filter(R_DSIBM, SUBREGION %in% c("Confluence","NE Suisun","NW Suisun","SE Suisun","Suisun Marsh","SW Suisun","Lower Sacramento River"))

ggplot() +
  geom_sf(data = WW_Delta)+
  geom_sf(data = smeltregions, aes(fill = SUBREGION), alpha = 0.3)+
  geom_sf(data =stations_all)+
  coord_sf(ylim = c(38, 38.25), xlim = c(-122.15, -121.7))


#######
#put all the continuous water quality data together############
Allcontwq = bind_rows(WQ20102016a, WQ20172024)  %>%
  mutate(Date = date(ObsDate), Month = month(ObsDate), DOY = yday(ObsDate), Year = year(ObsDate))
Allcontwqmean = Allcontwq %>%
  group_by(Date, Parameter, StationID, Month, DOY, Year) %>%
  summarize(Value = mean(Value, na.rm =T)) %>%
  left_join(cdecstations) %>%
  mutate(Source = "CDEC")

save(Allcontwq, Allcontwqmean, cdecstations, file = "ContinuousWQ.RData")
load("ContinuousWQ.RData")

#how much data by year and station?

WQdatasummary = group_by(Allcontwqmean, StationID, Year, Parameter) %>%
  summarize(n = n()) 

ggplot(WQdatasummary, aes(x = Year, y = n, fill = StationID)) + geom_col()+
  facet_wrap(~Parameter)

#we only have turbidity from a few stations, may want to fill in with turbidity from discrete measurements
library(discretewq)
WQdiscrete = wq(Sources = c("20mm", "Baystudy", "DJFMP", "DOP", "EDSM", 
                            "EMP", "FMWT", "STN", "USGS_CAWSC", "USGS_SFBS"), Start_year = 2010, End_year = 2022)

WQdiscretelong = WQdiscrete %>%
  select(Source, Station, Latitude, Longitude, Date, Secchi, Temperature, Conductivity, Salinity, TurbidityFNU, TurbidityNTU) %>%
  mutate(Turbidity = case_when(is.na(TurbidityFNU) ~ (1.1*TurbidityNTU-7)/0.74,
         TRUE ~ TurbidityFNU)) %>%
  pivot_longer(cols = c(Turbidity, Salinity, Temperature), names_to = "Parameter", values_to = "Value")

ggplot(WQdiscretelong, aes(x = Date, y = Value, color = Source)) +geom_point()+
  facet_wrap(~Parameter, scales = "free_y", nrow =3)
#oh, there are some wacky salinity values, need to follow up with Sam. 

saltest = filter(WQdiscrete, Salinity >50)
#huh, it looks like issues with the underlying dataset. That's a problem for another day thoguh 

#How can I create regional means that are fair?
WQdiscretelongf = filter(WQdiscretelong, !(Parameter == "Salinity" & Value >40),
                         !(Parameter == "Temperature" & Value >35), !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = F) %>%
  st_transform(crs = st_crs(Regions)) %>%
  st_join(Regions) %>%
  st_drop_geometry() %>%
  filter(!is.na(Region)) %>%
  rename(StationID = Station) %>%
  select(Source, StationID, Latitude, Longitude, Date, Parameter, Value, Region) %>%
  mutate(Year = year(Date), Month = month(Date), DOY = yday(Date),
  Parameter = recode(Parameter, "Salinity" = "salinity", 
         "Temperature" = "watertemperature", "Turbidity" = "turbidity"))

names(WQdiscretelongf)
names(Allcontwqmean)

AllWQ = bind_rows(WQdiscretelongf, Allcontwqmean) %>%
  filter(!(Parameter == "salinity" & Value >35))

ggplot(AllWQ, aes(x = Date, y = Value, color = Region))+
  facet_grid(Parameter~Region, scales = "free_y")+ geom_point()

#I guess I'll need to impute missing values?
AllWQmean = group_by(AllWQ, Date, Year, Month, DOY, Parameter, Region) %>%
  summarize(Value = mean(Value, na.rm =T))

ggplot(AllWQmean, aes(x = Date, y = Value, color = Region))+
  facet_grid(Parameter~Region, scales = "free_y")+ geom_point()

#expanded data frame with NA's
AllWQmean2 = left_join(allvals, AllWQmean)

#impute na's
library(imputeTS)

AllWQmean2 = AllWQmean2 %>%
  group_by(Region, Parameter) %>%
  arrange(Date) %>%
  mutate(ValueImputed = na_interpolation(Value, option = "linear"), DOY = yday(Date),
         daylight = daylength(38.15, DOY)$Daylength*60,)

ggplot(AllWQmean2, aes(x = Date, y = ValueImputed, color = Region)) + geom_line()+
  facet_grid(Parameter~Region, scales = "free_y")

foowe = filter(AllWQmean2, Year == 2015, Region == "NE Suisun", DOY >100)

ggplot(AllWQmean2, aes(x = DOY, y = ValueImputed, color = Region)) + geom_line()+
  facet_grid(Parameter~Year, scales = "free_y")
#That's better!. Still might have some issues with that one high tubidity value, but let's see how that goes. 

#I think the model wants it in wdie format?
#shoot, I need to turn them into arrays by year
Turbwide = filter(AllWQmean2, Parameter == "turbidity") %>%
  select(Region, Date, Parameter, Year, DOY, ValueImputed) %>%
  distinct() %>%
  arrange(Region) %>%
  pivot_wider(names_from = Region, values_from = ValueImputed)


Turbwidef = filter(Turbwide, DOY %in% c(153:305))

test = Turbwidef  %>%
  split(list(Turbwidef$Year))

Turb = array(unlist(test),dim=c(153,11, 15), 
             dimnames = list(c(153:305), names(Turbwidef),
                             unique(Turbwidef$Year)))

#OK! Now i just need to get rid of the id columns
Turbx = Turb[c(1:153), c(5:11), c(1:15)]
Turbx2 = apply(Turbx, c(2,3), as.numeric)


Tempwide = filter(AllWQmean2, Parameter == "watertemperature")%>%
  select(Region, Date, Parameter, Year, DOY, ValueImputed) %>%
  distinct() %>%
  arrange(Region) %>%
  pivot_wider(names_from = Region, values_from = ValueImputed)


Tempwidef = filter(Tempwide, DOY %in% c(153:305))

test = Tempwidef  %>%
  split(list(Tempwidef$Year))

Temp = array(unlist(test),dim=c(153,11, 15), 
             dimnames = list(c(153:305), names(Tempwidef),
                             unique(Tempwidef$Year)))

#OK! Now i just need to get rid of the id columns
Tempx = Temp[c(1:153), c(5:11), c(1:15)]
Tempx2 = apply(Tempx, c(2,3), as.numeric)

# I actually might not need salinity
Salwide = filter(AllWQmean2, Parameter == "salinity") %>%
  select(-Value) %>%
  pivot_wider(names_from = Region, values_from = ValueImputed)

save(AllWQmean2, Tempx2, Turbx2, file = "WaterQuality20102024.RData")

#################################################
#constant temperature
# calculate average temperature 
Tempave = mean(filter(AllWQmean2, Parameter == "watertemperature", DOY %in% c(153:305))$ValueImputed)

Tempwide_constant = filter(AllWQmean2, Parameter == "watertemperature")%>%
  select(Region, Date, Parameter, Year, DOY, ValueImputed) %>%
  distinct() %>%
  mutate(ValueImputed = Tempave) %>%
  arrange(Region) %>%
  pivot_wider(names_from = Region, values_from = ValueImputed)


Tempwidef_constant = filter(Tempwide_constant, DOY %in% c(153:305))

test_constant = Tempwidef_constant  %>%
  split(list(Tempwidef_constant$Year))

Temp_constant = array(unlist(test_constant),dim=c(153,11, 15), 
             dimnames = list(c(153:305), names(Tempwidef_constant),
                             unique(Tempwidef_constant$Year)))

#OK! Now i just need to get rid of the id columns
Tempx_constant = Temp_constant[c(1:153), c(5:11), c(1:15)]
Tempx2_constant = apply(Tempx_constant, c(2,3), as.numeric)

#################################################
#constant turbidity
# calculate average turbidity 
turbave = mean(filter(AllWQmean2, Parameter == "turbidity", DOY %in% c(153:305))$ValueImputed)

turbwide_constant = filter(AllWQmean2, Parameter == "turbidity")%>%
  select(Region, Date, Parameter, Year, DOY, ValueImputed) %>%
  distinct() %>%
  mutate(ValueImputed = turbave) %>%
  arrange(Region) %>%
  pivot_wider(names_from = Region, values_from = ValueImputed)


turbwidef_constant = filter(turbwide_constant, DOY %in% c(153:305))

test_constantv = turbwidef_constant  %>%
  split(list(turbwidef_constant$Year))

turb_constant = array(unlist(test_constantv),dim=c(153,11, 15), 
                      dimnames = list(c(153:305), names(turbwidef_constant),
                                      unique(turbwidef_constant$Year)))

#OK! Now i just need to get rid of the id columns
turbx_constant = turb_constant[c(1:153), c(5:11), c(1:15)]
turbx2_constant = apply(turbx_constant, c(2,3), as.numeric)

# I actually might not need salinity
Salwide = filter(AllWQmean2, Parameter == "salinity") %>%
  select(-Value) %>%
  pivot_wider(names_from = Region, values_from = ValueImputed)

save(AllWQmean2, Tempx2, Turbx2, turbx2_constant, Tempx2_constant, file = "WaterQuality20102024.RData")


##########################################################################################
#just mesozoops to put in BEM


#Why don't i have DOP on here?
# zoop_data_meso<-Zoopsynther(Data_type="Community", Sources=c("EMP", "STN", "20mm", "FMWT", "DOP"), 
#                             Time_consistency = FALSE, Years = c(2010:2022), Size_class = "Meso")

zoop_taxa<-read_csv(here("data/zoopstaxa.csv"), col_types=cols_only(Taxlifestage="c", IBMR="c", CarbonWeight_ug = "d"))%>%
  distinct()

#IBMR regions
Regions = deltamapr::R_DSIBM %>%
  filter(SUBREGION %in% c("Lower Sacramento River",  "Confluence", "Suisun Marsh",
                          "SW Suisun",  "SE Suisun", "NE Suisun", "NW Suisun")) %>%
  rename(Region = SUBREGION)

# zoopsallm = zoop_data_meso %>%
#   filter(!is.na(Longitude), !is.na(SalSurf)) %>%
#   mutate(Taxlifestage=recode(Taxlifestage, `Synchaeta bicornis Adult`="Synchaeta Adult", # Change some names to match to biomass conversion dataset
#                              `Pseudodiaptomus Adult`="Pseudodiaptomus forbesi Adult",
#                              `Acanthocyclops vernalis Adult`="Acanthocyclops_UnID Adult"))%>%
#   left_join(zoop_taxa, by="Taxlifestage")%>% # Add IBMR categories
#   mutate(BPUE = CPUE * CarbonWeight_ug/1000) %>%
#   select(SampleID, Station, Latitude, Longitude, SalSurf, Date, Year, IBMR, CPUE, BPUE)%>%
#   group_by(SampleID, Station, Latitude, Longitude, SalSurf, Date, Year, IBMR)%>%
#   summarise(CPUE=sum(CPUE, na.rm =T), BPUE = sum(BPUE, na.rm =T), .groups="drop")%>% # Sum each IBMR categories
#   st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
#   st_transform(crs=st_crs(Regions)) %>% 
#   st_join(Regions %>%
#             select(Region)) %>%
#   st_drop_geometry() %>% 
#   filter(!is.na(Region))%>%
#   mutate(doy=yday(Date), #Day of year
#          Month=month(Date), # Month
#          Year_fac=factor(Year), # Factor year for model random effect
#          Station_fac=factor(Station), # Factor station for model random effect
#          across(c(SalSurf, doy), list(s=~(.x-mean(.x))/sd(.x))), # Center and standardize predictors
#          CPUE_log1p=log(CPUE+1),
#          BPUE_log1p = log(BPUE+1)) # log1p transform BPUE for model, then convert

#put it in the right format for the model
#order of taxa is: limno, othcaljuv, pdiapjuv, othcalad, acartela, othflad, allcopnaup, daphnia, othcyc, other, eurytem, pdipfor


load("data/sfhazoops.RData")
zoopsallm = sfhazoops

daystaxa = merge(data.frame(Date = alldays), data.frame(IBMR = unique(zoopsallm$IBMR))) %>%
  merge(data.frame(Region = unique(zoopsallm$Region))) %>%
  mutate(Year = year(Date), doy = yday(Date), Month = month(Date))
# 
# zoopsmwide = zoopsallm %>%
#    group_by(Region, IBMR, Date, Month, Year, doy) %>%
#   summarize(BPUE = mean(BPUE)) %>%
#   filter(!is.na(IBMR)) %>%
#   full_join(daystaxa) %>%
#   group_by(Region, IBMR, Year) %>%
#   arrange(doy) %>%
#   mutate(BPUEImputed = na_locf(ts(BPUE))) %>% #repeat most recent value
#   ungroup() %>%
#   pivot_wider(id_cols = c(Region,  Date, Month, doy,  Year),
#                          names_from = IBMR, values_from = BPUEImputed) %>%
#   select(Region,  Date, Month, doy,  Year,
#          limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
#          daphnia, othcyc, other, eurytem, pdiapfor) %>%
#   mutate(Day = mday(Date)) 

#try doing monthly averages instead of LOCF
#shoudl this be geometric mean or median insted of mean?

zoopsmAvea = zoopsallm %>%
  group_by(Region, IBMR, Month, Year) %>%
  summarize(BPUE = mean(BPUE, na.rm =T)) %>%
  filter(!is.na(IBMR)) %>%
  full_join(daystaxa) %>%
  ungroup() %>%
mutate(BPUE = case_when(is.na(BPUE) & Year == 2024 & Region == "SW Suisun" ~ mean(filter(zoopsallm, Region == "NW Suisun", Year == 2024, Month %in% c(6:10))$BPUE),
                        is.na(BPUE) & Year == 2015 & Region == "NE Suisun" ~ mean(filter(zoopsallm, Region == "NW Suisun", Year == 2015, Month == 9)$BPUE),
                        is.na(BPUE) & Year == 2019 & Region == "NE Suisun" ~ mean(filter(zoopsallm, Region == "NW Suisun", Year == 2019, Month == 10)$BPUE),
                         TRUE ~ BPUE)) 

zoopsmAve = zoopsmAvea%>%
  pivot_wider(id_cols = c(Region,  Date, Month, doy,  Year),
              names_from = IBMR, values_from = BPUE) %>%
  select(Region,  Date, Month, doy,  Year,
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) %>%
  mutate(Day = mday(Date)) 

#we are randomly missingm September of 2015 and daphnia in october of 2019. Also a lot of 2024. 

#pull this out
# pseduo = zoopsallm %>%
#   group_by(Region, IBMR, Date, Month, Year, doy) %>%
#   summarize(BPUE = mean(BPUE)) %>%
#   full_join(daystaxa) %>%
#   filter(IBMR == "pdiapfor", Year == 2017, Region == "Suisun Marsh") %>%
#   ungroup() %>%
#   arrange(doy) %>%
#   mutate(BPUEImputed = na_locf(ts(BPUE)))


ggplot(filter(zoopsmAve, doy %in% c(200:300), Year == 2024),
       aes(x = Date, y = pdiapfor)) + geom_line()+ geom_point()+
  facet_wrap(~Region)

#I don't have any data from SW suisun ffor 2024. I'm going to copy the NW suisun data. 

zoopsmwidef = filter(zoopsmAve, doy %in% c(153:305))%>%
  arrange(Region)

summary(zoopsmwidef)
arg = filter(zoopsmwidef, is.na(pdiapfor))
anyNA(zoopsmwidef)
#also missing June of 2024. NOt sure we can do much about that.

summaryzoop = group_by(zoopsmwidef, Year, Region) %>%
  summarize(n())

test = zoopsmwidef %>%
    split(list(zoopsmwidef$Region, zoopsmwidef$Year))

zoop = array(unlist(test),dim=c(153,18,7, 15), 
             dimnames = list(c(153:305), names(zoopsmwidef), unique(zoopsmwidef$Region),
                                                              unique(zoopsmwidef$Year)))
  

#OK! Now i just need to get rid of the id columns
zoopx = zoop[c(1:153), c(6:17), c(1:7), c(1:15)]
zoopx2 = apply(zoopx, c(2,3,4), as.numeric)
PD.mn.array = zoopx2 #days by prey by strata by year!



save(zoopsmwide, zoopsmAve, zoopx2, file = "zoopsmwide.RData")

#
###########constant zooplankton ##############


zoops_constant = zoopsallm %>%
  group_by(IBMR) %>%
  summarize(BPUE = mean(BPUE, na.rm =T)) %>%
  filter(!is.na(IBMR)) %>%
  full_join(daystaxa) %>%
  group_by(Region, IBMR, Year) %>%
  arrange(doy) %>%
  pivot_wider(id_cols = c(Region,  Date, Month, doy,  Year),
              names_from = IBMR, values_from = BPUE) %>%
  select(Region,  Date, Month, doy,  Year,
         limno, othcaljuv, pdiapjuv, othcalad, acartela, othclad, allcopnaup, 
         daphnia, othcyc, other, eurytem, pdiapfor) %>%
  mutate(Day = mday(Date)) 


zoopsmwidef_constant = filter(zoops_constant, doy %in% c(153:305))%>%
  arrange(Region)


test_constant = zoopsmwidef_constant %>%
  split(list(zoopsmwidef_constant$Region, zoopsmwidef_constant$Year))

zoop_constant = array(unlist(test_constant),dim=c(153,18,7, 15), 
             dimnames = list(c(153:305), names(zoopsmwidef_constant), unique(zoopsmwidef_constant$Region),
                             unique(zoopsmwidef_constant$Year)))


#OK! Now i just need to get rid of the id columns
zoopx_constant = zoop_constant[c(1:153), c(6:17), c(1:7), c(1:15)]
zoopx2_constant = apply(zoopx_constant, c(2,3,4), as.numeric)
PD.mn.array_constant = zoopx2_constant #days by prey by strata by year!



save(zoops_constant,zoopsmwidef_constant , zoopx2_constant, file = "zoopsmwide_constant.RData")

########################################
#waer qulity summaries

#kick out growth values that came from bad salinity
sal = filter(AllWQmean2, Parameter == "salinity") %>%
  rename(Day = DOY, salinity = ValueImputed) %>%
  mutate(GoodSalinity = case_when(salinity <6 ~ T,
                                  salinity >= 6 ~ F),
         GoodSalinity2 = case_when(salinity <2 ~ T,
                                   salinity >= 2 ~ F),
         Region2 = str_remove_all(Region, " "))

ggplot(sal,  aes(x = Day, y = salinity, color = Region2)) +
  geom_point()+
  geom_smooth()+
  facet_wrap(~Year)

save(sal, file = "data/sal.Rdata")

#What was the average temperature in regions with good salinity?

Meantemps6 = AllWQmean2 %>%
  pivot_wider(id_cols = c(DOY, Date, Region, Year, Month), names_from = Parameter,
              values_from = ValueImputed) %>%
  filter(salinity <6) %>%
  group_by(Year, Month, DOY) %>%
  summarize(Temperature = mean(watertemperature)) %>%
  filter(Month %in% c(6:9)) %>%
  group_by(Year) %>%
  summarize(Temperature = mean(Temperature, na.rm = T))


Meantemps2 = AllWQmean2 %>%
  pivot_wider(id_cols = c(DOY, Date, Region, Year, Month), names_from = Parameter,
              values_from = ValueImputed) %>%
  filter(salinity <2) %>%
  group_by(Year, Month, DOY) %>%
  summarize(Temperature = mean(watertemperature)) %>%
  filter(Month %in% c(6:9)) %>%
  group_by(Year) %>%
  summarize(Temperature = mean(Temperature, na.rm = T))


Meantemps_nosalt = AllWQmean2 %>%
  pivot_wider(id_cols = c(DOY, Date, Region, Year, Month), names_from = Parameter,
              values_from = ValueImputed) %>%
  group_by(Year, Month, DOY) %>%
  summarize(Temperature = mean(watertemperature)) %>%
  filter(Month %in% c(6:9)) %>%
  group_by(Year) %>%
  summarize(Temperature = mean(Temperature, na.rm = T))



#add average X2 
load("data/DF.Rdata")

DF = Dayflow %>%
  filter(Month %in% c(6,7,8,9)) %>%
  group_by(Year) %>%
  summarize(X2 = mean(X2, na.rm =T))

save(DF, Meantemps_nosalt, Meantemps2, Meantemps6, file = "data/DF.RData")

#################################zoops time series ####################
library(RColorBrewer)
mypal = c(brewer.pal(8, "Dark2"), "grey", "firebrick", "blue", "darkgreen", "pink", "purple", "tan")
ggplot( zoopsallm, aes(x = Date, y = BPUE_log1p, fill = IBMR))+ geom_area()+
  facet_wrap(~Region)+
  scale_fill_manual(values = mypal)+
  coord_cartesian(ylim = c(0, 10))
# I b et tortanus is throwing it off

zoopsm = zoopsallm %>%
  group_by(Region, IBMR, Month, Year) %>%
  summarize(BPUE = mean(BPUE)) %>%
  mutate(Monthyear = Year + (Month-1)/12) %>%
  filter(!is.na(IBMR))

ggplot( zoopsm, aes(x = Monthyear, y = BPUE, fill = IBMR))+ geom_area()+
  facet_wrap(~Region)+
  scale_fill_manual(values = mypal)+
  coord_cartesian(ylim = c(0, 30))

ggplot( zoopsallm, aes(x = Date, y = BPUE, fill = IBMR))+ geom_area()+
  facet_wrap(~Region)+
  scale_fill_manual(values = mypal)+
  coord_cartesian(ylim = c(0, 50))

##########################################################
#pull out tortanus ####


zoopsallm_tort = zoop_data_meso %>%
  filter(!is.na(Longitude), !is.na(SalSurf)) %>%
  mutate(Taxlifestage=recode(Taxlifestage, `Synchaeta bicornis Adult`="Synchaeta Adult", # Change some names to match to biomass conversion dataset
                             `Pseudodiaptomus Adult`="Pseudodiaptomus forbesi Adult",
                             `Acanthocyclops vernalis Adult`="Acanthocyclops_UnID Adult"))%>%
  left_join(zoop_taxa, by="Taxlifestage")%>% # Add IBMR categories
  mutate(IBMR = case_when(Taxlifestage %in% c("Tortanus_UnID Adult", "Tortanus_UnID Juvenile") ~ "Tortanus",
                          TRUE ~ IBMR)) %>%
  mutate(BPUE = CPUE * CarbonWeight_ug/1000) %>%
  select(SampleID, Station, Latitude, Longitude, SalSurf, Date, Year, IBMR, CPUE, BPUE)%>%
  group_by(SampleID, Station, Latitude, Longitude, SalSurf, Date, Year, IBMR)%>%
  summarise(CPUE=sum(CPUE, na.rm =T), BPUE = sum(BPUE, na.rm =T), .groups="drop")%>% # Sum each IBMR categories
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=st_crs(Regions)) %>% 
  st_join(Regions %>%
            select(Region)) %>%
  st_drop_geometry() %>% 
  filter(!is.na(Region))%>%
  mutate(doy=yday(Date), #Day of year
         Month=month(Date), # Month
         Year_fac=factor(Year), # Factor year for model random effect
         Station_fac=factor(Station), # Factor station for model random effect
         across(c(SalSurf, doy), list(s=~(.x-mean(.x))/sd(.x))), # Center and standardize predictors
         CPUE_log1p=log(CPUE+1),
         BPUE_log1p = log(BPUE+1)) 

ggplot( zoopsallm_tort, aes(x = Date, y = BPUE_log1p, fill = IBMR))+ geom_area()+
  facet_wrap(~Region)+
  scale_fill_manual(values = mypal)+
  coord_cartesian(ylim = c(0, 10))
# I b et tortanus is throwing it off


ggplot( zoopsallm_tort, aes(x = Date, y = BPUE, fill = IBMR))+ geom_area()+
  facet_wrap(~Region)+
  scale_fill_manual(values = mypal)+
  coord_cartesian(ylim = c(0, 50))


ggplot( zoopsallm_tort, aes(x = Date, y = BPUE, fill = IBMR))+ geom_area(position = "fill")+
  facet_wrap(~Region)+
  scale_fill_manual(values = mypal)+
  coord_cartesian(ylim = c(0, 1))




ggplot( filter(zoopsallm_tort, Month %in% c(6:10)), aes(x = Year, y = BPUE, fill = IBMR))+ 
  geom_col(position = "fill")+
  facet_grid(Month~Region)+
  scale_fill_manual(values = mypal)+
  coord_cartesian(ylim = c(0, 1))

zoopstortsum =  filter(zoopsallm_tort, Month %in% c(6:10)) %>%
  group_by(IBMR, Month, Year, Region) %>%
  summarize(CPUE = mean(CPUE, na.rm =T), BPUE = mean(BPUE, na.rm = T))

ggplot(zoopstortsum,aes(x = Year, y = BPUE, fill = IBMR))+ 
  geom_col()+
  facet_grid(Month~Region)+
  scale_fill_manual(values = mypal)


ggplot( filter(zoopsallm_tort, IBMR == "limno"), aes(x = Date, y = BPUE, fill = IBMR))+ geom_line()+
  facet_wrap(~Region)+
  scale_fill_manual(values = mypal)+
  coord_cartesian(ylim = c(0, 5))
 
#Alll the zoops put togehter
yrs = read_csv("data/YearTypes_actions.csv")  %>%
  select(Year2, `Yr-type`) %>%
  distinct()

zoopsallave = group_by(zoopstortsum) %>%
  group_by( Month, Year, Region) %>%
  summarize(CPUE = sum(CPUE, na.rm =T), BPUE = sum(BPUE, na.rm = T))  %>%
  left_join(yrs, by = c("Year" = "Year2")) %>%
  filter(!is.na(`Yr-type`))

ggplot(zoopsallave, aes(x = `Yr-type`, y = BPUE)) +
  geom_boxplot()+
  facet_grid(Month~Region)+
  scale_y_log10()



#let's just do the monthly average biomass of each critter, monthly average temp, turb. 


zoopstortsum2 = mutate(zoopstortsum, Parameter = "Zooplankton", Value = BPUE) %>%
  filter(Month %in% c(6:10)) %>%
  group_by(Year, Parameter, Region, IBMR) %>%
  summarize(Value = mean(BPUE, na.rm =T))

Allwqmonth = AllWQmean %>%
  filter(Month %in% c(6:10)) %>%
  group_by(Year, Parameter, Region) %>%
  summarize(Value = mean(Value, na.rm = T))

allstuff = bind_rows(zoopstortsum2, Allwqmonth) %>%
  left_join(yrs, by = c("Year" = "Year2")) %>%
  left_join(DF)

ggplot(allstuff, aes(x = Year, y = Value, fill = IBMR)) +
  geom_col()+
  facet_grid(Parameter~Region, scales = "free")+
  scale_fill_manual(values = mypal)



ggplot(filter(allstuff, Parameter != "Zooplankton"), aes(x = Year, y = Value, color = `Yr-type`)) +
  geom_point()+
  facet_grid(Parameter~Region, scales = "free")+
  scale_color_manual(values = mypal)


ggplot(filter(allstuff, Parameter != "Zooplankton"), aes(x = X2, y = Value)) +
  geom_point(aes(color = `Yr-type`))+
  facet_grid(Parameter~Region, scales = "free")+
  scale_color_manual(values = mypal) +
  geom_smooth(method = "lm")

ggplot(filter(allstuff, Parameter == "Zooplankton", !is.na(IBMR)), aes(x = X2, y = Value)) +
  geom_point(aes(color = `Yr-type`))+
  facet_grid(IBMR~Region, scales = "free")+
  scale_color_manual(values = mypal) +
  geom_smooth(method = "lm")
