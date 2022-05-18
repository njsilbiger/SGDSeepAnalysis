# playing around

library(tidyverse)
library(seacarb)
library(broom)
library(lubridate)

# This is temporary until we get ther rest of the data
Cdata_orig<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/CarbonateChemistry/pHProbe_Data_calculated_NOTPOcorrect.csv") %>%
  mutate(datetime = mdy_hms(paste(as.character(Date), as.character(SamplingTime))))

locations<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/Sandwich_Locations_Final.csv")

# only run cases where both pH and TA are present of the code does not work
Cdata<- Cdata_orig %>%
  drop_na(TA,pH)

# calculate carbonate parameters
CO2<-carb(flag=8, Cdata$pH, Cdata$TA/1000000, S=Cdata$Salinity, 
          T=Cdata$TempInSitu, Patm=1, P=0, Pt=0, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")

#TA is divided by 1000000 because all calculations are in mol/kg in the seacarb package


#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
CO2all<-CO2 %>%
  mutate_at(vars(c("CO2","HCO3","CO3","DIC","ALK")), .funs = list(~.*1000000)) %>% # multiple everything by 1e6
  select(!c(pH, ALK, S, T,Patm, flag, P) )%>% #remove the columns we don't need
  bind_cols(Cdata,.) # bind it with the Cdata on the RHS 

# add back in all the data with missing TA and pH values
totaldata<-left_join(Cdata_orig, CO2all)


totaldata %>%
  left_join(locations) %>%
  ggplot(aes(x = TA*Salinity/36, y = DIC*Salinity/36, color = paste(Day_Night, Tide)))+
  geom_point()+
  facet_wrap(Location~Plate_Seep, scale = "free")

models<-totaldata %>%
  left_join(locations) %>%
  filter(Plate_Seep == "Plate") %>%
  mutate(TA_sal = TA*Salinity/36, # salinity normalize
         DIC_sal = DIC*Salinity/36) %>%
  nest(data = -c(Location,CowTagID)) %>%
  mutate(fit = map(data, ~lm(TA_sal~DIC_sal, data = .)),
         coefs = map(fit, tidy)) %>%
  select(!fit)%>%
  unnest(cols = coefs) %>%
  filter(term == "DIC_sal") %>%
  unnest(cols = data) %>%
  group_by(Location,CowTagID)%>%
  summarise_at(vars(pH, Salinity, estimate), .funs = ~min(.x,na.rm=TRUE))


models %>%
  ggplot(aes(x = Salinity, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location, scales = "free")

  
totaldata %>%
  left_join(locations) %>%
ggplot(aes(x = TA, y = DIC))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*CowTagID, scale = "free")

##### bring in the august and the march data together

AugData<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochem_wCarb.csv") %>%
  mutate(Season = "Dry") %>%
  relocate(Season, .after = Location)

MarchData <- totaldata %>%
  left_join(locations) %>%
  mutate(Season = "Wet",
         Date = mdy(Date)) %>%
  relocate(Season, .after = Location)

AllData<-bind_rows(AugData,MarchData)

AllData %>%
  ggplot(aes(x = TA*Salinity/36, y = DIC*Salinity/36, color = paste(Day_Night, Tide)))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(Location~Plate_Seep, scale = "free")

# Bring in the turb data from August
turbs<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Nutrients/Turb_NC.csv") %>%
  mutate(Season = "Dry") %>%
  select(CowTagID, Season, del15N, N_percent, C_N)


models2<-AllData %>%
  left_join(turbs)%>%
  filter(Plate_Seep == "Plate") %>%
  mutate(TA_sal = TA*Salinity/36, # salinity normalize
         DIC_sal = DIC*Salinity/36) %>%
  nest(data = -c(Location,CowTagID)) %>%
  mutate(fit = map(data, ~lm(TA_sal~DIC_sal, data = .)),
         coefs = map(fit, tidy)) %>%
  select(!fit)%>%
  unnest(cols = coefs) %>%
  filter(term == "DIC_sal") %>%
  unnest(cols = data) %>%
  group_by(Location,CowTagID)%>%
  summarise_at(vars(pH, Salinity, estimate, del15N, N_percent, C_N), .funs = ~min(.x,na.rm=TRUE))

models2 %>%
  ggplot(aes(x = del15N, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location, scales = "free")

# TA/DIC slope ~ del15N
mod15N<-lm(estimate~del15N, data = models2 %>%filter(Location =="Varari"))
anova(mod15N)

mod15N<-lm(estimate~del15N, data = models2 %>%filter(Location =="Cabral"))
anova(mod15N)

mod15N<-lm(estimate~del15N*Location, data = models2)
anova(mod15N)

AllData %>%
  mutate(TA_sal = TA*Salinity/36, # salinity normalize
         DIC_sal = DIC*Salinity/36) %>%
  ggplot(aes(x = TA_sal, y = DIC_sal))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*CowTagID, scale = "free")

