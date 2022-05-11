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
  ggplot(aes(x = TA, y = DIC, color = log(Salinity)))+
  geom_point()+
  facet_wrap(Location~Plate_Seep, scale = "free")

models<-totaldata %>%
  left_join(locations) %>%
  filter(Plate_Seep == "Plate") %>%
  nest(data = -c(Location,CowTagID)) %>%
  mutate(fit = map(data, ~lm(TA~DIC, data = .)),
         coefs = map(fit, tidy)) %>%
  select(!fit)%>%
  unnest(cols = coefs) %>%
  filter(term == "DIC") %>%
  unnest(cols = data) %>%
  group_by(Location,CowTagID)%>%
  summarise_at(vars(pH, Salinity, estimate), .funs = ~mean(.x,na.rm=TRUE))


models %>%
  ggplot(aes(x = pH, y = estimate))+
  geom_point()+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location, scales = "free")

  

ggplot(aes(x = TA, y = DIC))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*CowTagID, scale = "free")
