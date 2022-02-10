### TA script ####

####################################

# load libraries #################
library(tidyverse)
library(ggfortify)
library(lubridate)
library(ggforce)
library(viridis)
library(patchwork)
library(here)
library(ggridges)
library(seacarb)

# load the 24 hour chemistry data #####################
Data<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv")
turbdata<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Nutrients/Turb_NC.csv")

#### Visualize the TA data ####

## also filter out all the data from the first low tide where the water level was super high
remove2<-Data %>% filter(CowTagID=="V2", Tide =="Low", Day_Night=="Day", Date == ymd("2021-08-08"))

# make a set of density plots
Data %>%
  anti_join(remove2)%>%
  filter(Plate_Seep=="Plate") %>%
  mutate(Tide_Time = paste(Tide, Day_Night))%>%
  ggplot(aes(x = TA, y = Tide_Time)) +
    geom_density_ridges()+
    facet_wrap(~Location)

Data %>%
  anti_join(remove2)%>%
  filter(Plate_Seep=="Plate") %>%
  mutate(Tide_Time = paste(Tide, Day_Night))%>%
  ggplot(aes(x = Salinity, y = TA, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location)

### Calculate all carbonate parameters #############
#remove rows with NAs
Cdata<-Data %>%
  select(-c(Jamie_Plate_ID, Bottom_Plate_ID, Top_Plate_ID)) %>%
  mutate(Tide_Time = paste(Tide, Day_Night)) %>%
  drop_na(TA) # keep only complete cases
 

#calculate rest of carbonate params
CO2<-carb(flag=8, Cdata$pH, Cdata$TA/1000000, S=Cdata$Salinity, 
          T=Cdata$Temperature, Patm=1, P=0, Pt=Cdata$Phosphate_umolL/1000000, Sit=Cdata$Silicate_umolL/1000000, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")

#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propogation
er<-errors(flag=8, Cdata$pH, Cdata$TA/1000000, 
           S=Cdata$Salinity, T=Cdata$Temperature, 
           Patm=1, P=0,Pt=Cdata$Phosphate_umolL/1000000,
           Sit=Cdata$Silicate_umolL/1000000,evar1 = 0.01, evar2 = 5e-6) 

#average error for DIC based on pH and TA
mean(er$DIC*1000000)
sd(er$DIC*1000000)/sqrt(nrow(er))

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
CO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-CO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

Cdata[,c("CO2","HCO3","CO3","DIC","OmegaArag","OmegaCalcite","pCO2","fCO2")]<-
  CO2[,c("CO2","HCO3","CO3","DIC","OmegaAragonite","OmegaCalcite","pCO2","fCO2")]


### Pull out the endmembers ###
## FOr the seep, we will take the min salinity value as the end member, while we have well samples, I don't believe this is as representative because the water has not gone through the sediment yet.

## Use the well and spring as endmembers and try again 

Endmembers_HighSGD<-Cdata %>%
 # filter(Plate_Seep == "Seep") %>%
  filter(Plate_Seep == "Well"|Plate_Seep == "Spring") %>%
  group_by(Location) %>%
  filter(rank(Salinity, ties.method="first")==1) %>% # select the lowest salinity
  select(Location, Salinity_highend=Salinity, Temperature_highend = Temperature, TA_highend = TA, pHhighend = pH, Phosphate_umolL_highend = Phosphate_umolL, Silicate_umolL_highend = Silicate_umolL, NN_umolL_highend = NN_umolL, DIC_highend = DIC) 

# For now, take the MCR bottle samples from offshore. The last sample as from 8/12/2012... this is not good enough, but will hold until we get a good offshore sample.  

# calculate pH for endmember from TA and DIC
CO2_end<-carb(flag=15, 2385.3/1000000, 2029.6/1000000, S=36.224, 
          T=26.2962, Patm=1, P=0, Pt=0.13/1000000, Sit=0.38/1000000, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")


Endmembers_LowSGD<-tibble(Location = c("Varari", "Cabral"), Salinity_lowend = 36.224,Temperature_lowend=26.2962, TA_lowend=2385.3, pH_lowend=CO2_end$pH, Phosphate_umolL_lowend=0.13, Silicate_umolL_lowend=0.38, NN_umolL_lowend=0.13, DIC_lowend = 2029.6)

Endmembers<-left_join(Endmembers_HighSGD, Endmembers_LowSGD)

## join all the endmembers to calculate normalized values
#predicted TA based on mixing line
## use cristina's methods  C1 = Cmix + (Cmix – Csgd)(((Smix – 35.2)/(Ssgd – Smix))  

# Data_predictions<-Cdata %>%
#   left_join(Endmembers) %>%
#   mutate(TA.pred =  TA+(TA - TA_highend)*((Salinity - Salinity_lowend)/(Salinity_highend - Salinity)),
#          DIC.pred = DIC+(DIC - DIC_highend)*((Salinity - Salinity_lowend)/(Salinity_highend - Salinity)),
#          NN.pred = NN_umolL+(NN_umolL-NN_umolL_highend)*((Salinity - Salinity_lowend)/(Salinity_highend - Salinity)),
#          TA.diff = (TA_lowend - TA.pred)/2,
#          DIC.diff = DIC_lowend - DIC.pred) %>%
#   select(-ends_with("end")) # remove all the endmemmbers

Data_predictions<-Cdata %>%
  left_join(Endmembers) %>%
  mutate(TA.pred =  TA+(TA - TA_highend)*((Silicate_umolL - Silicate_umolL_lowend)/(Silicate_umolL_highend - Silicate_umolL)),
         DIC.pred = DIC+(DIC - DIC_highend)*((Silicate_umolL - Silicate_umolL_lowend)/(Silicate_umolL_highend - Silicate_umolL)),
         NN.pred = NN_umolL+(NN_umolL-NN_umolL_highend)*((Silicate_umolL - Silicate_umolL_lowend)/(Silicate_umolL_highend - Silicate_umolL)),
         TA.diff = (TA_lowend - TA.pred)/2,
         DIC.diff = DIC_lowend - DIC.pred) %>%
  select(-ends_with("end")) # remove all the endmemmbers

# Make some plots
Data_predictions %>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = Silicate_umolL, y = TA.diff, color = Tide_Time))+
  geom_point()+
  facet_wrap(~Location, scales = "free")

Data_predictions %>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = log(Silicate_umolL), y = DIC.diff, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

Data_predictions %>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = DIC.diff, y = TA.diff, color = Tide_Time, size = log(Silicate_umolL)))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")


Data_predictions %>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = DIC, y = TA, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")+
  theme_bw()


Data_predictions %>%
  filter(Plate_Seep == "Seep") %>%
  mutate(Tide = fct_relevel(Tide,c("Low","Mid","High")))%>%
  ggplot(aes(x = log(Silicate_umolL), y = pH, color = Day_Night))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(y = "pH")+
  facet_wrap(~Location, scales = "free")

Data_predictions %>%
  filter(Plate_Seep == "Plate") %>%
#  mutate(Tide = fct_relevel(Tide,c("Low","Mid","High")))%>%
  ggplot(aes(x = log(Silicate_umolL), y = DIC))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(y = "DIC")+
  facet_wrap(~Location*Day_Night, scales = "free")

Data_predictions %>%
  filter(Plate_Seep == "Seep") %>%
  #  mutate(Tide = fct_relevel(Tide,c("Low","Mid","High")))%>%
  ggplot(aes(x = log(Silicate_umolL), y = pH, color = Day_Night))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(y = "DIC")
 # facet_wrap(~Location, scales = "free") 

