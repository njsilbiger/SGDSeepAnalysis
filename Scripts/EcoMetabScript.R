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
  ggplot(aes(x = log(Silicate_umolL), y = TA.diff, color = Tide_Time))+
  geom_smooth(method = "lm")+
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
  ggplot(aes(x = DIC.diff, y = TA.diff, color = Tide))+
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



# Silicate versus everything in the seep

Data_predictions %>%
  filter(Plate_Seep == "Seep") %>%
  select(Location, Day_Night, Tide, Tide_Time,Silicate_umolL, pH, Salinity, TA, DIC, NN_umolL, Ammonia_umolL, Phosphate_umolL, pCO2)%>%
  pivot_longer(cols = pH:pCO2, names_to = "Params", values_to = "Values") %>%
 # filter(Silicate_umolL<400)%>%
  drop_na()%>%
  ggplot(aes(x = Silicate_umolL, y = Values, color = Location))+
  geom_smooth(method = "lm")+
  geom_point()+
  labs(title = "Seep data")+
 # coord_trans(y = "log",x = "log")+
  facet_wrap(~Params, scales = "free")



#### relationships among parameters
# Si vs NN
Data_predictions %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = log(NN_umolL)))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

# Si vs Salinity
Data_predictions %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = log(Salinity)))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

# Si vs Temperature
Data_predictions %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = Temperature))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

# Si vs pH
Data_predictions %>% ## In Varari you see the effect in the correct direction at low night, everything else is flat
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = pH))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

# Si vs Ammonium
Data_predictions %>% ## In Varari you see the effect in the correct direction at low night, everything else is flat
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = log(Ammonia_umolL)))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

# Si vs DIC
Data_predictions %>% ## In Varari you see the effect in the correct direction at low night, everything else is flat
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = DIC, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

# Si vs TA
Data_predictions %>% ## In Varari you see the effect in the correct direction at low night, everything else is flat
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = TA, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

#  Temperature vs delta DIC
Data_predictions %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = Temperature, y = DIC.diff))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location, scales = "free")


# NN vs delta DIC
Data_predictions %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(NN_umolL), y = DIC.diff))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location, scales = "free")

# Nh4 vs delta DIC # Day night might have an interaction
Data_predictions %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Ammonia_umolL), y = DIC.diff))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location, scales = "free")


# delta DIC vs pH
Data_predictions %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = DIC.diff, y = pH))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location)

# pH vs delta TA
Data_predictions %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = pH, y = TA.diff))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location, scales = "free")

#  Temperature vs delta TA
Data_predictions %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = Temperature, y = TA.diff))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")





## What is the distribution of average silicate at the plates

SIMeans<-Data_predictions %>%
  filter(Plate_Seep=="Plate")%>%
  group_by(Location, CowTagID)%>%
  summarise(meanSI = mean(Silicate_umolL, na.rm=TRUE)) %>% # cakculate means
  mutate(QRank = factor(ntile(meanSI,2))) %>% # put the plates into 3 groups by quartiles
  right_join(Data_predictions) # join in back with the datapredictions

SIMeans%>%
  ggplot(aes(x = meanSI, y = Location))+
  ggridges::geom_density_ridges()

SIMeans %>%
  filter(Plate_Seep=="Plate")%>%
ggplot(aes(x = DIC, y = TA, color = QRank, group = QRank))+
  geom_point()+
  geom_smooth(method= "lm")+
  facet_wrap(~Location)

mod1<-lm(TA.pred~DIC.pred*QRank, SIMeans[SIMeans$Location=="Varari" & SIMeans$Plate_Seep=="Plate",])
anova(mod1)

mod2<-lm(TA.pred~DIC.pred*QRank, SIMeans[SIMeans$Location=="Cabral" & SIMeans$Plate_Seep=="Plate",])
anova(mod2)

SIMeans %>%
  filter(Plate_Seep=="Plate")%>%
  ggplot(aes(x = DIC.pred, y = TA.pred, color = Tide, group = Tide))+
  geom_point()+
  geom_smooth(method= "lm")+
  facet_wrap(~Location)

mod3<-lm(TA.pred~DIC.pred*Tide, SIMeans[SIMeans$Location=="Varari" & SIMeans$Plate_Seep=="Plate",])
anova(mod3)
summary(mod3)

mod4<-lm(TA.pred~DIC.pred*Tide, SIMeans[SIMeans$Location=="Cabral" & SIMeans$Plate_Seep=="Plate",])
anova(mod4)
summary(mod4)


SIMeans %>%
  filter(Plate_Seep=="Plate")%>%
  mutate(QRank_raw = factor(ntile(log(Silicate_umolL),3))) %>% # put the plates into 3 groups by quartiles
  ggplot(aes(x = DIC.pred, y = TA.pred, color = QRank_raw, group = QRank_raw))+
  geom_point()+
  geom_smooth(method= "lm")+
  facet_wrap(~Location)
