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
library(ggtext)
library(lme4)
library(lmerTest)
library(broom)
library(oce)
library(margins)



# # load the 24 hour chemistry data #####################
# Data_Dry<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC2.csv") %>%
#   mutate(Season = "Dry") %>%
#   mutate(Date = mdy(Date))
# 
# # wet season
# Data_wet<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/CarbonateChemistry/pHProbe_Data_calculated_POcorrect.csv") %>%
#   mutate(Season = "Wet") %>%
#   rename(Temperature = TempInSitu,
#          Time = SamplingTime) %>%
#   mutate(DateTime = ymd_hms(paste(as.character(Date), as.character(Time)))) %>% # make a datetime column like the dry season
#   left_join(read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/Sandwich_Locations_Final.csv"))
# 
# ## Remove the one crazy Ammonium outlier
# remove2<-Data_Dry %>% filter(CowTagID=="V2", Tide =="Low", Day_Night=="Day", Date == ymd("2021-08-08"))
# removelow<- Data_Dry %>% # remove the not real low tide
#   filter(Date == ymd("2021-08-06") & Tide == "Low" & Plate_Seep == "Plate")
# turbdata<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Nutrients/Turb_NC.csv") %>%
#   mutate(Season = "Dry")
# 
# turb_wet<- read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/Nutrients/Turb_NC.csv") %>%
#   mutate(Season = "Wet")
# 
# turb_all<-bind_rows(turbdata, turb_wet)
# 
##---------------------------------------------

# read in the data from the pca script so that all data cleaning is the same
source(here("Scripts","PCA_Script.R"))


## load the distance to seep and depth data to calculate residence times
distance<-read_csv(here("Data","distance_depth.csv")) %>%
  mutate(Depth_m = adj_CT_depth_cm/100)%>%
  #select(CowTagID, dist_to_seep_m, Depth_m) 
  select(CowTagID, dist_to_crest,Depth_m)


# calculate average depth of the site
SiteDepth<-distance %>%
  left_join(Data)%>%
  group_by(Location) %>%
  summarise(meanDepth = mean(Depth_m, na.rm = TRUE)) %>%
  mutate(reeflength = ifelse(Location  == "Varari",189,100)) # or 272 for cabral if south to north


## load the current speed data
currentspeed<-read_csv(here("Data","currentspeed.csv")) %>%
  mutate(DateTime = ymd_hms(DateTime))


#### Visualize the TA data ####

# bring both seasons together and add in other important metrics for calculating metabolism
Data <- Data %>%
  #Data_Dry %>%
  #bind_rows(Data_wet) %>%
  left_join(distance) %>% # add in the distance and depth data
  left_join(currentspeed)%>%
  left_join(SiteDepth) %>%
 # select(!c(SeepCode:LoggerID,Top_Plate_ID:Jamie_Plate_ID)) %>% # remove unnecessary columns
  select(!c(Top_Plate_ID:Jamie_Plate_ID)) %>%
  mutate(residence_time_h = (1/CurrentSpeed_m_s_av)*reeflength*(1/60),
        # residence_time_h = (1/CurrentSpeed_m_s)*dist_to_seep_m*(1/60),  # calculate residence time of water moving from seep to cowtagID in hours
         rho = swRho(salinity = Salinity, temperature = Temperature, pressure = 0, eos="unesco")) # calculate seawater density in kg/m3

# make a set of density plots
Data %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep=="Plate") %>%
  mutate(Tide_Time = paste(Tide, Day_Night))%>%
  ggplot(aes(x = TA, y = Tide_Time, fill = Season)) +
    geom_density_ridges()+
    facet_wrap(~Location)

Data %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep=="Plate") %>%
  mutate(Tide_Time = paste(Tide, Day_Night))%>%
  ggplot(aes(x = Salinity, y = TA, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location* Season, scales = "free")

### Calculate all carbonate parameters #############
#remove rows with NAs
Cdata<-Data %>%
  #select(-c(Jamie_Plate_ID, Bottom_Plate_ID, Top_Plate_ID)) %>%
  mutate(Tide_Time = paste(Tide, Day_Night)) %>%
  drop_na(TA) # keep only complete cases
 

#calculate rest of carbonate params
CO2<-carb(flag=8, Cdata$pH, Cdata$TA/1000000, S=Cdata$Salinity, 
          T=Cdata$Temperature, Patm=1, P=0, 
       #   Pt=Cdata$Phosphate_umolL/1000000, Sit=Cdata$Silicate_umolL/1000000, 
          k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")

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

# From the Maunalua paper --- this only works for Varari... V13 is the endmember for upstream sample for delta TA... use average high tide salinity for "ocean sample"
# for each location and each season

# mutate(TA.pred = case_when(Location == 'Varari' & Season == 'Dry'~ TA+(TA - TA.seep.end)*((Salinity - High.Salinity)/(Seep.salinity - Salinity)),
                            #Location == 'Varari' & Season == 'Wet'~ TA+(TA - TA.seep.end)*((Salinity - High.Salinity)/(Seep.salinity - Salinity))
#                            )


### Calculate endmembers from the mixing line (remove the dug well and the street samples at Cabral)

Cdata %>% filter(Plate_Seep == "Seep"| Plate_Seep == "Spring") %>% filter( TA <5000 & TA >1500)

## calculate the endmembers from the springs
Vend_spring<-Cdata %>%
  filter(Plate_Seep == "Spring", CowTagID == "VSPRING") %>%
  select(Location,  TA_Spring = TA,DIC_Spring = DIC, Salinity_Spring = Salinity, Silicate_umolL_Spring = Silicate_umolL)%>%
  slice(1)

Vend_spring <- bind_rows(Vend_spring, Vend_spring) %>% # use the same end member for dry and wet
  mutate(Season = c("Dry","Wet"))
  

Cend_spring<-Cdata %>%
  filter(Plate_Seep == "Spring", Location == "Cabral", CowTagID %in% c("CSPRING", "CSPRING_ROAD") )%>% # was CSPRING_ROAD prior
  select(Location,Season, TA_Spring = TA, DIC_Spring = DIC, Salinity_Spring = Salinity, Silicate_umolL_Spring = Silicate_umolL)%>%
  slice(1:2)

Spring<-bind_rows(Vend_spring, Cend_spring)
## Calculate endmembers for the "ocean" which will be the average salinity or silicate during high tide by season

endmembers<- Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  drop_na(Silicate_umolL)%>%
   group_by(Location, Season, Tide) %>%
  # summarise_at(vars(Salinity, Silicate_umolL), .funs = list(median)) %>%
   summarise_at(vars(Salinity, Silicate_umolL), .funs = list(min, max)) %>%
   filter(Tide == "High") %>%
   rename(Salinity_High = Salinity_fn2, Silicate_umolL_High = Silicate_umolL_fn1) %>%
   select(Location, Season ,Salinity_High, Silicate_umolL_High) %>%
   left_join(Spring) # bring in the spring data


Cdata <- left_join(Cdata, endmembers)  # add the endmemebers to the dataframe to make it easier to calculate
  
## Select the VSpring  which is after it had gone through the sand. The two values are almost identica

## instead of using end members, use the regression between TA and salinity at the seep to make a mixing line and predict what
## TA should be at the plate given the salinity from the line. Then calculate deviations from the line to account for biology
## Varari is really clean and the same across seasons and Cabral, of course, is a mess.

# some plots of the relationship betwen TA, DIC, NN, and PO versus salinity at the seep
Cdata %>%
  filter(Plate_Seep == "Seep"| Plate_Seep == "Spring") %>% 
  filter( TA <5000 & TA >1500) %>% # remove the road samples
  ggplot(aes(x = Salinity, y = TA, color = Season))+
  geom_point(aes(color = Season))+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scale = "free")

Cdata %>%
  filter(Plate_Seep == "Seep"| Plate_Seep == "Spring") %>% 
  filter( TA <5000 & TA >1500) %>% # remove the road samples
  ggplot(aes(x = Salinity, y = pH, color = Season))+
  geom_smooth(method = "lm")+
  geom_point()+facet_wrap(~Location, scale = "free")

Cdata %>%
  filter(Plate_Seep == "Seep") %>% 
  filter( TA <5000 & TA >1500) %>% # remove the road samples
  ggplot(aes(x = Salinity, y = NN_umolL, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scale = "free")

Cdata %>%
  filter(Plate_Seep == "Seep") %>% 
  filter( TA <5000 & TA >1500) %>% # remove the road samples
  ggplot(aes(x = Salinity, y = Phosphate_umolL, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scale = "free")

# Create a TA mixing line--- random slope and intercept for each season
VarariMixModel<-lmer(TA~Salinity+(Salinity|Season), data = Cdata %>% filter(Location == "Varari", Plate_Seep == "Seep"))
Vco<-coef(VarariMixModel)

CabralMixModel<-lmer(TA~Salinity+(Salinity|Season), data = Cdata %>% filter(Location == "Cabral", Plate_Seep == "Seep"))
Cco<-coef(CabralMixModel)

# DIC mixing line
VarariMixModelDIC<-lmer(DIC~Salinity+(Salinity|Season), data = Cdata %>% filter(Location == "Varari", Plate_Seep == "Seep"))
VcoDIC<-coef(VarariMixModelDIC)

CabralMixModelDIC<-lmer(DIC~Salinity+(Salinity|Season), data = Cdata %>% filter(Location == "Cabral", Plate_Seep == "Seep"))
CcoDIC<-coef(CabralMixModelDIC)

# NN mixing line
VarariMixModelNN<-lm(NN_umolL~Salinity, data = Cdata %>% filter(Location == "Varari", Plate_Seep == "Seep"))
VcoNN<-coef(VarariMixModelNN)

CabralMixModelNN<-lm(NN_umolL~Salinity, data = Cdata %>% filter(Location == "Cabral", Plate_Seep == "Seep"))
CcoNN<-coef(CabralMixModelNN)

# PO mixing line
VarariMixModelPO<-lm (Phosphate_umolL~Salinity, data = Cdata %>% filter(Location == "Varari", Plate_Seep == "Seep"))
VcoPO<-coef(VarariMixModelPO)

CabralMixModelPO<-lm (Phosphate_umolL~Salinity, data = Cdata %>% filter(Location == "Cabral", Plate_Seep == "Seep"))
CcoPO<-coef(CabralMixModelPO)

## put them in the dataframe
Cdata <- Cdata %>% # add the predicted mixing line
  # mutate(TA.mix = case_when(Location  == "Varari" & Season == "Dry" ~Salinity*Vco$Season[1,2]+Vco$Season[1,1],
  #                           Location  == "Varari" & Season == "Wet" ~Salinity*Vco$Season[2,2]+Vco$Season[2,1],
  #                           Location  == "Cabral" & Season == "Dry" ~Salinity*Cco$Season[1,2]+Cco$Season[1,1],
  #                           Location  == "Cabral" & Season == "Wet" ~Salinity*Cco$Season[2,2]+Cco$Season[2,1]
  #                           )) %>%
  # mutate(DIC.mix = case_when(Location  == "Varari" & Season == "Dry" ~Salinity*VcoDIC$Season[1,2]+VcoDIC$Season[1,1],
  #                                    Location  == "Varari" & Season == "Wet" ~Salinity*VcoDIC$Season[2,2]+VcoDIC$Season[2,1],
  #                                    Location  == "Cabral" & Season == "Dry" ~Salinity*CcoDIC$Season[1,2]+CcoDIC$Season[1,1],
  #                                    Location  == "Cabral" & Season == "Wet" ~Salinity*CcoDIC$Season[2,2]+CcoDIC$Season[2,1]
  #         )) %>%
   mutate(TA.mix = TA+(TA - TA_Spring)*((Silicate_umolL - Silicate_umolL_High)/(Silicate_umolL_Spring - Silicate_umolL)), # account for SGD mixing
          DIC.mix = DIC+(DIC - DIC_Spring)*((Silicate_umolL - Silicate_umolL_High)/(Silicate_umolL_Spring - Silicate_umolL))                                  
           ) 


# For Varari, Use Ta from V13 as the "incoming water"... ignore Cabral for now

incoming<-Cdata %>%
  filter(CowTagID == "V13") %>%
  select(Location,Season, Date, Tide, Day_Night, TA_incom = TA.mix,DIC_incom = DIC.mix)


## calculate average TA and DIC from offshore water/ VRC and CRS
crestwater<-Cdata %>%
  filter(Plate_Seep == "Offshore", DIC<2060) %>%
  group_by(Location, TimeBlock) %>%
  #group_by(Location, Day_Night) %>%
           #, Day_Night, Tide) %>%
  summarise(TA.offshore = mean(TA, na.rm = TRUE),
            DIC.offshore = mean(DIC, na.rm = TRUE)) %>%
  drop_na()

Cdata<-Cdata %>%
  left_join(incoming) %>%
  left_join(crestwater)%>%
         #   ifelse(Location  == "Varari",Salinity*Vco[2]+Vco[1],Salinity*Cco[2]+Cco[1]), # mixing line TA
         # TA.diff = TA.mix-TA, # observed - predicted to see how biology changes TA above what is expected by mixing
         # DIC.mix = ifelse(Location  == "Varari",Salinity*VcoDIC[2]+VcoDIC[1],Salinity*CcoDIC[2]+CcoDIC[1]),
         # DIC.diff = DIC.mix-DIC,
        mutate(
          # TA.diff = TA.mix- TA,
          # DIC.diff = DIC.mix - DIC,
          TA.diff = TA.offshore - TA.mix, #TA_incom --- take average from offshore water 2360
          DIC.diff = DIC.offshore - DIC.mix, # DIC_income 2000
          #  NN.mix = ifelse(Location  == "Varari",Salinity*VcoNN[2]+VcoNN[1],Salinity*CcoNN[2]+CcoNN[1]),
         # NN.diff =  NN.mix - NN_umolL,
         # PO.mix = ifelse(Location  == "Varari",Salinity*VcoPO[2]+VcoPO[1],Salinity*CcoPO[2]+CcoPO[1]),
         # PO.diff =  PO.mix -Phosphate_umolL,
         NEC = ((rho*meanDepth*TA.diff)/(2*residence_time_h))/1000, #mmol m-2 hr-1 0.575 is the average depth of the site
         NEP = ((rho*meanDepth*(DIC.diff - (TA.diff/2)))/residence_time_h)/1000, #mmol m-2 hr-1, rho = 1025
         NEC.proxy = TA.diff/2,
         NEP.proxy = DIC.diff - (TA.diff/2)
        # NEC = ((rho*0.575*TA.diff)/(2*residence_time_h))/1000, #mmol m-2 hr-1 0.575 is the average depth of the site
        # NEP = ((rho*0.575*(DIC.diff - (TA.diff/2)))/residence_time_h)/1000 #mmol m-2 hr-1
         )

Cdata %>%
  filter(CowTagID %in% c("V12", "V15")) %>%
  select(CowTagID, Tide, Time, Season, TA.mix, TA_incom, DIC_incom, DIC.mix, Salinity, NEC, NEP, residence_time_h)

# Endmembers_HighSGD<-Cdata %>%
#  # filter(Plate_Seep == "Seep") %>%
#   filter(Plate_Seep == "Well"|Plate_Seep == "Spring") %>%
#   group_by(Location) %>%
#   filter(rank(Salinity, ties.method="first")==1) %>% # select the lowest salinity
#   select(Location, Salinity_highend=Salinity, Temperature_highend = Temperature, TA_highend = TA, pHhighend = pH, Phosphate_umolL_highend = Phosphate_umolL, Silicate_umolL_highend = Silicate_umolL, NN_umolL_highend = NN_umolL, DIC_highend = DIC) 
# 
# # For now, take the MCR bottle samples from offshore. The last sample as from 8/12/2012... this is not good enough, but will hold until we get a good offshore sample.  
# 
# # calculate pH for endmember from TA and DIC
# CO2_end<-carb(flag=15, 2385.3/1000000, 2029.6/1000000, S=36.224, 
#           T=26.2962, Patm=1, P=0, Pt=0.13/1000000, Sit=0.38/1000000, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
# 
# 
# Endmembers_LowSGD<-tibble(Location = c("Varari", "Cabral"), Salinity_lowend = 36.224,Temperature_lowend=26.2962, TA_lowend=2385.3, pH_lowend=CO2_end$pH, Phosphate_umolL_lowend=0.13, Silicate_umolL_lowend=0.38, NN_umolL_lowend=0.13, DIC_lowend = 2029.6)
# 
# Endmembers<-left_join(Endmembers_HighSGD, Endmembers_LowSGD)

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

# Data_predictions<-Cdata %>%
#   left_join(Endmembers) %>%
#   mutate(TA.pred =  TA+(TA - TA_highend)*((Silicate_umolL - Silicate_umolL_lowend)/(Silicate_umolL_highend - Silicate_umolL)),
#          DIC.pred = DIC+(DIC - DIC_highend)*((Silicate_umolL - Silicate_umolL_lowend)/(Silicate_umolL_highend - Silicate_umolL)),
#          NN.pred = NN_umolL+(NN_umolL-NN_umolL_highend)*((Silicate_umolL - Silicate_umolL_lowend)/(Silicate_umolL_highend - Silicate_umolL)),
#          TA.diff = (TA_lowend - TA.pred)/2,
#          DIC.diff = DIC_lowend - DIC.pred) %>%
#   select(-ends_with("end")) # remove all the endmemmbers

# Make some plots
Cdata %>%
  #anti_join(removelow)%>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = log(Silicate_umolL), y = NEC.proxy, color = Tide_Time))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>%
  filter(Plate_Seep == "Plate") %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  ggplot(aes(x = log(Silicate_umolL), y = NEP.proxy, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = NEP.proxy, y = NEC.proxy, color = Tide, shape = Day_Night))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = NEP, y = NEC, color = Tide, shape = Day_Night))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = NEP.proxy, y = NEC.proxy, color = Tide,
             shape = Day_Night))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

# make a model to see if the slopes are different
TADICdiffmod<-lm(TA.diff/2~DIC.diff*Tide , data  = Cdata %>% anti_join(remove_varari)%>%
                   anti_join(remove_cabral)%>% filter(Location  == "Varari", Plate_Seep == "Plate", Season == "Dry"))
anova(TADICdiffmod)
summary(TADICdiffmod)

#Easy coefficients
TADICdiffcoef<-lm(TA.diff/2~Tide/DIC.diff -1, data  = Cdata %>% anti_join(remove_varari)%>%
                    anti_join(remove_cabral)%>% filter(Location  == "Varari", Plate_Seep == "Plate", Season == "Dry"))
coef(TADICdiffcoef)

# Have a random effect for day and night to allow the intercept to change 
TADICdiffmod_rand<-lmer(TA.diff/2~DIC.diff*Tide*Season + (1|Day_Night) , data  = Cdata %>% anti_join(remove_varari)%>%
                          anti_join(remove_cabral)%>% filter(Location  == "Varari", Plate_Seep == "Plate"))

anova(TADICdiffmod_rand)
summary(TADICdiffmod_rand)

TADICdiffcoef_rand<-lmer(TA.diff/2~Season/Tide/DIC.diff + (1|Day_Night) , data  = Cdata %>% anti_join(remove_varari)%>%
                           anti_join(remove_cabral)%>% filter(Location  == "Varari", Plate_Seep == "Plate"))
coef(TADICdiffcoef_rand)


Cdata %>%
  filter(Plate_Seep == "Plate") %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  ggplot(aes(x = DIC, y = TA, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")+
  theme_bw()

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  drop_na(NEC)%>%
  ggplot(aes(x = pH, y = NEC, color = Tide))+
  geom_point(aes(shape = Day_Night))+
#  geom_hline(yintercept = 0)+
#  geom_vline(xintercept = 0)+
  geom_smooth(method = "lm")+
 facet_wrap(~Location*Season*Tide, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  drop_na(NEC)%>%
  ggplot(aes(x = pH, y = NEC.proxy, color = Tide))+
  geom_point(aes(shape = Day_Night))+
  #  geom_hline(yintercept = 0)+
  #  geom_vline(xintercept = 0)+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season*Tide, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  drop_na(NEC)%>%
  ggplot(aes(x = NEC.proxy, y = NEC, color = Tide))+
  geom_point()+
  #  geom_hline(yintercept = 0)+
  #  geom_vline(xintercept = 0)+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season*Tide, scales = "free")

# Silicate versus everything in the seep

Cdata %>%
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
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = log(NN_umolL)))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

# Si vs Salinity
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(y = log(Silicate_umolL), x = Salinity))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

# Si vs Temperature
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = Temperature))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

# Si vs pH
Cdata %>% ## In Varari you see the effect in the correct direction at low night, everything else is flat
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = pH, color = Tide, shapre = Day_Night))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

# Si vs Ammonium
Cdata %>% ## In Varari you see the effect in the correct direction at low night, everything else is flat
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = log(Ammonia_umolL)))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

# Si vs DIC
Cdata %>% ## In Varari you see the effect in the correct direction at low night, everything else is flat
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = DIC, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

# Si vs TA
Cdata %>% ## In Varari you see the effect in the correct direction at low night, everything else is flat
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Silicate_umolL), y = TA, color = Tide_Time))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

#  Temperature vs delta DIC
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = Temperature, y = DIC.diff, color = Day_Night))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = Temperature, y = NEP.proxy))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location*Season, scales = "free")


# NN vs delta DIC
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(NN_umolL), y = NEP.proxy, shape = Tide_Time))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location*Season, scales = "free")

# Nh4 vs delta DIC # Day night might have an interaction
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = log(Ammonia_umolL), y = NEP.proxy))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location*Season, scales = "free")


# delta DIC vs pH
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = NEP.proxy, y = pH, color = Season))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location, scales = "free")

Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = NEP, y = pH, color = Season))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location*Season, scales = "free")

# pH vs delta TA
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = pH, y = NEC.proxy, color = Season))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location, scales = "free")

Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = pH, y = NEC, color = Season))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location, scales = "free")

#  Temperature vs delta TA
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  ggplot(aes(x = Temperature, y = NEC.proxy, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

## What is the distribution of average silicate at the plates

SIMeans<-Cdata %>%
  filter(Plate_Seep=="Plate")%>%
  group_by(Location, CowTagID)%>%
  summarise(meanSI = mean(Silicate_umolL, na.rm=TRUE)) %>% # cakculate means
  mutate(QRank = factor(ntile(meanSI,2))) %>% # put the plates into 3 groups by quartiles
  right_join(Cdata) # join in back with the datapredictions

SIMeans%>%
  ggplot(aes(x = meanSI, y = Location, fill = Season))+
  ggridges::geom_density_ridges()

SIMeans %>%
  filter(Plate_Seep=="Plate")%>%
ggplot(aes(x = DIC.diff, y = TA.diff, color = QRank, group = QRank))+
  geom_point()+
  geom_smooth(method= "lm")+
  facet_wrap(~Location*Season, scales = "free")

mod1<-lm(TA.diff~DIC.diff*QRank*Season, SIMeans[SIMeans$Location=="Varari" & SIMeans$Plate_Seep=="Plate",])
anova(mod1)

mod2<-lm(TA.diff~DIC.diff*QRank*Season, SIMeans[SIMeans$Location=="Cabral" & SIMeans$Plate_Seep=="Plate",])
anova(mod2)

SIMeans %>%
  filter(Plate_Seep=="Plate")%>%
  ggplot(aes(x = DIC.diff, y = TA.diff, color = Tide, group = Tide))+
  geom_point()+
  geom_smooth(method= "lm")+
  facet_wrap(~Location*Season, scales = "free")

SIMeans %>%
  filter(Plate_Seep=="Plate")%>%
  ggplot(aes(x = NEP, y = NEC, color = Tide, group = Tide))+
  geom_point()+
  geom_smooth(method= "lm")+
  facet_wrap(~Location*Season, scales = "free")

mod3<-lm(TA.diff~DIC.diff*Tide*Season, SIMeans[SIMeans$Location=="Varari" & SIMeans$Plate_Seep=="Plate",])
anova(mod3)
summary(mod3)

mod4<-lm(TA.diff~DIC.diff*Tide*Season, SIMeans[SIMeans$Location=="Cabral" & SIMeans$Plate_Seep=="Plate",])
anova(mod4)
summary(mod4)


SIMeans %>%
  filter(Plate_Seep=="Plate")%>%
  mutate(QRank_raw = factor(ntile(log(Silicate_umolL),3))) %>% # put the plates into 3 groups by quartiles
  ggplot(aes(x = DIC.diff, y = TA.diff, color = QRank_raw, group = QRank_raw))+
  geom_point()+
  geom_smooth(method= "lm")+
  facet_wrap(~Location*Season, scales = "free")


######## plot relationship between the first and second low day ###

lowVarari<-Cdata %>%
  filter(Tide == "Low", 
         TimeBlock =="Morning", 
         Plate_Seep =="Plate",
         Location == "Varari",
         Season == "Dry") %>% 
  mutate( SGDpres = ifelse(Date == "8/6/2021", "SGD supressed", "SGD present")) 


PTADIC<-lowVarari %>%
  ggplot(aes(x = DIC.diff, y = TA.diff/2, color = SGDpres, fill = SGDpres))+
  geom_hline(aes(yintercept = 0), lty = 2)+
  geom_vline(aes(xintercept = 0), lty = 2)+
  geom_point()+
  geom_smooth(method = "lm") +
#  xlim(2000, 2105)+
 # scale_size_continuous(trans = "log10")+
  labs(#title = "Daytime low tides at Varari",
        #title = "Data collected between 6:40 - 7:40am",
       color = " ",
       x = "&Delta; DIC",
       y = "&Delta; TA/2") + 
  geom_curve(
 #   aes(x = 2020, y = 2380, xend = 2035, yend = 2360),
    aes(x = -20, y = 10, xend = -10, yend = 5),
    
    curvature = -0.5,
    arrow = arrow(
        length = unit(0.03, "npc"), 
      type="closed" # Describes arrow head (open or closed)
    ),
    colour = "#9CC3D5FF",
    size = 1.2,
    angle = 90 # Anything other than 90 or 0 can look unusual
  )+
  geom_curve(
   # aes(x = 2060, y = 2405, xend = 2075, yend = 2385),
    aes(x = 13, y = 18, xend = 23, yend = 13),
    
    curvature = 0.5,
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed" # Describes arrow head (open or closed)
    ),
    colour = "#0063B2FF",
    size = 1.2,
    angle = 90 # Anything other than 90 or 0 can look unusual
  )+
  annotate("text",x = -20, y = 12, label = "SGD present", size = 8)+
  #annotate("text",x = 2020, y = 2385, label = "SGD present", size = 8)+
  annotate("text",x = 13, y = 20, label = "SGD suppressed", size = 8)+
  #annotate("text",x = 2060, y = 2410, label = "SGD suppressed", size = 8)+
  scale_color_manual(values = c("#9CC3D5FF","#0063B2FF"))+
  scale_fill_manual(values = c("#9CC3D5FF","#0063B2FF"))+
  theme_bw()+
  theme(legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_markdown(size = 18),
        axis.title.y = element_markdown(size = 18),
        axis.text = element_markdown(size = 16))

# box plot or density plot of the distribution of silicate 
lowVarari %>%
  ggplot(aes(x = Silicate_umolL, fill = SGDpres))+
  geom_density(alpha = 0.5)+
  scale_color_manual(values = c("#9CC3D5FF","#0063B2FF"))+
  scale_fill_manual(values = c("#9CC3D5FF","#0063B2FF"))+
  theme_bw()

pbox<-lowVarari %>%
  ggplot(aes(x = fct_reorder(SGDpres, NN_umolL, .desc = FALSE), y = NN_umolL, fill = SGDpres))+
  geom_boxplot(alpha = 0.5)+
  geom_jitter(width = 0.1, shape = 21)+
  theme_bw()+
  scale_color_manual(values = c("#9CC3D5FF","#0063B2FF"))+
  scale_fill_manual(values = c("#9CC3D5FF","#0063B2FF"))+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.y = element_markdown(size = 16),
        axis.text.y = element_text(size = 14))+
  labs(x = "",
       y = "Nitrate (&mu;mol L<sup>-1</sup>)")


## Make a plot with an inset
PTADIC + annotation_custom(ggplotGrob(pbox), xmin = 15,
                           xmax = 45, ymin = -20, ymax = 5)
ggsave(here("Output","LowTideTADIC.png"), width = 10, height = 8)

# run an ANCOVA to see of the slopes are different
mod.low<-  lm(TA.diff~DIC.diff*SGDpres, data = lowVarari)
anova(mod.low)

# t-test for difference in NN concentration between the two lows
mod.lowNN<-  lm(NN_umolL~SGDpres, data = lowVarari)
anova(mod.lowNN)


### All TA vs DIC for Varari
# Cdata %>%
#   filter(Plate_Seep =="Plate",
#          Location == "Varari",
#          Date != ymd("2021-08-06")) %>%
#   ggplot(aes(x = DIC.diff, y = TA.diff, color = Tide))+
#   geom_point(aes(shape = Day_Night))+
#   geom_smooth(method = "lm", aes(fill = Tide)) +
#   #facet_wrap(~Season)
#   # scale_size_continuous(trans = "log10")+
#   labs(
#        #title = "Data collected between 6:40 - 7:40am",
#        color = "Tide",
#        x = "DIC normalized to silicate",
#        y = "TA normalized to silicate") + 
#   scale_shape_manual(values = c(22,16))+
#   scale_colour_hue(l = 45)+
#   scale_fill_hue(l = 45)+
#   geom_curve(
#     aes(x = 1990, y = 2380, xend = 2000, yend = 2360),
#     curvature = -0.5,
#     arrow = arrow(
#       length = unit(0.03, "npc"), 
#       type="closed" # Describes arrow head (open or closed)
#     ),
#     colour = "grey",
#     size = 1.2,
#     angle = 90 # Anything other than 90 or 0 can look unusual
#   )+
#   geom_curve(
#     aes(x = 2040, y = 2405, xend = 2055, yend = 2385),
#     curvature = 0.5,
#     arrow = arrow(
#       length = unit(0.03, "npc"), 
#       type="closed" # Describes arrow head (open or closed)
#     ),
#     colour = "grey",
#     size = 1.2,
#     angle = 90 # Anything other than 90 or 0 can look unusual
#   )+
#   geom_curve(
#     aes(x = 1940, y = 2270, xend = 1925, yend = 2285),
#     curvature = 0.5,
#     arrow = arrow(
#       length = unit(0.03, "npc"), 
#       type="closed" # Describes arrow head (open or closed)
#     ),
#     colour = "grey",
#     size = 1.2,
#     angle = 90 # Anything other than 90 or 0 can look unusual
#   )+geom_curve(
#     aes(x = 2000, y = 2285, xend = 1985, yend = 2300),
#     curvature = 0.5,
#     arrow = arrow(
#       length = unit(0.03, "npc"), 
#       type="closed" # Describes arrow head (open or closed)
#     ),
#     colour = "grey",
#     size = 1.2,
#     angle = 90 # Anything other than 90 or 0 can look unusual
#   )+
#   annotate("text",x = 1990, y = 2385, label = "High Tide (Noon)", size = 8)+
#   annotate("text",x = 2040, y = 2410, label = "Low Tide (Dawn)", size = 8)+
#   annotate("text",x = 1940, y = 2265, label = "Low Tide (Dusk)", size = 8)+
#   annotate("text",x = 2000, y = 2280, label = "High Tide (Midnight)", size = 8)+
#   theme_bw()+
#   
#   theme(legend.position="none",
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 18),
#         axis.text = element_text(size = 16))
# ggsave(here("Output","TADICallV.png"), width = 10, height = 10)

plot_all<-Cdata %>%
  filter(Plate_Seep =="Plate",
         Location == "Varari",
         Date != ymd("2021-08-06"), 
         Season == "Dry") %>%
  ggplot(aes(x = DIC.diff, y = TA.diff/2, color = Tide))+
  geom_hline(aes(yintercept = 0), lty = 2)+
  geom_vline(aes(xintercept = 0), lty = 2)+
  geom_point(aes(shape = Day_Night))+
  geom_smooth(method = "lm", aes(fill = Tide)) +
  #facet_wrap(~Season)
  # scale_size_continuous(trans = "log10")+
  labs(
    #title = "Data collected between 6:40 - 7:40am",
    color = "Tide",
    x = "&Delta; DIC",
    y = "&Delta; TA/2") + 
  scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  geom_curve(
    aes(x = -100, y = -6, xend = -90, yend = -14),
    curvature = -0.5,
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed" # Describes arrow head (open or closed)
    ),
    colour = "grey",
    size = 1.2,
    angle = 90 # Anything other than 90 or 0 can look unusual
  )+
  geom_curve(
    aes(x = -25, y = -22, xend = -40, yend = -12),
    curvature = 0.5,
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed" # Describes arrow head (open or closed)
    ),
    colour = "grey",
    size = 1.2,
    angle = 90 # Anything other than 90 or 0 can look unusual
  )+
  annotate("text",x = -100, y = -4, label = "Low Tide", size = 8)+
  annotate("text",x = -25, y = -28, label = "High Tide", size = 8)+
    theme_bw()+
  
  theme(legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_markdown(size = 18),
        axis.title.y = element_markdown(size = 18),
        axis.text = element_text(size = 16))
#ggsave(here("Output","TADICallV.png"), width = 10, height = 10)


HLmod<- lm(TA.diff/2~DIC.diff*Tide,data = Cdata %>%
               filter(Plate_Seep =="Plate",
                      Location == "Varari",
                      Date != ymd("2021-08-06"),
                      Season == "Dry") )


anova(HLmod)
summary(HLmod)



# to make a coefficients plots
HLmodnested<- lm(TA.diff/2~Tide/DIC.diff-1,data = Cdata %>%
             filter(Plate_Seep =="Plate",
                    Location == "Varari",
                    Date != ymd("2021-08-06")) )
coeffs<-tidy(HLmodnested)[3:4,]  %>% # just take the slopes
  mutate(#Day_Night = c("Day","Night","Dawn","Dusk"),
         #Day_Night2 = c("Day","Night","Day","Night"),
         Tide = c("High Tide","Low Tide"))
        # NiceLabels  = paste(Tide,":",Day_Night))

ggplot(coeffs, aes(x=estimate, y = Tide, color = Tide))+
  geom_point(size = 5)+
  geom_errorbarh(aes(xmin =estimate-std.error, xmax = estimate+std.error ),height = 0.1 )+
#  geom_vline(aes(xintercept = 0), lty = 2)+
#  scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  labs(x = "Coefficient",
       y = "")+
  theme_bw()+
  theme(legend.position="none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))
  


## boxplot with all of them
plotN_all<-Cdata %>%
  filter(Plate_Seep =="Plate",
         Location == "Varari",
         Date != ymd("2021-08-06"),
         Season == "Dry") %>%
  ggplot(aes(x = Tide, y = NN_umolL, color = Tide))+
  geom_boxplot(alpha = 0.5, aes( fill = Tide))+
  geom_jitter(width = 0.1)+
 # scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  theme_bw()+
  theme(legend.position = "none",
   axis.text.x = element_blank(),
   axis.text.y = element_text(size = 16),
   axis.title.y = element_markdown(size = 18))+
  labs(x = "",
       y = "Nitrate (&mu;mol L<sup>-1</sup>)")

## Make a plot with an inset
plot_all + annotation_custom(ggplotGrob(plotN_all), xmin = -50,
                           xmax = 10, ymin = -15, ymax = -40)
ggsave(here("Output","HighLowTADIC.png"), width = 10, height = 8)



### same plots with high tide included
Cdata %>%
  filter(#Tide == "Low", 
         TimeBlock %in% c("Morning", "Afternoon"), 
         Plate_Seep =="Plate",
         Location == "Varari",
         Season == "Dry") %>%
  mutate(SGDpres = case_when(Date == "8/6/2021" & Tide == "Low" ~"SGD supressed",
                             Date == "8/8/2021" & Tide == "Low" ~"SGD present",
                             Tide == "High" ~"High Tide",
                             
    
  )) %>%
  ggplot(aes(x = NEP, y = NEC, color = SGDpres, fill = SGDpres))+
  geom_point()+
  geom_smooth(method = "lm")


## boxplot
Cdata %>%
  filter(#Tide == "Low", 
    TimeBlock %in% c("Morning", "Afternoon"), 
    Plate_Seep =="Plate",
    Location == "Varari",
    Season == "Dry") %>%
  mutate(SGDpres = case_when(Date == "8/6/2021" & Tide == "Low" ~"SGD supressed",
                             Date == "8/8/2021" & Tide == "Low" ~"SGD present",
                             Tide == "High" ~"High Tide",
                             
                             
  )) %>%
ggplot(aes(x = SGDpres, y = NN_umolL, fill = SGDpres))+
  geom_boxplot(alpha = 0.5)+
  geom_jitter(width = 0.1, shape = 21)+
  theme_bw()+
#  scale_color_manual(values = c("#9CC3D5FF","#0063B2FF"))+
#  scale_fill_manual(values = c("#9CC3D5FF","#0063B2FF"))+
  theme(legend.position = "none")+
       # axis.text.x = element_blank())+
  labs(x = "",
       y = "Nitrate (umol L-1)")


###### Look at TA/DIC relationships at each cowtag

models<-Cdata %>%
  filter(Plate_Seep == "Plate") %>%
        # Location == "Varari") %>%
  
  nest(data = -c(Location,CowTagID, Season)) %>%
  mutate(fit = map(data, ~lm(NEC~NEP, data = .)),
        # fit = map(data, ~lm(TA.diff/2~DIC.diff, data = .)),
         coefs = map(fit, tidy)) %>%
  select(!fit)%>%
  unnest(cols = coefs) %>%
 #filter(term == "DIC.diff") %>%
  filter(term == "NEP") %>%
  unnest(cols = data) %>%
  group_by(Location,CowTagID, Season, estimate)%>%
  #summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL), .funs = ~(max(.x,na.rm=TRUE)))%>%
#  summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL, del15N, N_percent, Phosphate_umolL, Temperature, Ammonia_umolL), .funs = ~(mean(.x,na.rm=TRUE)))%>%
#  summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL, del15N, N_percent, Phosphate_umolL, Temperature, Ammonia_umolL), .funs = ~(sd(.x,na.rm=TRUE)/mean(.x,na.rm=TRUE)))%>%
  summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL, Phosphate_umolL, Temperature, Ammonia_umolL), .funs = ~(max(.x,na.rm=TRUE)-min(.x,na.rm=TRUE)))%>%
#  summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL, del15N, N_percent, Phosphate_umolL, Temperature), .funs = ~(mean(.x,na.rm=TRUE)))%>%
  
#  summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL, del15N, N_percent), .funs = ~(max(.x,na.rm=TRUE)-min(.x,na.rm=TRUE)))%>%
  #summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL, del15N, N_percent), .funs = ~(mean(.x,na.rm=TRUE)))%>%
  ungroup() %>%
  left_join(turb_all) # join in the turb data


# models across both seasons

models_both<-Cdata %>%
  filter(Plate_Seep == "Plate") %>%
   # join in the turb data
  nest(data = -c(Location,CowTagID)) %>%
  mutate(fit = map(data, ~lm(NEC.proxy~NEP.proxy, data = .)),
         # fit = map(data, ~lm(TA.diff/2~DIC.diff, data = .)),
         coefs = map(fit, tidy)) %>%
  select(!fit)%>%
  unnest(cols = coefs) %>%
  #filter(term == "DIC.diff") %>%
  filter(term == "NEP.proxy") %>%
  unnest(cols = data) %>%
  group_by(Location,CowTagID, estimate)%>%
  #summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL), .funs = ~(max(.x,na.rm=TRUE)))%>%
  summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL,  Phosphate_umolL, Temperature, Ammonia_umolL), .funs = ~(sd(.x,na.rm=TRUE)/mean(.x,na.rm=TRUE)))%>%
  #summarise_at(vars(pH, Salinity, NN_umolL, Silicate_umolL, del15N, N_percent, Phosphate_umolL, Temperature, Ammonia_umolL), .funs = ~(var(.x,na.rm=TRUE)))%>%
  ungroup() %>%
  left_join(turb_all) 
  

models %>%
  ggplot(aes(x = NN_umolL, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location*Season, scales = "free")

mod_NN<-lm(estimate ~ NN_umolL*Season*Location, data = models)
anova(mod_NN)

models_both %>%
  ggplot(aes(x = Silicate_umolL, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location, scales = "free")

mod_silicate<-lm(estimate ~ Silicate_umolL*Location, data = models_both)
anova(mod_silicate)


models %>%
  ggplot(aes(x = Temperature, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location*Season, scales = "free")

models_both %>%
  ggplot(aes(x = Phosphate_umolL, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location, scales = "free")

mod_Phosphate_umolL<-lm(estimate ~ Phosphate_umolL*Location, data = models )
anova(mod_Phosphate_umolL)

models %>%
  ggplot(aes(x = (NN_umolL)/Phosphate_umolL, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location*Season, scales = "free")

models %>%
  ggplot(aes(x = pH, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location*Season, scales = "free")

models %>%
  ggplot(aes(x = Salinity, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location*Season, scales = "free")

models %>%
  ggplot(aes(x = N_percent, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location*Season, scales = "free")

models %>%
  ggplot(aes(x = del15N, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location*Season, scales = "free")

# run a model of the slope versus NN range
mod_site<-lm(estimate ~ (NN_umolL+pH)*Season, data = models %>% filter(Location == "Varari"))
anova(mod_site)


ggplot(Cdata %>% filter(Location  == "Varari"), aes(x = NEP.proxy, y = NEC.proxy))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~CowTagID, scales = "free")

ggplot(Cdata %>% filter(Location  == "Cabral"), aes(x = NEP.proxy, y = NEC.proxy))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~CowTagID*Season, scales = "free")

### join in the benthic data
source(here("Scripts","BenthicData.R"))

models2 <-models %>%
  left_join(Benthic.Cover_Categories) %>%
  mutate(logratio = log((TotalAlgae+1)/(TotalCalc+1)))

models2 %>%
  ggplot(aes(x = logratio, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location*Season, scales = "free")

models2_both <-models_both %>%
  left_join(Benthic.Cover_Categories) %>%
  mutate(logratio = log((TotalAlgae+1)/(TotalCalc+1)))

models2_both %>%
  ggplot(aes(x = TotalCoral, y = estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_label(aes(label = CowTagID))+
  facet_wrap(~Location, scales = "free")


Cdata %>%
  #filter(NEP<60)%>%
  left_join(turb_all) %>%
  left_join(Benthic.Cover_Categories) %>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = log(NN_umolL), y = NEP.proxy,color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location,scales = "free")

moda<-lm(NEP.proxy~log(NN_umolL)*Season, data = Cdata %>%
#           filter(NEP<60)%>%
           left_join(turb_all) %>%
           left_join(Benthic.Cover_Categories) %>%
           filter(Location == "Varari", Plate_Seep == "Plate"))

anova(moda)


Cdata %>%
 # filter(NEP<60)%>%
  left_join(turb_all) %>%
  left_join(Benthic.Cover_Categories) %>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = Salinity, y = pH, shape = Tide_Time))+
  geom_point(aes(color = Season))+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scale = "free")


Cdata %>%
#  filter(NEP<60)%>%
  left_join(turb_all) %>%
  left_join(Benthic.Cover_Categories) %>%
  filter( Plate_Seep == "Plate") %>%
  ggplot(aes(x = pH, y = NEC.proxy, color = Season)) +
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location,scales = "free")


modb<-lm(NEC.proxy~pH*Season, data = Cdata %>%
#           filter(NEP<60)%>%
           filter(Location == "Varari", Plate_Seep == "Plate"))

anova(modb)

## Calculate summaries
AllDataSummary<- Cdata %>%
  left_join(Benthic.Cover_Categories) %>%
  group_by(Location, Plate_Seep, CowTagID)%>%
  summarise_at(vars(Salinity:Ammonia_umolL, pCO2), .funs =  function(x)(max(x, na.rm = TRUE) - min(x,na.rm = TRUE))) %>%
  left_join(Benthic.Cover_Categories) %>%
  left_join(turb_all) %>%
  mutate(logratio = log((TotalAlgae +1)/(TotalCalc+1)))


AllDataSummary %>%
  filter(Plate_Seep == "Plate",
        # logratio > 0
        ) %>% # one outlier
  ggplot(aes(x = logratio, y = pH))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

modpHratio<-lm(pH~logratio, data = 
                 AllDataSummary %>%
                 filter(Location  == "Varari", 
                        Plate_Seep == "Plate"
                    #    logratio > 0)
               ) )

modpHratioc<-lm(pH~logratio, data = 
                 AllDataSummary %>%
                 filter(Location  == "Cabral", 
                        Plate_Seep == "Plate"
                        #    logratio > 0)
                 ) )
#removed the one outlier 
anova(modpHratioc)

## boxplots of NEC and NEP
Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  droplevels()%>%
ggplot(aes(x = Tide, y = NEP.proxy, color = Day_Night))+
  geom_boxplot()+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  droplevels()%>%
  ggplot(aes(x = Tide, y = NEC.proxy, color = Day_Night))+
  geom_boxplot()+
  facet_wrap(~Location*Season, scales = "free")

# average rate at each location by mean or cov in silicate

CoV <- function(x){
  var(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
}

testing<-Cdata %>%
  filter(Plate_Seep == "Plate", NEC < 5) %>%
  group_by(Location,  CowTagID, Season)%>%
  summarise_at(vars(Salinity:Ammonia_umolL, pCO2, NEC, NEP,NEC.proxy, NEP.proxy, MarineHumic_Like:Lignin_Like), .funs =  c("mean" = mean, "CoV" = CoV, "sum" = sum)) %>%
  droplevels()%>%
  left_join(turb_all) %>%
  left_join(Benthic.Cover_Categories)
  
testing %>%  
  ggplot(aes(x = pH_mean, y = NEC.proxy_mean, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")

testing %>%  
ggplot(aes(x = log(Silicate_umolL_CoV), y = NEC.proxy_mean, color = Season))+
   geom_point()+
  geom_smooth(method = "lm")+
   facet_wrap(~Season*Location, scales = "free")

testing %>%  
  ggplot(aes(x = log(Silicate_umolL_mean), y = MarineHumic_Like_mean, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")

testing %>%  
  ggplot(aes(y = TotalCoral, x = Silicate_umolL_mean, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")

testing %>%  
  ggplot(aes(x = TotalCoral, y = NEC.proxy_mean, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")

testing %>%  
  ggplot(aes(x = Silicate_umolL_mean, y = NEC.proxy_mean, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")

testing %>%  
  ggplot(aes(y = NN_umolL_mean, x = Silicate_umolL_mean, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")


testing %>%  
  ggplot(aes(x = log(NN_umolL_mean), y = NEP_sum, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")

testing %>%  
  ggplot(aes(x = log(NN_umolL_mean), y = NEC.proxy_mean, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")

moda<-lm(NEC.proxy_mean ~ pH_mean*Season, data = testing%>% filter(Location == "Varari"))
anova(moda)  

moda<-lm(NEC.proxy_mean ~ log(Silicate_umolL_CoV)*Season, data = testing%>% filter(Location == "Varari"))
anova(moda) 
 
### look at fDOM
Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate",
         VisibleHumidic_Like <0.05) %>%
  ggplot(aes(y = MarineHumic_Like+ VisibleHumidic_Like, x = NEP.proxy, color = Day_Night, shape = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Net community production")+
  facet_wrap(~Tide*Location, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate",
         VisibleHumidic_Like <0.05) %>%
  ggplot(aes(y = Tryptophan_Like+Tyrosine_Like, x = NEC.proxy, color = Day_Night, shape = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Net community calcification")+
  facet_wrap(~Tide*Location, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate",
         VisibleHumidic_Like <0.05) %>%
  ggplot(aes(x = Day_Night, y = Tryptophan_Like+Tyrosine_Like, fill = Tide))+
  geom_boxplot()+
  facet_wrap(~Season*Location, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate",
         VisibleHumidic_Like <0.05) %>%
  ggplot(aes(x = Day_Night, y = MarineHumic_Like+VisibleHumidic_Like, fill = Tide))+
  geom_boxplot()+
  facet_wrap(~Season*Location, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate",
         VisibleHumidic_Like <0.05) %>%
  ggplot(aes(y = MarineHumic_Like, x =Silicate_umolL, color = Day_Night))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(y = "Humics")+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  select(Location, Tide, Silicate_umolL, TA, pH, NN_umolL, Phosphate_umolL, Ammonia_umolL) %>%
  pivot_longer(Silicate_umolL:Ammonia_umolL) %>%
  ggplot(aes(x = value, fill = Location)) +
  geom_density(alpha = 0.5)+
  scale_x_continuous(trans = "log10")+
  theme_bw() +
  facet_wrap(~name, scale = "free")

ggsave(here("Output", "desnityplotssite.png"), width = 8, height = 4)

Cdata %>% anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = log(Silicate_umolL), y = Tryptophan_Like, color = Day_Night))+ 
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>% anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = log(Silicate_umolL), y = VisibleHumidic_Like+ MarineHumic_Like, color = Day_Night))+ 
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>% anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = log(Silicate_umolL), y = Tryptophan_Like+Tyrosine_Like, color = Day_Night))+ 
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>% anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Seep"|Plate_Seep == "Spring" ) %>%
  ggplot(aes(x = Silicate_umolL, y = MarineHumic_Like+VisibleHumidic_Like, color = Season))+ 
  geom_point() +
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  geom_smooth(method = "lm")+
  facet_wrap(~Location)

testing %>%  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  ggplot(aes(x = TotalCoral, y = Tyrosine_Like_mean, color = log(NN_umolL_mean)))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")

#### Buffering models

ModelData<-Cdata %>%
 anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  select(Location, Tide, Day_Night, TimeBlock, Season, pH, Salinity, Silicate_umolL, NEP, NEP.proxy, NEC, NEC.proxy) %>%
  mutate(Silicate_umolL = log(Silicate_umolL)) %>% 
  mutate_at(vars("pH":"NEC.proxy"),  function(x) scale(x, center = TRUE)[,1]) 

## Effect of pH on NEC
Cdata %>% anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = pH, y = NEC.proxy, color = Season, shape = Tide))+ 
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Tide, scales = "free")

mod.pH<-lm(NEC.proxy~pH*Season*Tide, data = ModelData %>%
             filter(Location == "Varari"))
anova(mod.pH)


v1<-summary(margins(mod.pH, at = list(Season = c("Dry","Wet"), Tide = c("High","Low"))))[1:4,] %>%
  tibble()%>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Tide, color = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(0.3, 1), guide = "none")+
  labs(#title = "Varari",
       y = " ",
       x = " ")+
  theme_bw()

mod.pHC<-lm(NEC.proxy~pH*Season*Tide, data = ModelData %>%
              filter(Location == "Cabral"))
anova(mod.pHC)

c1<-summary(margins(mod.pHC, at = list(Season = c("Dry","Wet"), Tide = c("High","Low"))))[1:4,] %>%
  tibble()%>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Tide, color = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(0.3, 1), guide = "none")+
  labs(#title = "Cabral",
       y = "Effect of pH on NEC",
       x = " ")+
  theme_bw()

# TA distribution
v2<-Cdata %>%
  filter(Plate_Seep == "Plate", Location == "Varari") %>%
  ggplot()+
  geom_density_ridges(aes(y = Tide, x = TA, fill = Season), alpha = 0.5, show.legend = FALSE, scale = 1)+
  scale_x_continuous(limits = c(2000,2500))+
  labs(y = "Varari", x = " ")+
  theme_bw()
  

c2<-Cdata %>%
  filter(Plate_Seep == "Plate", Location == "Cabral") %>%
  ggplot()+
  geom_density_ridges(aes(y = Tide, x = TA, fill = Season), alpha = 0.5, show.legend = FALSE, scale = 1)+
  scale_x_continuous(limits = c(2000,2500))+
  labs(y = "Cabral")+
  theme_bw()

## effect of NEP and Salinity on pH

mod.NEP<-lm(pH~(Silicate_umolL+NEP.proxy)*Season*Tide, data = ModelData %>%
              filter(Location == "Varari"))
anova(mod.NEP)

v3<-summary(margins(mod.NEP, at = list(Season = c("Dry","Wet"), Tide = c("High","Low"))))[1:4,] %>%
  tibble() %>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Tide, color = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(1, 0.3), guide = "none")+
  labs(#title = "Cabral",
    y = " ", x = " ")+
  theme_bw()

v4<-summary(margins(mod.NEP, at = list(Season = c("Dry","Wet"), Tide = c("High","Low"))))[9:12,] %>%
  tibble() %>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Tide, color = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(0.3, 1), guide = "none")+
  labs(#title = "Cabral",
    y = " ", x = "")+
  theme_bw()

mod.NEPc<-lm(pH~(log(Silicate_umolL)+NEP.proxy)*Season*Tide, data = ModelData %>%
               filter(Location == "Cabral"))
anova(mod.NEPc)

c3<-summary(margins(mod.NEPc, at = list(Season = c("Dry","Wet"), Tide = c("High","Low"))))[1:4,] %>%
  tibble() %>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Tide, color = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(1, 0.3), guide = "none")+
  labs(#title = "Cabral",
    y = "Effect of NEP on pH",
    x = " ")+
  theme_bw()

c4<-summary(margins(mod.NEPc, at = list(Season = c("Dry","Wet"), Tide = c("High","Low"))))[9:12,] %>%
  tibble() %>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Tide, color = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(0.3, 1), guide = "none")+
  labs(#title = "Cabral",
    y = "Effect of SGD on pH", x = " ")+
  theme_bw()

((v2+v4+v3+v1)+ plot_layout(ncol = 4))/((c2+c4+c3+c1) + plot_layout(ncol = 4))+plot_layout(guides = "collect")

### not divided by tide 
Cdata %>% anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = pH, y = NEC.proxy, color = Season))+ 
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

mod.pH<-lm(NEC.proxy~pH*Season, data = ModelData %>%
             filter(Location == "Varari"))
anova(mod.pH)


v1<-summary(margins(mod.pH, at = list(Season = c("Dry","Wet"))))[1:2,] %>%
  tibble()%>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(1, 0.3), guide = "none")+
  labs(#title = "Varari",
    y = " ",
    x = " ")+
  theme_bw()

mod.pHC<-lm(NEC.proxy~pH*Season, data = ModelData %>%
              filter(Location == "Cabral"))
anova(mod.pHC)

c1<-summary(margins(mod.pHC, at = list(Season = c("Dry","Wet"))))[1:2,] %>%
  tibble()%>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(0.3, 1), guide = "none")+
  labs(#title = "Cabral",
    y = "Effect of pH on NEC",
    x = " ")+
  theme_bw()

# TA distribution
v2<-Cdata %>%
  filter(Plate_Seep == "Plate", Location == "Varari") %>%
  ggplot()+
  geom_density_ridges(aes(y = Season, x = TA, fill = Tide), alpha = 0.5,  scale = 1)+
  scale_x_continuous(limits = c(2000,2500))+
  labs(y = "Varari", x = " ")+
  theme_bw()


c2<-Cdata %>%
  filter(Plate_Seep == "Plate", Location == "Cabral") %>%
  ggplot()+
  geom_density_ridges(aes(y = Season, x = TA, fill = Tide), alpha = 0.5,  scale = 1)+
  scale_x_continuous(limits = c(2000,2500))+
  labs(y = "Cabral")+
  theme_bw()

## effect of NEP and Salinity on pH

mod.NEP<-lm(pH~(Silicate_umolL+NEP.proxy)*Season, data = ModelData %>%
              filter(Location == "Varari"))
anova(mod.NEP)

v3<-summary(margins(mod.NEP, at = list(Season = c("Dry","Wet"))))[1:2,] %>%
  tibble() %>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x =  Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(1, 0.3), guide = "none")+
  labs(#title = "Cabral",
    y = " ", x = " ")+
  theme_bw()

v4<-summary(margins(mod.NEP, at = list(Season = c("Dry","Wet"))))[5:6,] %>%
  tibble() %>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c( 1, 0.3), guide = "none")+
  labs(#title = "Cabral",
    y = " ", x = "")+
  theme_bw()

mod.NEPc<-lm(pH~(Silicate_umolL+NEP.proxy)*Season, data = ModelData %>%
               filter(Location == "Cabral"))
anova(mod.NEPc)

c3<-summary(margins(mod.NEPc, at = list(Season = c("Dry","Wet"))))[1:2,] %>%
  tibble() %>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x =  Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(1, 0.3), guide = "none")+
  labs(#title = "Cabral",
    y = "Effect of NEP on pH",
    x = " ")+
  theme_bw()

c4<-summary(margins(mod.NEPc, at = list(Season = c("Dry","Wet"))))[5:6,] %>%
  tibble() %>%
  mutate(sig = ifelse(p <= 0.05, "yes", "no")) %>%
  ggplot(aes(x = Season, y = AME, alpha = sig))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  scale_y_continuous(limits = c(-1,1.5))+
  coord_flip()+
  scale_alpha_manual(values = c(0.3, 1), guide = "none")+
  labs(#title = "Cabral",
    y = "Effect of SGD on pH", x = " ")+
  theme_bw()

((v2+v4+v3+v1)+ plot_layout(ncol = 4))/((c2+c4+c3+c1) + plot_layout(ncol = 4))+plot_layout(guides = "collect")
