# playing around

### The august data here still has the two low day values

library(tidyverse)
library(seacarb)
library(broom)
library(lubridate)
library(here)

# This is temporary until we get the rest of the data
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
mod15N<-lm(estimate~del15N*Location, data = models2 )
anova(mod15N)

AllData %>%
  mutate(TA_sal = TA*Salinity/36, # salinity normalize
         DIC_sal = DIC*Salinity/36) %>%
  ggplot(aes(x = TA_sal, y = DIC_sal))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*CowTagID, scale = "free")

# read in the benthic data and look for patterns between TA/DIC

source(here("Scripts","BenthicData.R"))
# join the data
AllData<-AllData %>%
left_join(Benthic.Cover_Categories)

models2 <-models2 %>%
  left_join(Benthic.Cover_Categories)

models2 %>%
ggplot(aes(x = log((TotalCalc+1)/(TotalAlgae+1)), y = estimate, color = Location))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location)

BenthicTADIC<-lm(estimate ~ log((TotalCalc+1)/(TotalAlgae+1))*Location, data = models2)
anova(BenthicTADIC)

## Turn the corals and algae to icons
p1<-models2 %>%
  ggplot(aes(y = log((TotalCalc+1)/(TotalAlgae+1)), x = del15N))+
  geom_point(aes(color = Location, size = N_percent))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, lty = 2)+
  scale_size_binned("%N (Nutrient Loading)") +
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = "log ratio of calcifiers to fleshy algae",
       color = "",
       x = "del15N (Nutrient Source)")+
  annotate("text", x = 4.75, y = 0.5, label = "Calcifier-dominated")+
  annotate("text", x = 4.75, y = -0.5, label = "Fleshy algal-dominated")+
  theme_bw()+
  theme(legend.direction = "horizontal",
        legend.position = c(.18, .1))

ggsave(here("Output","Benthic15N.pdf"), width = 10, height = 10)

p2<-models2 %>%
  ggplot(aes(y = log((TotalCalc+1)/(TotalAlgae+1)), x = N_percent, color = Location))+
  geom_point()+
  geom_smooth(method = "lm", data = models2%>%filter(Location == "Varari")) +
  geom_hline(yintercept = 0, lty = 2)+
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = "log ratio of calcifiers to fleshy algae",
       color = "",
       x = "%N (Nutrient Loading)")+
 # annotate("text", x = 1.1, y = 0.5, label = "Calcifier-dominated")+
 # annotate("text", x = 1.1, y = -0.5, label = "Fleshy algal-dominated")+
  theme_bw()+
  theme(legend.position = "none")


modBenthic15N<-lm(log((TotalCalc+1)/(TotalAlgae+1))~del15N*Location, data = models2)
anova(modBenthic15N)

modBenthicNpercent<-lm(log((TotalCalc+1)/(TotalAlgae+1))~N_percent*Location, data = models2)
anova(modBenthicNpercent)

modBenthicNpercent<-lm(log((TotalCalc+1)/(TotalAlgae+1))~N_percent, data = models2 %>% filter(Location == "Varari"))
anova(modBenthicNpercent)

p1+p2


####3 Calculate summaries of all Data 
AllDataSummary<- AllData %>%
  group_by(Location, Plate_Seep, CowTagID)%>%
  summarise_at(vars(Salinity:Ammonia_umolL, pCO2), .funs =  function(x)(max(x, na.rm = TRUE) - min(x,na.rm = TRUE))) %>%
  left_join(Benthic.Cover_Categories) %>%
  mutate(logratio = log((TotalCalc +1)/(TotalAlgae+1)))
  

# Varari only for hendrikje of log community versus pH range
AllDataSummary %>%
  filter(Plate_Seep == "Plate", Location == "Varari") %>%
  ggplot(aes(x=logratio, y = pH ))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(y = "pH Range",
       x = "log(Calcifiers/Algae) of Benthic Data",
       title = "Varari")+
  geom_vline(xintercept = 0, lty = 2)+
  annotate("text", x = -1.25, y = 0.2, label = "More Algae-Dominated")+
  annotate("segment", x = -2, xend = -3, y = .2, yend = .2,
           arrow = arrow( angle = 45, length = unit(.2,"cm")))+
 # annotate("text", x = 0.3, y = 0.3, label = "Calcifier-Dominated")+
  theme_bw()

ggsave(here("Output", "CommunityvspH.pdf"), width = 8, height = 8)

#Model of community versus pH range
modpHlog<-lm(pH ~ logratio, data = AllDataSummary %>% filter(Plate_Seep == "Plate" & Location == "Varari"))
anova(modpHlog)


AllDataSummary %>%
  filter(Plate_Seep == "Plate") %>%
ggplot(aes(x = pH, y = TA, color = Location))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")

#Model of TA versus pH range
modpHTA<-lm(TA ~ pH*Location, data = AllDataSummary)
anova(modpHTA)


AllDataSummary %>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = logratio, y = TA, color = Location))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location, scales = "free")
