#### Look at long term trends in seep data ###
#### By Nyssa Silbiger ###
#### Created on 12/15/22 ####
############################




## load libraries ####-----------
library(here)
library(tidyverse)
library(lubridate)
library(forecast)
library(patchwork)
library(viridis)
library(janitor)
library(ggh4x)
library(ggtext)
library(patchwork)


### read in all the data
CondPath<-here("Data", "Varari", "CT")
files <- dir(path = CondPath,pattern = ".csv", full.names = TRUE)

CT_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  mutate(Site = "Varari",
         date = if_else(is.na(date),Date,date)) %>% # dealing with different outputs
  select(date,TempInSitu, Salinity_psu)

WLPath<-here("Data", "Varari", "WL")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

WL_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date, Depth) %>%
  mutate(Site = "Varari")

pHPath<-here("Data", "Varari", "pH")
files <- dir(path = pHPath,pattern = ".csv", full.names = TRUE)

pH_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,pH = pH_total)%>%
  mutate(Site = "Varari")

LUXPath<-here("Data", "Varari", "LUX")
files <- dir(path = LUXPath,pattern = ".csv", full.names = TRUE)

LUX_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,Lux)%>%
  mutate(Site = "Varari")

PARPath<-here("Data", "Varari", "PAR")
files <- dir(path = PARPath,pattern = ".csv", full.names = TRUE)

PAR_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,PAR)%>%
  mutate(Site = "Varari") %>%
  mutate(date = ceiling_date(date, unit = "minutes")) # the par data isn't always on the minute mark

DOPath<-here("Data", "Varari", "DO")
files <- dir(path = DOPath,pattern = ".csv", full.names = TRUE)

DO_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,DO_mg_L)%>%
  mutate(Site = "Varari")

waves<-read_csv(here("Data","IslandData","WestSideADCPMCR.csv")) %>%
  filter(datetime > ymd("2021-08-03")) %>%
  mutate(date = datetime)%>%
  select(!datetime)

# bring everything together
## missing data in join.. check the june data

AllVarari<-CT_Varari %>%
  full_join(WL_Varari)%>%
  full_join(pH_Varari) %>%
  full_join(LUX_Varari) %>%
  full_join(PAR_Varari)%>%
  full_join(DO_Varari)%>%
  full_join(waves)%>% # add in wave height data from LTER 6 offshore
  relocate(Site, .before = date)  # move the site column
  # mutate(PAR_calc = case_when(is.na(PAR) & Season == "Dry" ~ 16778.33+(0-16778.33)*exp(1)^(-exp( -13.29572)*Lux), # calculate PAR when it is missing
  #                             !is.na(PAR)~PAR,
  #                             is.na(PAR) & Season =="Wet" ~ 1555.751648 +(0-1555.751648 )*exp(1)^(-exp( -10.777624)*Lux)))  %>%
 # select(-Lux,-PAR) %>%  # remove Lux and original PAR
#  mutate(PAR_calc = ifelse(PAR_calc<0,0,PAR_calc)) # if PAR is negative make it 0


## make some plots
AllVarari %>%
  ggplot(aes(x = date, y = Salinity_psu))+
  geom_line()

AllVarari %>%
 # drop_na(Significant_wave_height, Depth) %>%
  ggplot(aes(x = Depth, y = Salinity_psu, color = log(Lux+1)))+
  geom_point()


AllVarari %>%
  filter(Depth < 1.5) %>%
  # drop_na(Significant_wave_height, Depth) %>%
  ggplot(aes(x = Depth, y = Salinity_psu, color = log(Lux+1)))+
  geom_point()


AllVarari %>%
   drop_na(pH, Depth) %>%
  ggplot(aes(x = Depth, y = pH))+
  geom_point()


AllVarari %>%
  drop_na(pH, Salinity_psu) %>%
  ggplot(aes(x = Salinity_psu, y = TempInSitu))+
  geom_point()

### Need to line up March waveheight data

AllVarari %>% ## values are not lining up
  filter(date > ymd_hms("2022-06-01 00:00:00")) %>%
  pivot_longer(cols = TempInSitu:Significant_wave_height, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  #  geom_vline(data = Varari_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params, scales = "free", ncol = 2)+

  theme_bw() +
  labs(title = "Varari Sled")

