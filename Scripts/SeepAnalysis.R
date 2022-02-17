### Look for patterns among data collected at the seep versus WW3 data
## Created by Nyssa Silbiger
## Edited on 10/4/2021

## load libraries ####-----------
library(here)
library(tidyverse)
library(lubridate)
library(forecast)
library(patchwork)
library(viridis)

## Read in the different datasets

## Tide predictions------------
tides<-read_tsv(here("Data","IslandData","TidePredictions.txt"), skip = 13) %>%
  mutate(date = ymd_hms(paste(Date,Time))) %>%
  select(date,tideheight = Pred)

## Weather data (wind, rain, waves from windguru)#####
weather<-read_csv(here("Data","IslandData","weather.csv")) %>%
  left_join(tides)

# CT----------------------------
CondPath<-here("Data", "Varari", "CT")
files <- dir(path = CondPath,pattern = ".csv", full.names = TRUE)

CT_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,TempInSitu, Salinity_psu)%>%
  mutate(Site = "Varari")

# Cabral
CondPath<-here("Data", "Cabral", "CT")
files <- dir(path = CondPath,pattern = ".csv", full.names = TRUE)

CT_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,TempInSitu, Salinity_psu)%>%
  mutate(Site = "Cabral")


#Water Level-----
WLPath<-here("Data", "Varari", "WL")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

WL_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date, Depth)%>%
  mutate(Site = "Varari")

WLPath<-here("Data", "Cabral", "WL")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

#Cabral
WL_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date, Depth)%>%
  mutate(Site = "Cabral")


## pH ----------------
pHPath<-here("Data", "Varari", "pH")
files <- dir(path = pHPath,pattern = ".csv", full.names = TRUE)

pH_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,pH = pH_total)%>%
  mutate(Site = "Varari")

#Cabral
pHPath<-here("Data", "Cabral", "pH")
files <- dir(path = pHPath,pattern = ".csv", full.names = TRUE)

pH_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,pH = pH_total)%>%
  mutate(Site = "Cabral")


## Lux ----------------
LUXPath<-here("Data", "Varari", "LUX")
files <- dir(path = LUXPath,pattern = ".csv", full.names = TRUE)

LUX_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,Lux)%>%
  mutate(Site = "Varari")

#Cabral
LUXPath<-here("Data", "Cabral", "LUX")
files <- dir(path = LUXPath,pattern = ".csv", full.names = TRUE)

LUX_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,Lux)%>%
  mutate(Site = "Cabral")

## PAR -------------------
PARPath<-here("Data", "Varari", "PAR")
files <- dir(path = PARPath,pattern = ".csv", full.names = TRUE)

PAR_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,PAR)%>%
  mutate(Site = "Varari")

#Cabral
PARPath<-here("Data", "Cabral", "PAR")
files <- dir(path = PARPath,pattern = ".csv", full.names = TRUE)

PAR_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,PAR)%>%
  mutate(Site = "Cabral")


## DO ---------------
DOPath<-here("Data", "Varari", "DO")
files <- dir(path = DOPath,pattern = ".csv", full.names = TRUE)

DO_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,DO_mg_L)%>%
  mutate(Site = "Varari")

#Cabral
DOPath<-here("Data", "Cabral", "DO")
files <- dir(path = DOPath,pattern = ".csv", full.names = TRUE)

DO_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,DO_mg_L)%>%
  mutate(Site = "Cabral")

### Join everything together

AllVarari<-CT_Varari %>%
  left_join(WL_Varari)%>%
  left_join(pH_Varari) %>%
  left_join(LUX_Varari) %>%
  left_join(PAR_Varari)%>%
  left_join(DO_Varari)%>%
  relocate(Site, .before = date) %>% # move the site column
  mutate(PAR_calc = ifelse(is.na(PAR), # if PAR is missing, calculate it from LuX or else leave it the same
                           16778.33+(-0.5003277-16778.33)*exp(1)^(-exp( -13.29572)*Lux),PAR)) %>%
  select(-Lux,-PAR) # remove Lux and original PAR

Allcabral<-CT_Cabral %>%
  left_join(WL_Cabral)%>%
  left_join(pH_Cabral) %>%
  left_join(LUX_Cabral) %>%
  left_join(PAR_Cabral)%>%
  left_join(DO_Cabral)%>%
  relocate(Site, .before = date) %>% # move the site column
  mutate(PAR_calc = ifelse(is.na(PAR), # if PAR is missing, calculate it from LuX or else leave it the same
                           16778.33+(-0.5003277-16778.33)*exp(1)^(-exp( -13.29572)*Lux),PAR)) %>%
  select(-Lux,-PAR) # remove Lux and original PAR


#### fill in the missing times with NA for easier plotting

AllVarari<-data.frame(date = seq(AllVarari$date[1], AllVarari$date[nrow(AllVarari)], by = "1 min"))  %>%
  full_join(AllVarari, by = "date")
  
AllCabral<-data.frame(date = seq(Allcabral$date[1], Allcabral$date[nrow(Allcabral)], by = "1 min"))  %>%
  full_join(Allcabral, by = "date")


## Add in sampling times for Varari and Cabral
Varari_sample<-tibble(datetime = ymd_hms(c("2021-08-05 11:57:00", "2021-08-05 00:00:00",
                         "2021-08-08 18:30:00", "2021-08-06 06:40:00",
                         "2021-08-08 07:30:00" )))
Cabral_sample<-tibble(datetime =ymd_hms(c("2021-08-09 07:00:00", "2021-08-09 13:00:00",
                         "2021-08-09 01:10:00","2021-08-09 19:00:00",
                         "2021-08-10 07:00:00")))

## Simple plot

AllVarari %>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
#  geom_vline(data = Varari_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw() +
  labs(title = "Varari Sled")
ggsave(here('Output',"Varari_timeseries.pdf"), width = 8, height = 5)

# Just during the sampling times
AllVarari %>%
  filter(date >= ymd("2021-08-05"), date <= ymd("2021-08-09"))%>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  geom_vline(data = Varari_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw() +
  labs(title = "Varari Sled")

# Cabral
AllCabral %>%
  filter(date >= ymd("2021-08-09"), date <= ymd("2021-08-10"))%>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  geom_vline(data = Cabral_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw() +
  labs(title = "Cabral Sled")

# Just during sampling times
AllCabral %>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  geom_vline(data = Cabral_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw() +
  labs(title = "Cabral Sled")


# take hour averages

AllVarari_onehour<-AllVarari %>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  mutate(date = floor_date(date,"hour"))%>% # round to the lowest hour
  group_by(Site,Params, date)%>%
  summarise(Values = mean(Values,na.rm = TRUE))%>% # take the hourly average
  ungroup() %>%
  pivot_wider(names_from = Params, values_from = Values) %>%
  left_join(weather) # join in the weather data

# plot with weather
AllVarari_onehour %>%
  pivot_longer(cols = Depth:tideheight, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw()
ggsave(here("Output","Varari_ts_hourly.pdf"), width = 8)

AllVarari_onehour<-AllVarari %>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  mutate(date = floor_date(date,"hour"))%>% # round to the lowest hour
  group_by(Site,Params, date)%>%
  summarise(Values = mean(Values,na.rm = TRUE))%>% # take the hourly average
  ungroup() %>%
  pivot_wider(names_from = Params, values_from = Values) %>%
  left_join(weather) # join in the weather data

AllCabral_onehour<-AllCabral %>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  mutate(date = floor_date(date,"hour"))%>% # round to the lowest hour
  group_by(Site,Params, date)%>%
  summarise(Values = mean(Values,na.rm = TRUE))%>% # take the hourly average
  ungroup() %>%
  pivot_wider(names_from = Params, values_from = Values) %>%
  left_join(weather) # join in the weather data

# plot with weather
AllCabral_onehour %>%
  pivot_longer(cols = Depth:tideheight, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw()
ggsave(here("Output","Cabral_ts_hourly.pdf"), width = 8)


## plot waves vs depth and tide vs depth
WD_TP<-AllVarari_onehour %>%
  bind_rows(AllCabral_onehour)%>%
  drop_na(Site)%>%
  ggplot(aes(x = tideheight, y = Depth, color = waves))+
  geom_point()+
  geom_smooth()+
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75,1,1.25))+
  ylim(0,1.5)+
  labs(x = 'Tide Predictions (m)',
       y = 'Water Depth (m)')+
  theme_bw()+ 
  guides(colour=guide_colourbar(barwidth=15,label.position = "top", 
                                            title.position = "top",
                                            direction = "horizontal"))+
  scale_color_viridis_c("Wave Height (m)")+
  theme(legend.position="top")+
  facet_wrap(~Site, ncol = 1)

WD_Wave<-AllVarari_onehour %>%
  bind_rows(AllCabral_onehour)%>%
  drop_na(Site)%>%
  ggplot(aes(x = waves, y = Depth, color = tideheight))+
  geom_point()+
  geom_smooth()+
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75,1,1.25))+
  ylim(0,1.5)+
  labs(x = 'Sig Wave height (m)',
       y = "")+
       #y = 'Water Depth (m)')+
  guides(colour=guide_colourbar(barwidth=15,label.position = "top", 
                                title.position = "top",
                                direction = "horizontal"))+
  theme_bw()+
  theme( axis.text.y = element_blank(),
         legend.position = 'top')+
  scale_color_viridis("Tide Predictions (m)",option = "plasma")+
  facet_wrap(~Site, ncol = 1)
  

WD_TP + WD_Wave
ggsave(filename = here("Output","WaterDepth_tide_wave.pdf"), width = 8)

### Depth versus different parameters
AllVarari_onehour %>%
  pivot_longer(cols = c("pH","Salinity_psu","TempInSitu", "DO_mg_L"), names_to = "Params", values_to = "Values")%>%
  ggplot(aes(x = Depth, y  = Values, color = waves))+
  geom_point()+
  facet_wrap(~Params, scales="free_y")+
  theme_bw()+
  scale_color_viridis_c(option = "plasma")


ggplot(AllVarari, aes(x = pH, y = DO_mg_L, col = TempInSitu))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_c(option = "plasma")


### Plot pH data hand collected from seep on pH HOBO data ####
DiscreteData<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv")


AllVarari %>%
  filter(date >= ymd("2021-08-05"), date <= ymd("2021-08-09"))%>%
  ggplot()+
  geom_line(aes(x = date, y = pH))+
  geom_point(data = DiscreteData%>% filter(Plate_Seep =="Seep", Location=="Varari"), aes(x = DateTime, y = pH), color = "red", size =2)+
  theme_bw() +
  labs(title = "Varari Sled")


AllCabral %>%
  filter(date >= ymd("2021-08-08"), date <= ymd("2021-08-10"))%>%
   ggplot(aes(x = date, y = pH))+
  geom_line()+
  geom_point(data = DiscreteData%>% filter(Plate_Seep =="Seep", Location=="Cabral"), aes(x = DateTime, y = pH), color = "red", size =2)+
  theme_bw() +
  labs(title = "Cabral Sled")
