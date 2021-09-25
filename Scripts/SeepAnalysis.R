### Look for patterns among data collected at the seep versus WW3 data
## Created by Nyssa Silbiger
## Edited on 9/23/2021

## load libraries ####-----------
library(here)
library(tidyverse)
library(lubridate)
library(forecast)
library(patchwork)

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

#Water Level-----
WLPath<-here("Data", "Varari", "WL")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

WL_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date, Depth)%>%
  mutate(Site = "Varari")

## pH ----------------
pHPath<-here("Data", "Varari", "pH")
files <- dir(path = pHPath,pattern = ".csv", full.names = TRUE)

pH_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,pH)%>%
  mutate(Site = "Varari")

## Lux ----------------
LUXPath<-here("Data", "Varari", "LUX")
files <- dir(path = LUXPath,pattern = ".csv", full.names = TRUE)

LUX_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,Lux)%>%
  mutate(Site = "Varari")

## PAR -------------------
PARPath<-here("Data", "Varari", "PAR")
files <- dir(path = PARPath,pattern = ".csv", full.names = TRUE)

PAR_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,PAR)%>%
  mutate(Site = "Varari")


### Join everything together

AllVarari<-CT_Varari %>%
  left_join(WL_Varari)%>%
  left_join(pH_Varari) %>%
  left_join(LUX_Varari) %>%
  left_join(PAR_Varari)%>%
  relocate(Site, .before = date) %>% # move the site column
  mutate(PAR_calc = ifelse(is.na(PAR), # if PAR is missing, calculate it from LuX or else leave it the same
                           16778.33+(-0.5003277-16778.33)*exp(1)^(-exp( -13.29572)*Lux),PAR)) %>%
  select(-Lux,-PAR) # remove Lux and original PAR

#### fill in the missing times with NA for easier plotting

AllVarari<-data.frame(date = seq(AllVarari$date[1], AllVarari$date[nrow(AllVarari)], by = "1 min"))  %>%
  full_join(AllVarari, by = "date")
  

## Simple plot

AllVarari %>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw() +
  labs(title = "Varari Sled")
ggsave(here('Output',"Varari_timeseries.pdf"), width = 6)

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
ggsave(here("Output","Varari_ts_hourly.pdf"), width = 6)


## plot waves vs depth and tide vs depth
WD_TP<-AllVarari_onehour %>%
  ggplot(aes(x = tideheight, y = Depth, color = waves))+
  geom_point()+
  geom_smooth()+
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75,1,1.25))+
  labs(x = 'Tide Predictions (m)',
       y = 'Water Depth (m)')+
  theme_bw()+ 
  guides(colour=guide_colourbar(barwidth=15,label.position = "top", 
                                            title.position = "top",
                                            direction = "horizontal"))+
  scale_color_viridis_c("Wave Height (m)")+
  theme(legend.position="top")

WD_Wave<-AllVarari_onehour %>%
  ggplot(aes(x = waves, y = Depth, color = tideheight))+
  geom_point()+
  geom_smooth()+
  scale_y_continuous(breaks = c(0, 0.25,0.5,0.75,1,1.25))+
  labs(x = 'Sig Wave height (m)',
       y = "")+
       #y = 'Water Depth (m)')+
  guides(colour=guide_colourbar(barwidth=15,label.position = "top", 
                                title.position = "top",
                                direction = "horizontal"))+
  theme_bw()+
  theme( axis.text.y = element_blank(),
         legend.position = 'top')+
  scale_color_viridis("Tide Predictions (m)",option = "plasma")
  

WD_TP + WD_Wave
ggsave(filename = here("Output","WaterDepth_tide_wave.pdf"), width = 10)

### Depth versus different parameters
AllVarari_onehour %>%
  pivot_longer(cols = c("pH","Salinity_psu","TempInSitu"), names_to = "Params", values_to = "Values")%>%
  ggplot(aes(x = Depth, y  = Values, color = waves))+
  geom_point()+
  facet_wrap(~Params, scales="free_y")+
  theme_bw()+
  scale_color_viridis_c(option = "plasma")
