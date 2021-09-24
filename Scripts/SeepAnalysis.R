### Look for patterns among data collected at the seep versus WW3 data
## Created by Nyssa Silbiger
## Edited on 9/23/2021

## load libraries ####-----------
library(here)
library(tidyverse)
library(lubridate)
library(forecast)

## Read in the different datasets

## Weather data (wind, rain, waves from windguru)#####
weather<-read_csv(here("Data","IslandData","weather.csv"))

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

### Join everything together

AllVarari<-CT_Varari %>%
  left_join(WL_Varari)%>%
  left_join(pH_Varari) %>%
  left_join(LUX_Varari) %>%
  relocate(Site, .before = date) # move the site column

#### fill in the missing times with NA for easier plotting

AllVarari<-data.frame(date = seq(AllVarari$date[1], AllVarari$date[nrow(AllVarari)], by = "1 min"))  %>%
  full_join(AllVarari, by = "date")
  

## Simple plot

AllVarari %>%
  pivot_longer(cols = TempInSitu:Lux, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw() +
  labs(title = "Varari Sled")
ggsave(here('Output',"Varari_timeseries.pdf"), width = 6)

# take hour averages

AllVarari_onehour<-AllVarari %>%
  pivot_longer(cols = TempInSitu:Lux, names_to = "Params", values_to = "Values") %>%
  mutate(date = floor_date(date,"hour"))%>% # round to the lowest hour
  group_by(Site,Params, date)%>%
  summarise(Values = mean(Values,na.rm = TRUE))%>% # take the hourly average
  ungroup() %>%
  pivot_wider(names_from = Params, values_from = Values) %>%
  left_join(weather) # join in the weather data

# plot with weather
AllVarari_onehour %>%
  pivot_longer(cols = Depth:waves, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw()
ggsave(here("Output","Varari_ts_hourly.pdf"), width = 6)
  