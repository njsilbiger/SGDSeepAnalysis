#### All the sled data ####
### By Nyssa Silbiger ###
### 2023-08-20 ####

##load libraries ###

library(here)
library(tidyverse)
library(lubridate)

###

## some of the files are character dates and some are UTC... read in separately 
CondPath<-here("Data", "AllSled", "CT","Date_chr")
files <- dir(path = CondPath,pattern = ".csv", full.names = TRUE)

CT_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(Date,TempInSitu, Salinity_psu) %>%
  mutate(Date = mdy_hm(Date))

# Not the UTC files
CondPath<-here("Data", "AllSled", "CT","Date_UTC")
files <- dir(path = CondPath,pattern = ".csv", full.names = TRUE)

CT_Varari2<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(Date,TempInSitu, Salinity_psu) 

CT_Varari<-CT_Varari %>%
  bind_rows(CT_Varari2)%>%
  arrange(Date)

### Water level ###
#Water Level-----
WLPath<-here("Data", "AllSled", "WL")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

WL_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(Date = date, Depth)

## NOw the chr ones
WLPath<-here("Data", "AllSled", "WL","WL_chr")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

WL_Varari2<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(Date = date, Depth)%>%
  mutate(Date = mdy_hm(Date))

# Bring them together
WL_Varari<-WL_Varari %>%
  bind_rows(WL_Varari2)%>%
  arrange(Date)

### Bring all together ####

All_Date<-CT_Varari %>%
  full_join(WL_Varari) %>%
  mutate(DateTime = Date,
         Date = date(Date)) %>%
  group_by(Date)%>%
  mutate(Sal_norm = Salinity_psu -max(Salinity_psu)) %>% # normalize to daily max
  ungroup()

All_Date %>%
  drop_na(Depth)%>%
  droplevels()%>%
  ggplot(aes(x = log(Depth), y = log(Salinity_psu)))+
  geom_point()+
  labs(y = "Difference from max daily salinity (psu)",
       x = "Water Depth (m)")+
  facet_wrap(~Date, scales = "free")


Data_summary <-
  All_Date %>%
  group_by(Date)%>%
  summarise_at(vars(TempInSitu:Depth), .funs = list(mean = mean, max = max, min = min, var = var))

Data_summary %>%
  ggplot(aes(x = log(Depth_mean), y =log(Salinity_psu_var)))+
  geom_point()

Data_summary %>%
  ggplot(aes(x = log(TempInSitu_min), y =log(Salinity_psu_min)))+
  geom_point()+
  geom_smooth(method = "lm", formula = "y~poly(x,2)")
