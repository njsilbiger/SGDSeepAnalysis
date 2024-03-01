#### All the sled data ####
### By Nyssa Silbiger ###
### 2023-08-20 ####

##load libraries ###

library(here)
library(tidyverse)
library(lubridate)

###

# bring in the tide data
TidePath<-here("Data","AllSled","AllTides")
files<-dir(path = TidePath, pattern = ".txt", full.names = TRUE)

Tides<-files %>%
  set_names()%>%
  map_df(~read_tsv(.,skip = 13), .id = "filename")%>%
  mutate(Date = ymd_hms(paste(Date,Day))) %>%
  select(Date,tideheight = Time)

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
  left_join(Tides)%>%
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
  ggplot(aes(x = Depth_min, y =Salinity_psu_min))+
  geom_point()

Data_summary %>%
  ggplot(aes(y = log(TempInSitu_min), x =log(Salinity_psu_min)))+
  geom_point()+
  geom_smooth(method = "lm", formula = "y~poly(x,2)")


All_Date %>%
  ggplot(aes(x = Salinity_psu, y = TempInSitu))+
  geom_point()

# create a 10 min moving average
All_Date<-All_Date %>%
  mutate(Sal_10 = zoo::rollmean(Salinity_psu, 60, na.pad = TRUE),
         Temp_10 = zoo::rollmean(TempInSitu, 60, na.pad = TRUE),
         Depth_10 = zoo::rollmean(Depth, 60, na.pad = TRUE),
         Tide_10 = zoo::rollmean(tideheight, 60, na.pad = TRUE),
         month = month(Date),
         season = case_when(month %in% c(12,1,2,3,4,5)~"Rainy",
                            month %in% c(6:11)~"Dry"
         )
         ) %>%
  mutate(rise_fall = ifelse(Tide_10 > lag(Tide_10), "Rising","Falling")) # add a column if the tide is rising or falling


All_Date %>%
  filter(Date == ymd("2021-08-04")) %>%
  ggplot()+
  geom_line(aes(x = DateTime, y = TempInSitu))+
  geom_line(aes(x = DateTime, y = Temp_10), color = "red")

All_Date %>% 
 # filter(Date == ymd("2021-08-06")) %>%
  ggplot(aes(x = Temp_10, y = Sal_10, color = rise_fall))+
  geom_point()+
  facet_wrap(season~rise_fall, scales = "free_x")

All_Date%>%
  ggplot(aes(x = season, y = Salinity_psu))+
  geom_boxplot()

All_Date%>%
  ggplot(aes(x = season, y = Depth))+
  geom_boxplot()
