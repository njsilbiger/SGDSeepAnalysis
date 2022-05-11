### Look for patterns among data collected at the seep versus WW3 data
## Created by Nyssa Silbiger
## Edited on 5/11/2022

## load libraries ####-----------
library(here)
library(tidyverse)
library(lubridate)
library(forecast)
library(patchwork)
library(viridis)
library(janitor)
library(ggh4x)

## Read in the different datasets

## Tide predictions------------
# August
tideAug<-read_tsv(here("Data","IslandData","TidePredictions.txt"), skip = 13) %>%
  mutate(date = ymd_hms(paste(Date,Time)),
         Season = "Dry") %>%
  select(date,tideheight = Pred, Season)

# March-April
tideMarch<-read_tsv(here("Data","IslandData","TidePredictions2022.txt"), skip = 13) %>%
  mutate(date = ymd_hms(paste(Date,Time)),
         Season = "Wet") %>%
  select(date,tideheight = Pred, Season)

tides<-bind_rows(tideAug,tideMarch)

## Weather data (wind, rain, waves from windguru)#####
weather<-read_csv(here("Data","IslandData","weather.csv")) %>%
  bind_rows(read_csv(here("Data","IslandData","weather2022.csv")) )%>%
  left_join(tides)%>%
  drop_na(Season)

# CT----------------------------
CondPath<-here("Data", "Varari", "CT")
files <- dir(path = CondPath,pattern = ".csv", full.names = TRUE)

CT_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,TempInSitu, Salinity_psu)%>%
  mutate(Site = "Varari",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))

# Cabral
CondPath<-here("Data", "Cabral", "CT")
files <- dir(path = CondPath,pattern = ".csv", full.names = TRUE)

CT_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,TempInSitu, Salinity_psu)%>%
  mutate(Site = "Cabral",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))


#Water Level-----
WLPath<-here("Data", "Varari", "WL")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

WL_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date, Depth)%>%
  mutate(Site = "Varari",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))

WLPath<-here("Data", "Cabral", "WL")
files <- dir(path = WLPath,pattern = ".csv", full.names = TRUE)

#Cabral
WL_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date, Depth)%>%
  mutate(Site = "Cabral",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))


## pH ----------------
pHPath<-here("Data", "Varari", "pH")
files <- dir(path = pHPath,pattern = ".csv", full.names = TRUE)

pH_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,pH = pH_total)%>%
  mutate(Site = "Varari",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))

#Cabral
pHPath<-here("Data", "Cabral", "pH")
files <- dir(path = pHPath,pattern = ".csv", full.names = TRUE)

pH_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,pH = pH_total)%>%
  mutate(Site = "Cabral",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))


## Lux ----------------
LUXPath<-here("Data", "Varari", "LUX")
files <- dir(path = LUXPath,pattern = ".csv", full.names = TRUE)

LUX_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,Lux)%>%
  mutate(Site = "Varari",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))

#Cabral
LUXPath<-here("Data", "Cabral", "LUX")
files <- dir(path = LUXPath,pattern = ".csv", full.names = TRUE)

LUX_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,Lux)%>%
  mutate(Site = "Cabral",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))

## PAR -------------------
PARPath<-here("Data", "Varari", "PAR")
files <- dir(path = PARPath,pattern = ".csv", full.names = TRUE)

PAR_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,PAR)%>%
  mutate(Site = "Varari", 
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet")) %>%
  mutate(date = ceiling_date(date, unit = "minutes")) # the par data isn't always on the minute mark

#Cabral
PARPath<-here("Data", "Cabral", "PAR")
files <- dir(path = PARPath,pattern = ".csv", full.names = TRUE)

PAR_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,PAR)%>%
  mutate(Site = "Cabral",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))


## DO ---------------
DOPath<-here("Data", "Varari", "DO")
files <- dir(path = DOPath,pattern = ".csv", full.names = TRUE)

DO_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,DO_mg_L)%>%
  mutate(Site = "Varari",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))

#Cabral
DOPath<-here("Data", "Cabral", "DO")
files <- dir(path = DOPath,pattern = ".csv", full.names = TRUE)

DO_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  select(date,DO_mg_L)%>%
  mutate(Site = "Cabral",
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))

# Radon ---------------------
RnPath<-here("Data", "Varari", "Radon")
files <- dir(path = RnPath,pattern = ".csv", full.names = TRUE)

Rn_Varari<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  clean_names()%>% # the columns names are messy
  select(date = full_date,Rn_dpm_L = radon_in_water_dpm_l)%>%
  mutate(Site = "Varari",
         date = mdy_hm(date),
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))

# Cabral
RnPath<-here("Data", "Cabral", "Radon")
files <- dir(path = RnPath,pattern = ".csv", full.names = TRUE)

Rn_Cabral<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename")  %>% # map everything to a dataframe and put the id in a column called filename
  clean_names()%>% # the columns names are messy
  select(date = full_date,Rn_dpm_L = radon_in_water_dpm_l)%>%
  mutate(Site = "Cabral",
         date = mdy_hm(date),
         Season = ifelse(date < ymd("2022-01-02"),"Dry","Wet"))


### Join everything together

AllVarari<-CT_Varari %>%
  left_join(WL_Varari)%>%
  left_join(pH_Varari) %>%
  left_join(LUX_Varari) %>%
  left_join(PAR_Varari)%>%
  left_join(DO_Varari)%>%
  left_join(Rn_Varari)%>%
  relocate(Site, .before = date) %>% # move the site column
  mutate(PAR_calc = case_when(is.na(PAR) & Season == "Dry" ~ 16778.33+(0-16778.33)*exp(1)^(-exp( -13.29572)*Lux), # calculate PAR when it is missing
                              !is.na(PAR)~PAR,
                              is.na(PAR) & Season =="Wet" ~ 1555.751648 +(0-1555.751648 )*exp(1)^(-exp( -10.777624)*Lux)))  %>%
  select(-Lux,-PAR) %>%  # remove Lux and original PAR
  mutate(PAR_calc = ifelse(PAR_calc<0,0,PAR_calc)) %>% # if PAR is negative make it 0
  write_csv(here("Data","Varari","AllVarariSeepData.csv"))

Allcabral<-CT_Cabral %>%
  left_join(WL_Cabral)%>%
  left_join(pH_Cabral) %>%
  left_join(LUX_Cabral) %>%
  left_join(PAR_Cabral)%>%
  left_join(DO_Cabral)%>%
  left_join(Rn_Cabral)%>%
  relocate(Site, .before = date) %>% # move the site column
  mutate(PAR_calc = case_when(is.na(PAR) & Season == "Dry" ~ 16778.33+(0-16778.33)*exp(1)^(-exp( -13.29572)*Lux),
                              !is.na(PAR)~PAR,
                              is.na(PAR) & Season =="Wet" ~ 1555.751648 +(0-1555.751648 )*exp(1)^(-exp( -10.777624)*Lux))) %>%
  select(-Lux,-PAR) %>% # remove Lux and original PAR
  mutate(PAR_calc = ifelse(PAR_calc<0,0,PAR_calc)) %>% # if PAR is negative make it 0
  write_csv(here("Data","Cabral","AllCabralSeepData.csv"))


#### fill in the missing times with NA for easier plotting
# Need to split up Dry and Wet and bring back together or it will be huge

# Dry first
AllVarari_Dry <- AllVarari %>%
  filter(Season =="Dry") 
AllVarari_Dry<-data.frame(date = seq(AllVarari_Dry$date[1], AllVarari_Dry$date[nrow(AllVarari_Dry)], by = "1 min"))  %>%
  full_join(AllVarari_Dry, by = "date") %>%# fill in NA by 1 min
  mutate(Season = "Dry",
         Site = "Varari") # these got filled with NAs
# Wet
AllVarari_Wet <-AllVarari %>%
  filter(Season =="Wet") # we did every 2 mins in the Wet Season
AllVarari_Wet<-data.frame(date = seq(AllVarari_Wet$date[1], AllVarari_Wet$date[nrow(AllVarari_Wet)], by = "2 min"))  %>%
  full_join(AllVarari_Wet, by = "date") %>%# fill in NA by 2 min
  mutate(Season = "Wet",
         Site = "Varari") # these got filled with NAs

# bring together
AllVarari<-bind_rows(AllVarari_Dry, AllVarari_Wet)%>%
  relocate(Season, .after = Site)

# AllVarari<-data.frame(date = seq(AllVarari$date[1], AllVarari$date[nrow(AllVarari)], by = "1 min"))  %>%
#   full_join(AllVarari, by = "date")
 
# Cabral
# Dry
AllCabral_Dry <- Allcabral %>%
  filter(Season =="Dry") 
AllCabral_Dry<-data.frame(date = seq(AllCabral_Dry$date[1], AllCabral_Dry$date[nrow(AllCabral_Dry)], by = "1 min"))  %>%
  full_join(AllCabral_Dry, by = "date") %>% # fill in NA by 1 min
  mutate(Season = "Dry",
         Site = "Cabral") # these got filled with NAs

# Wet
AllCabral_Wet <- Allcabral %>%
  filter(Season =="Wet") 
AllCabral_Wet<-data.frame(date = seq(AllCabral_Wet$date[1], AllCabral_Wet$date[nrow(AllCabral_Wet)], by = "2 min"))  %>% # 2 mins in the Wet Season
  full_join(AllCabral_Wet, by = "date") %>%# fill in NA by 1 min
  mutate(Season = "Wet",
         Site = "Cabral") # these got filled with NAs

# AllCabral<-data.frame(date = seq(AllCabral$date[1], AllCabral$date[nrow(AllCabral)], by = "1 min"))  %>%
#   full_join(AllCabral, by = "date")

AllCabral<-bind_rows(AllCabral_Dry, AllCabral_Wet) %>%
  relocate(Season, .after = Site)

## Add in sampling times for Varari and Cabral
Varari_sample<-tibble(datetime = ymd_hms(c("2021-08-05 11:57:00", "2021-08-05 00:00:00",
                         "2021-08-08 18:30:00", "2021-08-06 06:40:00", "2021-08-08 07:30:00",
                         "2022-03-21 5:00:00","2022-03-21 8:00:00","2022-03-21 11:00:00",
                         "2022-03-21 14:00:00","2022-03-21 17:00:00","2022-03-21 20:00:00",
                         "2022-03-21 21:45:00","2022-03-21 23:00:00","2022-03-22 2:00:00",
                         "2022-03-28 17:00:00","2022-03-28 18:00:00","2022-03-28 19:00:00",
                         "2022-03-29 7:00:00","2022-03-29 8:00:00","2022-03-29 9:00:00",
                         "2022-03-29 11:00:00")))


Cabral_sample<-tibble(datetime =ymd_hms(c("2021-08-09 07:00:00", "2021-08-09 13:00:00",
                         "2021-08-09 01:10:00","2021-08-09 19:00:00","2021-08-10 07:00:00",
                         "2022-03-30 04:00:00", "2022-03-30 07:00:00", "2022-03-30 10:00:00", 
                         "2022-03-30 13:00:00", "2022-03-30 16:00:00", "2022-03-30 19:00:00",
                         "2022-03-30 22:00:00", "2022-03-31 01:00:00")))

## Simple plot

AllVarari %>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
#  geom_vline(data = Varari_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params*Season, scales = "free", ncol = 2)+
  facetted_pos_scales( # make y limits the same by panel
    y = rep(list(
      scale_y_continuous(limits = c(0, 2.5)),
      scale_y_continuous(limits = c(0, 18)),
      scale_y_continuous(limits = c(0, 2000)),
      scale_y_continuous(limits = c(7.0, 8.2)),
      scale_y_continuous(limits = c(0, 6)),
      scale_y_continuous(limits = c(20, 38)),
      scale_y_continuous(limits = c(24, 32))
    ), each = 2)
  )+
  theme_bw() +
  labs(title = "Varari Sled")
ggsave(here('Output',"Varari_timeseries.pdf"), width = 8, height = 10)

# Just during the sampling times
AllVarari %>%
  filter(date >= ymd("2021-08-04") & date <= ymd("2021-08-09") |
           date >= ymd("2022-03-21") & date <= ymd("2022-03-30"))%>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
 # geom_point()+ # since radon is only every 6 mins the line doesnt work... need to average out above
  geom_line()+
  geom_vline(data = Varari_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params*Season, scales = "free", ncol = 2)+
  theme_bw() +
  labs(title = "Varari Sled")

# Just during sampling times
AllCabral %>%
  filter(date >= ymd("2021-08-09"), date <= ymd("2021-08-10")|
           date >= ymd("2022-03-29"), date <= ymd("2022-04-01")  )%>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
#  geom_point()+
  geom_vline(data = Cabral_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params*Season, scales = "free", ncol = 2)+
  facetted_pos_scales( # make y limits the same by panel
    y = rep(list(
      scale_y_continuous(limits = c(0, 1)),
      scale_y_continuous(limits = c(0, 18)),
      scale_y_continuous(limits = c(0, 2000)),
      scale_y_continuous(limits = c(7.6, 8.2)),
      scale_y_continuous(limits = c(0, 6)),
      scale_y_continuous(limits = c(34, 40)),
      scale_y_continuous(limits = c(26, 32))
    ), each = 2)
  )+
  theme_bw() +
  labs(title = "Cabral Sled")
ggsave(here('Output',"Cabral_timeseries.pdf"), width = 8, height = 10)

#Cabral all
AllCabral %>%
  pivot_longer(cols = TempInSitu:PAR_calc, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  geom_vline(data = Cabral_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params*Season, scales = "free", ncol = 2)+
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
  left_join(weather) %>%# join in the weather data
  relocate(Season, .after = date)

# plot with weather
AllVarari_onehour %>%
  pivot_longer(cols = Depth:tideheight, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  facet_wrap(~Params*Season, scales = "free", ncol = 2)+
  facetted_pos_scales( # make y limits the same by panel
    y = rep(list(
      scale_y_continuous(limits = c(0, 1.5)),
      scale_y_continuous(limits = c(0, 16)),
      scale_y_continuous(limits = c(0, 1500)),
      scale_y_continuous(limits = c(7.0, 8.2)),
      scale_y_continuous(limits = c(0, 2)),
      scale_y_continuous(limits = c(0, 6)),
      scale_y_continuous(limits = c(20, 38)),
      scale_y_continuous(limits = c(24, 32)),
      scale_y_continuous(limits = c(0, 0.3)),
      scale_y_continuous(limits = c(0, 3.0)),
      scale_y_continuous(limits = c(0, 30))
    ), each = 2)
  )+
  theme_bw()
ggsave(here("Output","Varari_ts_hourly.pdf"), width = 6, height = 15 )

# Sampling times
AllVarari_onehour %>%
  filter(date >= ymd("2021-08-04") & date <= ymd("2021-08-09") |
           date >= ymd("2022-03-21") & date <= ymd("2022-03-30"))%>%
  pivot_longer(cols = Depth:tideheight, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  #geom_vline(data = Varari_sample, aes(xintercept = datetime), color = "red")+
  facet_wrap(~Params*Season, scales = "free", ncol = 4)+
  theme_bw()


# Cabral
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
  facet_wrap(~Params*Season, scales = "free", ncol = 2)+
  facetted_pos_scales( # make y limits the same by panel
    y = rep(list(
      scale_y_continuous(limits = c(0, 0.6)),
      scale_y_continuous(limits = c(0, 16)),
      scale_y_continuous(limits = c(0, 1000)),
      scale_y_continuous(limits = c(7.5, 8.2)),
      scale_y_continuous(limits = c(0, 2)),
      scale_y_continuous(limits = c(0, 6)),
      scale_y_continuous(limits = c(34, 40)),
      scale_y_continuous(limits = c(24, 32)),
      scale_y_continuous(limits = c(0, 0.3)),
      scale_y_continuous(limits = c(0, 3.0)),
      scale_y_continuous(limits = c(0, 30))
    ), each = 2)
  )+
  theme_bw()
ggsave(here("Output","Cabral_ts_hourly.pdf"), width = 6, height = 15)

# Sampling times
AllCabral_onehour %>%
  filter(date >= ymd("2021-08-09"), date <= ymd("2021-08-10")|
           date >= ymd("2022-03-29"), date <= ymd("2022-04-01") )%>%
  pivot_longer(cols = Depth:tideheight, names_to = "Params", values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  facet_wrap(~Params*Season, scales = "free", ncol = 4)+
  theme_bw()

## plot waves vs depth and tide vs depth
WD_TP<-AllVarari_onehour %>%
  bind_rows(AllCabral_onehour)%>%
  drop_na(Site)%>%
  ggplot(aes(x = tideheight, y = Depth, color = waves))+
  geom_point()+
  geom_smooth(method = "lm")+
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
  facet_wrap(~Site, ncol = 1, scales = "free_y")+
  facetted_pos_scales( # make y limits the same by panel
    y = list(
      scale_y_continuous(limits = c(0, 0.6)),
      scale_y_continuous(limits = c(0, 1.5))
    )
  )

WD_Wave<-AllVarari_onehour %>%
  bind_rows(AllCabral_onehour)%>%
  drop_na(Site)%>%
  ggplot(aes(x = waves, y = Depth, color = tideheight))+
  geom_point()+
  geom_smooth(method = "lm")+
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
  facet_wrap(~Site, ncol = 1, scales = "free_y")+
  facetted_pos_scales( # make y limits the same by panel
    y = list(
      scale_y_continuous(limits = c(0, 0.6)),
      scale_y_continuous(limits = c(0, 1.5))
    )
  )
  

WD_TP + WD_Wave
ggsave(filename = here("Output","WaterDepth_tide_wave.pdf"), width = 8)

### Depth versus different parameters
AllVarari_onehour %>%
  pivot_longer(cols = c("pH","Salinity_psu","TempInSitu", "DO_mg_L", "Rn_dpm_L"), names_to = "Params", values_to = "Values")%>%
  drop_na()%>%
  ggplot(aes(x = Depth, y  = Values, color = waves))+
  geom_point()+
  facet_wrap(~Params*Season, scales="free", ncol = 2)+
  theme_bw()+
  scale_color_viridis_c(option = "plasma")


ggplot(AllVarari, aes(x = pH, y = DO_mg_L, col = TempInSitu))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_c(option = "plasma") +
  facet_wrap(~Season)


AllCabral_onehour %>%
  pivot_longer(cols = c("pH","Salinity_psu","TempInSitu", "DO_mg_L", "Rn_dpm_L"), names_to = "Params", values_to = "Values")%>%
  drop_na()%>%
  ggplot(aes(x = Depth, y  = Values, color = waves))+
  geom_point()+
  facet_wrap(~Params*Season, scales="free", ncol = 2)+
  theme_bw()+
  scale_color_viridis_c(option = "plasma")

# looking at relationships with radon
AllCabral %>%
  drop_na()%>%
  filter(Rn_dpm_L>10)%>% # really low values look like outliers
ggplot(aes(x = log(Rn_dpm_L), y = log(TempInSitu)))+
  geom_point()

### Plot pH data hand collected from seep on pH HOBO data ####
DiscreteData<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv") 

# This is temporary until we get ther rest of the data
Discrete_wet<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/CarbonateChemistry/pHProbe_Data_calculated_NOTPOcorrect.csv") %>%
  mutate(datetime = mdy_hms(paste(as.character(Date), as.character(SamplingTime))))


AllVarari %>%
  filter(date >= ymd("2021-08-05"), date <= ymd("2021-08-09"))%>%
  ggplot()+
  geom_line(aes(x = date, y = pH))+
  geom_point(data = DiscreteData%>% filter(Plate_Seep =="Seep", Location=="Varari"), aes(x = DateTime, y = pH), color = "red", size =2)+
  geom_line(data = DiscreteData%>% filter(Plate_Seep =="Seep", Location=="Varari"), aes(x = DateTime, y = pH), color = "red")+
  theme_bw() +
  labs(title = "Varari Sled")

AllVarari %>%
  filter(date >= ymd("2022-03-21"), date <= ymd("2022-04-01"))%>%
  ggplot()+
  geom_line(aes(x = date, y = pH))+
  geom_point(data = Discrete_wet%>% filter(CowTagID =="VSEEP"), aes(x = datetime, y = pH), color = "red", size =2)+
  geom_line(data = Discrete_wet%>% filter(CowTagID =="VSEEP"), aes(x = datetime, y = pH), color = "red")+
  theme_bw() +
  labs(title = "Varari Sled")

AllCabral %>%
  filter(date >= ymd("2021-08-08"), date <= ymd("2021-08-10"))%>%
   ggplot(aes(x = date, y = pH))+
  geom_line()+
  geom_point(data = DiscreteData%>% filter(Plate_Seep =="Seep", Location=="Cabral"), aes(x = DateTime, y = pH), color = "red", size =2)+
  geom_line(data = DiscreteData%>% filter(Plate_Seep =="Seep", Location=="Cabral"), aes(x = DateTime, y = pH), color = "red")+
  theme_bw() +
  labs(title = "Cabral Sled")

AllCabral %>%
  filter(date >= ymd("2022-03-29"), date <= ymd("2022-04-01"))%>%
  ggplot()+
  geom_line(aes(x = date, y = pH))+
  geom_point(data = Discrete_wet%>% filter(CowTagID =="CSEEP"), aes(x = datetime, y = pH), color = "red", size =2)+
  geom_line(data = Discrete_wet%>% filter(CowTagID =="CSEEP"), aes(x = datetime, y = pH), color = "red")+
  theme_bw() +
  labs(title = "Cabral Sled")
