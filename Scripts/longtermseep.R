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
