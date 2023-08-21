#### Tahiti "La Source" Plot ####


library(tidyverse)
library(here)
library(ggrepel)
library(lubridate)


## read in the data ###
LaSource<-read_csv(here("Data","Tahiti","Calibrated_CT_330_La_source_2023_0625.csv"))%>%
  mutate(Site = "La Source")

Nordhoff<-read_csv(here("Data","Tahiti","Calibrated_CT_322_Nordhoff_20230625.csv"))%>%
  mutate(Site = "Nordhoff")

AllData <-
  bind_rows(LaSource, Nordhoff)

AllData %>%
  mutate(label = if_else(Date == max(Date), as.character(Site), NA_character_)) %>% 
  ggplot(aes(x = Date, y = Salinity_psu, color = Site))+
  geom_line(linewidth = 1.2)+
  scale_color_manual(values = c("#244474","#907a71"))+
  labs(y = "Salinity (psu)",
       x = "")+
  geom_text(data = AllData %>% filter(Date == last(Date)), aes(label = Site, x = Date + minutes(80), y = Salinity_psu, color = Site), size = 4) +
  guides(color = FALSE)+
theme_bw()+
  theme(text = element_text(size = 14))
ggsave(here("Output","LaSource.png"), width = 10, height = 6)
