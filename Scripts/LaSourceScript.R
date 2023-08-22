#### Tahiti "La Source" Plot ####


library(tidyverse)
library(here)
library(ggrepel)
library(lubridate)
library(patchwork)


## read in the data ###
LaSource<-read_csv(here("Data","Tahiti","Calibrated_CT_330_La_source_2023_0625.csv"))%>%
  mutate(Site = "La Source")

Nordhoff<-read_csv(here("Data","Tahiti","Calibrated_CT_322_Nordhoff_20230625.csv"))%>%
  mutate(Site = "Nordhoff")

AllData <-
  bind_rows(LaSource, Nordhoff)

p1<-AllData %>%
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



# Now with temperature
p2<-AllData %>%
  mutate(label = if_else(Date == max(Date), as.character(Site), NA_character_)) %>% 
  ggplot(aes(x = Date, y = TempInSitu, color = Site))+
  geom_line(linewidth = 1.2)+
  scale_color_manual(values = c("#244474","#907a71"))+
  labs(y = expression("Temperature " ( degree*C)),
       x = "")+
  geom_text(data = AllData %>% filter(Date == last(Date)), aes(label = Site, x = Date + minutes(80), y = TempInSitu, color = Site), size = 4) +
  guides(color = FALSE)+
  theme_bw()+
  theme(text = element_text(size = 14))
  
#ggsave(here("Output","LaSource_temp.png"), width = 10, height = 6)


p1/p2
ggsave(here("Output","LaSource.png"), width = 10, height = 6)

## plot relationship between Temp and Salinity
p3<-AllData %>%
 ggplot(aes(y  = TempInSitu, x = Salinity_psu, color = Site))+
  geom_point(alpha = 0.2)+
  annotate("text",x = 30, y = 27, label = "La Source", color = "#244474", size = 5)+
  annotate("text",x = 34, y = 28.35, label = "Nordhoff", color = "#907a71",
           size = 5)+
  scale_color_manual(values = c("#244474","#907a71"))+
  geom_smooth(method = "lm", formula = "y~poly(x,2)")+
  labs(x = "Salinity (psu)",
       y = expression("Temperature " ( degree*C)))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

p1/p2|p3
ggsave(here("Output","LaSource_all.png"), width = 14, height = 6)
