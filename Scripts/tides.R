### Make a figure of the sampling times by tide
### By Nyssa Silbiger ###
### updated 9/14/2023 ######

## Tide predictions------------ https://tidesandcurrents.noaa.gov/noaatidepredictions.html?id=1732417
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


## Add in sampling times for Varari and Cabral
Varari_sample<-tibble(date = ymd_hms(c("2021-08-05 11:57:00", "2021-08-05 00:00:00",
                                       "2021-08-08 18:30:00", 
                                       "2021-08-08 07:30:00","2021-08-06 19:00:00",
                                       "2021-08-04 18:45:00", "2021-08-04 23:51:00",
                                       "2021-08-05 02:51:00", "2021-08-05 06:45:00", 
                                       "2021-08-05 08:51:00", "2021-08-06 08:30:00",
                                       "2022-03-21 5:00:00", "2022-03-21 8:00:00",
                                       "2022-03-21 11:00:00", "2022-03-21 14:00:00",
                                       "2022-03-21 17:00:00","2022-03-21 20:00:00",
                                        "2022-03-21 21:45:00","2022-03-21 23:00:00",
                                       "2022-03-22 2:00:00", "2022-03-28 17:00:00",
                                       "2022-03-28 18:00:00","2022-03-28 19:00:00",
                                        "2022-03-29 7:00:00","2022-03-29 8:00:00",
                                       "2022-03-29 9:00:00", "2022-03-29 11:00:00")),
                      Plate_Seep = c("All Sites", "All Sites",
                                     "All Sites", 
                                     "All Sites", "Seep Only","Seep Only", "Seep Only", "Seep Only",
                                     "Seep Only", "Seep Only", "Seep Only","Seep Only","All Sites","Seep Only",
                                     "All Sites","Seep Only","All Sites", "Seep Only","Seep Only","All Sites",
                                     "Seep Only","Seep Only","Seep Only", "Seep Only","Seep Only","Seep Only","Seep Only"))

Cabral_sample<-tibble(date =ymd_hms(c("2021-08-09 07:00:00", 
                                      "2021-08-09 13:00:00",
                                      "2021-08-09 01:10:00",
                                      "2021-08-09 19:00:00",
                                      "2021-08-10 07:00:00",
                                          "2021-08-09 04:00:00", "2021-08-09 10:00:00",
                                          "2021-08-09 13:00:00", "2021-08-09 16:10:00","2021-08-09 19:00:00",
                                          "2021-08-09 22:00:00", "2021-08-10 01:20:00", "2021-08-10 04:20:00",
                                          "2021-08-10 07:00:00","2022-03-30 04:00:00", "2022-03-30 07:00:00", 
                                          "2022-03-30 10:00:00", "2022-03-30 13:00:00", "2022-03-30 16:00:00", 
                                          "2022-03-30 19:00:00", "2022-03-30 22:00:00", "2022-03-31 01:00:00")),
                      Plate_Seep = c("All Sites", "All Sites",
                                     "All Sites","All Sites","Seep Only",
                                     "Seep Only", "Seep Only", "Seep Only",
                                     "Seep Only","Seep Only","Seep Only","Seep Only", "Seep Only","Seep Only",
                                     "Seep Only", "All Sites", "Seep Only", 
                                     "All Sites", "Seep Only", "All Sites",
                                     "Seep Only", "All Sites"))


Varari_sample<-Varari_sample %>%
  mutate(date = round_date(date, unit = "hour"))%>% # tides are at hourly 
  left_join(tides) %>%
  mutate(Location = "Varari")

Cabral_sample<-Cabral_sample %>%
  mutate(date = round_date(date, unit = "hour"))%>%
  left_join(tides) %>%
  mutate(Location = "Cabral") 

SamplingTimes<-bind_rows(Varari_sample, Cabral_sample)%>%
  filter(Plate_Seep == "All Sites") %>%
  mutate(Season = ifelse(Season == "Dry", "Dry Season","Wet Season")) 

tides_all<-tides %>%
  filter(between(date, ymd_hms("2021-08-04 06:00:00"),ymd_hms("2021-08-10 06:00:00"))|between(date, ymd_hms("2022-03-21 6:00:00"),ymd_hms("2022-03-31 06:00:00"))
           ) %>% ### filter to just the sampline times
  mutate(Time = format(date, format = "%H:%M:%S"),
    sunrise_sunset = case_when(Time == "06:00:00"~ "Sunrise",
                                Time == "18:00:00"~"Sunset")) %>%
  mutate(Season = ifelse(Season == "Dry", "Dry Season","Wet Season")) 

tides_all %>%
  ggplot(aes(x = date, y = tideheight))+
  geom_line(color = "black")+
  geom_point(data = SamplingTimes, aes(x = date, y = tideheight, color = Location), size = 3)+
  geom_vline(data = tides_all %>% filter(sunrise_sunset  == "Sunrise"),aes(xintercept = date), color = "#FFBF00", alpha = 0.5)+
  geom_vline(data = tides_all %>% filter(sunrise_sunset  == "Sunset"),aes(xintercept = date), color = "grey", alpha = 0.5)+
    scale_color_manual(values = c("#16697A","#D64550"))+
  labs(x = "",
       y = "Tide height (m)")+
  facet_wrap(~Season, scales = "free_x", nrow = 2)+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "none",
        line = element_blank())

ggsave(here("Output","tideplot.pdf"), width = 6, height = 8)
