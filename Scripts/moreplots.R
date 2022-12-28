
#More Plots!

Data %>%
  filter( Location == "Varari", Plate_Seep == "Plate", Season == "Dry") %>%
  mutate(SGDsupressed = if_else(DateTime == ymd_hms("2021-08-06 06:40:00"), "yes","no")) %>%
  ggplot(aes(x = pH, y =paste(Tide,SGDsupressed, Day_Night), fill = SGDsupressed))+
  geom_density_ridges()





Cdata %>%
  mutate(NewDay_Night = case_when(
    Tide == "Low" & Day_Night == "Day" ~ "Night",
    Tide == "Low" & Day_Night == "Night" ~ "Day",
    Tide == "High" & Day_Night == "Day" ~ "Day",
    Tide == "High" & Day_Night == "Night" ~ "Night"
      ) 
           
           ) %>%
  filter(Plate_Seep == "Plate",
         DateTime != ymd_hms("2021-08-06 06:40:00")) %>%
  ggplot(aes(x = DIC, y = TA, color = Tide, shape = NewDay_Night)) +
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Season*Location, scales = "free")


Cdata %>%
  mutate(NewDay_Night = case_when(
    Tide == "Low" & Day_Night == "Day" ~ "Night",
    Tide == "Low" & Day_Night == "Night" ~ "Day",
    Tide == "High" & Day_Night == "Day" ~ "Day",
    Tide == "High" & Day_Night == "Night" ~ "Night"
  ) 
  
  ) %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(y = TA.diff/2, x = DIC.diff - (TA.diff/2), shape = Tide, color = NewDay_Night))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location*Season, scales = "free")


Cdata %>%
  mutate(NewDay_Night = case_when(
    Tide == "Low" & Day_Night == "Day" ~ "Night",
    Tide == "Low" & Day_Night == "Night" ~ "Day",
    Tide == "High" & Day_Night == "Day" ~ "Day",
    Tide == "High" & Day_Night == "Night" ~ "Night"
  ) 
  
  ) %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(y = TA.diff/2, x = paste(NewDay_Night, Tide)))+
  geom_boxplot()+
  facet_wrap(~Location*Season, scales = "free")



Cdata %>%
  mutate(NewDay_Night = case_when(
    Tide == "Low" & Day_Night == "Day" ~ "Night",
    Tide == "Low" & Day_Night == "Night" ~ "Day",
    Tide == "High" & Day_Night == "Day" ~ "Day",
    Tide == "High" & Day_Night == "Night" ~ "Night"
  ) 
  
  ) %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(y = DIC.diff - (TA.diff/2), NewDay_Night, fill =  Tide))+
  geom_boxplot()+
  facet_wrap(~Location*Season, scales = "free")


Cdata %>%
  mutate(NewDay_Night = case_when(
    Tide == "Low" & Day_Night == "Day" ~ "Night",
    Tide == "Low" & Day_Night == "Night" ~ "Day",
    Tide == "High" & Day_Night == "Day" ~ "Day",
    Tide == "High" & Day_Night == "Night" ~ "Night"
  ) 
  
  ) %>%
  filter(NewDay_Night == "Night", Tide == "Low", Location == "Varari") %>%
  mutate(SGDsupressed = if_else(DateTime == ymd_hms("2021-08-06 06:40:00"), "yes","no")) %>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(y = (TA.diff/2), SGDsupressed))+
  geom_boxplot()

Cdata <- Cdata%>%
  mutate(NewDay_Night = case_when(
    Tide == "Low" & Day_Night == "Day" ~ "Night",
    Tide == "Low" & Day_Night == "Night" ~ "Day",
    Tide == "High" & Day_Night == "Day" ~ "Day",
    Tide == "High" & Day_Night == "Night" ~ "Night"
  ) 
  
  ) 
  

Cdata %>%
filter(NewDay_Night == "Night", Tide == "Low", Location == "Varari") %>%
  mutate(SGDsupressed = if_else(DateTime == ymd_hms("2021-08-06 06:40:00"), "yes","no")) %>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(y = DIC.diff - (TA.diff/2), SGDsupressed))+
  geom_boxplot()



Cdata %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = log(NN_umolL), y = DIC.diff - (TA.diff/2), color = NewDay_Night))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")


Cdata %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = NN_umolL, y = (TA.diff/2), color = NewDay_Night, shape =  Tide))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~Location*Season, scales = "free")


Cdata %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = log(Silicate_umolL), y = DIC.diff - (TA.diff/2), color = NewDay_Night))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(x = log(Silicate_umolL), y = (TA.diff/2), color = NewDay_Night))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")

Cdata %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes(y = VisibleHumidic_Like, x = DIC.diff - (TA.diff/2), color = NN_umolL))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Net community production")+
  facet_wrap(~NewDay_Night*Location, scales = "free")

Cdata %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate", HIX<10) %>%
  ggplot(aes(y = VisibleHumidic_Like, x = DIC.diff - (TA.diff/2), color = Tide))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "N_umol")+
  facet_wrap(~NewDay_Night*Location, scales = "free")


Cdata %>%
  anti_join(removelow)%>%
  filter(Plate_Seep == "Plate") %>%
  ggplot(aes (x =log(Silicate_umolL), y = VisibleHumidic_Like, color = NewDay_Night))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Location*Season, scales = "free")
