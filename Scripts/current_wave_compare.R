## plot significant wave height versus current speed.
library(tidyverse)


current1<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/ADCP/210804/210804_summary.csv")

current2<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/ADCP/210821/210821_summary.csv") %>%
  mutate(DateTime = DateTime+minutes(3))

current3<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/ADCP/220320/220320_summary.csv")%>%
  mutate(DateTime = DateTime+minutes(1))

current<-bind_rows(current1, current2, current3)

## read in Seep Data

AllVarari<-read_csv(here("Data","Varari","AllVarariSeepData.csv"))

## join with All Varari
currentdata<-AllVarari %>%
 # drop_na(Significant_wave_height)%>%
  mutate(DateTime = if_else(minute(date) == 39, date,date-minutes(1))) %>%
  left_join(current) %>%
  drop_na(Speed)
  
currentdata %>%
ggplot(aes(x = Significant_wave_height, y = Speed))+
  geom_point()+
  geom_smooth(method = "lm")

currentdata %>%
  ggplot(aes(x = Significant_wave_height, y = Depth))+
  geom_point()+
  geom_smooth(method = "lm")

currentdata %>%
  ggplot(aes(x = Depth, y = Speed))+
  geom_point(aes(color = Season))+
  geom_smooth(method = "lm")

mod1<-lm(Speed ~ Significant_wave_height, data = currentdata)
anova(mod1)
summary(mod1)


## predict current speed from significant wave height for the missing data
# y = 0.14489x-0.05347

mod2<-lm(Depth ~ Significant_wave_height, data = currentdata)
anova(mod2)
summary(mod2)

mod3<-lm(Speed ~ Depth, data = currentdata)
anova(mod3)
summary(mod3)


## residence time thought since I am looking at difference from mixing line and water is always coming from the direction of the seep... use distance from seep to point as the length parameter, average site depth as depth (use actual depth too to test), and current speed at time of collection to calculate water residence time

#residence time is ((1/current velocity (m/s))*distance to seep)/60  = hours

## extract sampling times 

## Add in sampling times for Varari and Cabral
Varari_sample<-tibble(DateTime = ymd_hms(c("2021-08-05 11:57:00", "2021-08-05 00:00:00",
                                           "2021-08-08 18:30:00", "2021-08-06 06:40:00", "2021-08-08 07:30:00",
                                           "2022-03-21 5:00:00","2022-03-21 8:00:00","2022-03-21 11:00:00",
                                           "2022-03-21 14:00:00","2022-03-21 17:00:00","2022-03-21 20:00:00",
                                           "2022-03-21 21:45:00","2022-03-21 23:00:00","2022-03-22 2:00:00",
                                           "2022-03-28 17:00:00","2022-03-28 18:00:00","2022-03-28 19:00:00",
                                           "2022-03-29 7:00:00","2022-03-29 8:00:00","2022-03-29 9:00:00",
                                           "2022-03-29 11:00:00")),
                      matchdate = ymd_hms(c("2021-08-05 11:57:00", "2021-08-05 00:01:00", # to get it to match with the times on the ADCP
                                            "2021-08-08 18:31:00", "2021-08-06 06:41:00", "2021-08-08 07:31:00",
                                            "2022-03-21 5:01:00","2022-03-21 8:01:00","2022-03-21 11:01:00",
                                            "2022-03-21 14:01:00","2022-03-21 17:01:00","2022-03-21 20:01:00",
                                            "2022-03-21 21:45:00","2022-03-21 23:01:00","2022-03-22 2:01:00",
                                            "2022-03-28 17:01:00","2022-03-28 18:01:00","2022-03-28 19:01:00",
                                            "2022-03-29 7:01:00","2022-03-29 8:01:00","2022-03-29 9:01:00",
                                            "2022-03-29 11:01:00"
                        
                      )))


Cabral_sample<-tibble(DateTime =ymd_hms(c("2021-08-09 07:00:00", "2021-08-09 13:00:00",
                                          "2021-08-09 01:10:00","2021-08-09 19:00:00","2021-08-10 07:00:00",
                                          "2022-03-30 04:00:00", "2022-03-30 07:00:00", "2022-03-30 10:00:00", 
                                          "2022-03-30 13:00:00", "2022-03-30 16:00:00", "2022-03-30 19:00:00",
                                          "2022-03-30 22:00:00", "2022-03-31 01:00:00")))

## join with speed

## calculate current speed from the with the missing values

waves<-read_csv(here("Data","IslandData","WestSideADCPMCR.csv")) %>%
  filter(datetime > ymd("2021-08-07") & datetime < ymd("2021-09-06")) %>%
  mutate(Season = "Dry",
         date = datetime)%>%
  select(!datetime)

# at 2021-08-08 18:31:00 the wave height is 1.7075
# at 2021-08-08 07:31:00 the wave height is 2.285
missing<- tibble(DateTime = ymd_hms("2021-08-08 18:31:00","2021-08-08 07:31:00"),
                 waveheight = c(1.7075,2.285)) %>%
  mutate(Speed =0.14489*waveheight-0.05347 )%>%
  select(!waveheight)



currentspeed<-current %>%
  bind_rows(missing)%>%
  rename(matchdate = DateTime) %>%
  right_join(Varari_sample) %>%
  select(DateTime, CurrentSpeed_m_s = Speed) 


write_csv(currentspeed,here("Data","currentspeed.csv"))


