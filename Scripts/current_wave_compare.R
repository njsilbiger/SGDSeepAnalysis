## plot significant wave height versus current speed.
library(tidyverse)


current1<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/ADCP/210804/210804_summary.csv")

current2<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/ADCP/210821/210821_summary.csv") %>%
  mutate(DateTime = DateTime+minutes(3))

current3<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/ADCP/220320/220320_summary.csv")%>%
  mutate(DateTime = DateTime+minutes(1))

current<-bind_rows(current1, current2, current3)

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

mod2<-lm(Depth ~ Significant_wave_height, data = currentdata)
anova(mod2)
summary(mod2)

mod3<-lm(Speed ~ Depth, data = currentdata)
anova(mod3)
summary(mod3)


## residence time thought since I am looking at difference from mixing line and water is always coming from the direction of the seep... use distance from seep to point as the length parameter, average site depth as depth (use actual depth too to test), and current speed at time of collection to calculate water residence time

#residence time is ((1/current velocity (m/s))*distance to seep)/60  = hours