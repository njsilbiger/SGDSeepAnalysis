## Clean windguru data
## Created by Nyssa Silbiger
## Edited on 9/23/2020
############

### Load libraries #####

library(here)
library(tidyverse)
library(lubridate)

#### Read in data ####
# wind
wind<-read_csv(here("Data","IslandData","wind2022.csv")) %>%
  mutate_if(is.numeric, as.character)%>% # for some reason its showing as a character
  pivot_longer(names_to = "time", cols = `00h`:`23h`, values_to = "wind") %>%
  mutate(wind = as.numeric(wind),
         date = dmy_h(paste(`...1`,time))) %>%
  select(date,wind)
  
# rain ################3
rain<-read_csv(here("Data","IslandData","rain2022.csv")) %>%
  mutate_if(is.numeric, as.character)%>% # for some reason its showing as a character
  pivot_longer(names_to = "time", cols = `00h`:`23h`, values_to = "rain") %>%
  mutate(rain = as.numeric(rain),
         rain = ifelse(is.na(rain),0,rain), # if na make a 0
         date = dmy_h(paste(`...1`,time))) %>%
  select(date,rain)

# waves

waves<-read_csv(here("Data","IslandData","waves2022.csv")) %>%
  mutate_if(is.numeric, as.character)%>% # for some reason its showing as a character
  pivot_longer(names_to = "time", cols = `00h`:`23h`, values_to = "waves") %>%
  mutate(waves = as.numeric(waves),
         date = dmy_h(paste(`...1`,time))) %>%
  select(date,waves)

### join everything ####
weather<-rain %>%
  left_join(wind)%>%
  left_join(waves)

## export it
write_csv(weather, here("Data","IslandData","weather2022.csv"))


weather %>%
  pivot_longer(cols = rain:waves, names_to = 'Params',values_to = "Values") %>%
  ggplot(aes(x = date, y = Values))+
  geom_line()+
  facet_wrap(~Params, scales = "free_y")+
  theme_bw()
