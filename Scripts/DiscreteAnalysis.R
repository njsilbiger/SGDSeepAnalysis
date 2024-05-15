### discrete Seep data ####


library(tidyverse)
library(here)


##### read in data ############

Data<-read_csv(here("Data","Varari","DiscreteSeep.csv"))


#### plot some relationships

Data %>% select(salinity:ammonia_umolL, pH, Seep_Spring, Year_Month)%>%
  pivot_longer(cols = TA:pH)%>%
  filter(Seep_Spring == "Seep")%>%
  ggplot(aes(x = salinity,y = value))+
  geom_point(aes(color = Year_Month))+
  geom_smooth(method = "lm")+
#  coord_trans(y = "log")+
  theme_bw()+
  facet_wrap(~name, scales = "free")
  
