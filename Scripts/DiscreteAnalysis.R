### discrete Seep data ####


library(tidyverse)
library(here)


##### read in data ############

Data<-read_csv(here("Data","Varari","DiscreteSeep.csv"))


#### plot some relationships

Data %>% select(salinity:ammonia_umolL, pH, TN, TP, Seep_Spring, Year_Month)%>%
   mutate(TON = TN - (NN_umolL+ammonia_umolL),
         TOP = TP-phosphate_umolL) %>% # calculate TON and TOP for the few samples we have
  pivot_longer(cols = c(TA:TP, TON, TOP))%>%
  filter(Seep_Spring == "Seep")%>%
  ggplot(aes(x = salinity,y = value))+
  geom_point(aes(color = Year_Month))+
  geom_smooth(method = "lm")+
#  coord_trans(y = "log")+
  theme_bw()+
  facet_wrap(~name, scales = "free")

# calculate TON and TOP for the few samples we have

Data %>%
  mutate(TON = TN - NN_umolL-ammonia_umolL,
         TOP = TP-phosphate_umolL)
  
