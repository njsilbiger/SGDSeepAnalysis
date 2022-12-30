### SEM script for Moorea reef data ###
### By Nyssa Silbiger ####
### Created on 12/29/2022 ####

#### Load libraries ############
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(ggrepel)
library(patchwork)
library(brms)
library(tidybayes)
library(bayesplot)
library(stringr)
library(broom)



#### Read in the data ####
# read in the data from the pca script so that all data cleaning is the same
source(here("Scripts","EcoMetabScript.R"))

rm(list= ls()[!(ls() %in% c("Cdata","Data", "Datalog", "remove_varari","remove_cabral","remove_vararilog","remove_cabrallog","turb_all"))])

# extract the model data
ModelData<-Cdata %>%
  anti_join(remove_varari)%>%
  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  select(Location, Tide, Day_Night, TimeBlock, Season, pH, Salinity, Silicate_umolL, NEP, NEP.proxy, NEC, NEC.proxy, NN_umolL, Phosphate_umolL, Ammonia_umolL, Temperature, VisibleHumidic_Like, MarineHumic_Like, Tryptophan_Like, Tyrosine_Like, TA) %>%
  mutate(Silicate_umolL = log(Silicate_umolL), # log transform all the nutrient data
         NN_umolL = log(NN_umolL),
         Phosphate_umolL = log(Phosphate_umolL),
         Ammonia_umolL = log(Ammonia_umolL),
         Humics = VisibleHumidic_Like+MarineHumic_Like,
         Proteinaceous = Tryptophan_Like +Tyrosine_Like) %>% 
  mutate_at(vars("pH":"Proteinaceous"),  function(x) scale(x, center = TRUE)[,1]) 

### Create the models ####

SGDmod<-bf(Silicate_umolL~Salinity)
NNmod<-bf(NN_umolL ~Silicate_umolL) ### N and P are highly correlated 
Pmod<-bf(Phosphate_umolL ~Silicate_umolL)
NH4mod<-bf(Ammonia_umolL ~Silicate_umolL) # N and NH4 are weakly correlated
NEPmod<-bf(NEP.proxy~Day_Night*(Temperature+NN_umolL)) ## N or P or NH4?
pHmod<-bf(pH~ NEP.proxy+Silicate_umolL)
NECmod<-bf(NEC.proxy~pH+Temperature)
Humicsmod<-bf(Humics~NEP.proxy*Day_Night+Silicate_umolL)
Protmod<-bf(Proteinaceous~NEC.proxy+Day_Night)

# Function to run Bayesian SEM and make the posterior predictive checks and plot marginal effects
RunSEM<-function(site){
  
  
  fit_brms <- brm(SGDmod+
                    NNmod+
                    NEPmod+
                    pHmod+
                    NECmod+
                    Humicsmod+
                    Protmod+
                    set_rescor(FALSE),
    data=ModelData[ModelData$Location == site,]
    ,cores=4, chains = 3)
  # calculate LOO (leave one out) diagnostics
  SGD_loo<-loo(fit_brms, reloo = TRUE) # looks good!
  
   p1<-pp_check(fit_brms, resp="SilicateumolL") +
     scale_color_manual(values=c("red", "black"))+
     ggtitle("Silicate")
  p2<-pp_check(fit_brms, resp="NNumolL") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("N+N")
  p3<-pp_check(fit_brms, resp="NEPproxy") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("NEP")

  p4<-pp_check(fit_brms, resp="pH") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("pH")
  
  p5<-pp_check(fit_brms, resp="NECproxy") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("NEC")
  
  p6<-pp_check(fit_brms, resp="Humics") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("Humics")
  
  p7<-pp_check(fit_brms, resp="Proteinaceous") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("Proteinaceous")
  # p6<-pp_check(fit_brms, resp="excessco2umolkgstd") +
  #   scale_color_manual(values=c("red", "black"))+
  #   ggtitle("Excess CO2")
  
  p1+ p2+p3+p4+p5+p6+p7+
    #p6+
    plot_layout(guides = "collect") +
    plot_annotation(title = 'Posterior Predictive Checks', tag_levels = "A")
  ggsave(here("Output",paste(site,"Posteriorchecks.pdf")), width = 5, height = 5)
  
  # plot the results
  # Model 1
  # Silicate ~ Silicate
  
  R<-conditional_effects(fit_brms, "Silicate_umolL", resp = "rnbqm3std", method = "predict", resolution = 1000)
  R1<-R$rnbqm3std.rnbqm3std_water_level_m_std %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center"),
           lower = lower__*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center"),
           upper = upper__*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")) %>%
    mutate(WaterLevel = water_level_m_std*attr (SGDData$water_level_m_std,"scaled:scale")+attr(SGDData$water_level_m_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = WaterLevel, y = estimate), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = WaterLevel,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = "grey")+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = water_level_m, y = rn_bq_m3, color = DayNight)) +
    xlab("Water level (m)")+
    ylab(expression(atop("Radon", paste("(bq m"^3,")"))))+
    ggtitle("Model 1")+
   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  ## Model 2
  # N ~ SGD*DayNight
  R<-conditional_effects(fit_brms, "rn_bq_m3_std:DayNight", resp = "noxumstd", method = "predict", resolution = 1000)
  R2<-R$`noxumstd.noxumstd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center"),
           lower = lower__*attr(SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center"),
           upper = upper__*attr(SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center")) %>%
    mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = nox_u_m, color = DayNight)) +
    ylab(expression(atop("Nitrate + Nitrite", paste("(",mu, "mol L"^-1,")"))))+
    xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
    ggtitle("Model 2")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 3
  # Temperature ~ SGD*DayNight
  
  R<-conditional_effects(fit_brms, "rn_bq_m3_std:DayNight", resp = "temperaturestd", method = "predict", resolution = 1000)
  R3<-R$`temperaturestd.temperaturestd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center"),
           lower = lower__*attr(SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center"),
           upper = upper__*attr(SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center")) %>%
    mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = temperature, color = DayNight)) +
    ylab(expression(atop("Temperature",paste("(", degree, "C)"))))+
    xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
    ggtitle("Model 3")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # # Model 4a
  # # pH ~ DayNight*(AOU + SGD)
  # R<-conditional_effects(fit_brms, "rn_bq_m3_std:DayNight", resp = "phinstd", method = "predict", resolution = 1000)
  # R4<-R$`phinstd.phinstd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
  #   mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
  #   geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = p_h_in, color = DayNight)) +
  #   ylab("pH")+
  #   xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  #   ggtitle("Model 4")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  # 
  # Model 4a
  # pH ~ DayNight*(AOU + SGD)
  # R<-conditional_effects(fit_brms, "rn_bq_m3_std", resp = "phinstd", method = "predict", resolution = 1000)
  # R4<-R$phinstd.phinstd_rn_bq_m3_std %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
  #   mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = radon, y = estimate), lwd = 1)+
  #   geom_ribbon(aes(x = radon,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = p_h_in)) +
  #   ylab("pH")+
  #   xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  #   ggtitle("Model 4")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 4a
  # pCO2 ~ DayNight*(AOU + SGD)
  R<-conditional_effects(fit_brms, "rn_bq_m3_std", resp = "pcoinuatmstd", method = "predict", resolution = 1000)
  R4<-R$pcoinuatmstd.pcoinuatmstd_rn_bq_m3_std %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center"),
           lower = lower__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center"),
           upper = upper__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center")) %>%
    mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = radon, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = radon,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = p_co_in_uatm)) +
    ylab(expression(paste("pCO"[2]," ", mu, "atm")))+
    xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
    ggtitle("Model 4")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 4b
  # pH ~ DayNight*(AOU + SGD)
  # R<-conditional_effects(fit_brms, "aou_umol_l_std:DayNight", resp = "phinstd", method = "predict", resolution = 1000)
  # R4b<-R$`phinstd.phinstd_aou_umol_l_std:DayNight` %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
  #   mutate(aou = aou_umol_l_std*attr (SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = aou, y = estimate, color = DayNight), lwd = 1)+
  #   geom_ribbon(aes(x = aou,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = p_h_in, color = DayNight)) +
  #   ylab("pH")+
  #   xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ggtitle("Model 4")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # pH ~ DayNight*(AOU + SGD)
  # R<-conditional_effects(fit_brms, "aou_umol_l_std", resp = "phinstd", method = "predict", resolution = 1000)
  # R4b<-R$phinstd.phinstd_aou_umol_l_std %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
  #   mutate(aou = aou_umol_l_std*attr (SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = aou, y = estimate), lwd = 1)+
  #   geom_ribbon(aes(x = aou,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = p_h_in)) +
  #   ylab("pH")+
  #   xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ggtitle("Model 4")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  R<-conditional_effects(fit_brms, "aou_umol_l_std", resp = "pcoinuatmstd", method = "predict", resolution = 1000)
  R4b<-R$pcoinuatmstd.pcoinuatmstd_aou_umol_l_std %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center"),
           lower = lower__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center"),
           upper = upper__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center")) %>%
    mutate(aou = aou_umol_l_std*attr (SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = aou, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = aou,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = p_co_in_uatm)) +
    ylab(expression(paste("pCO"[2]," ", mu, "atm")))+
    xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ggtitle("Model 4")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  
  # Model 5
  # AOU ~ DayNight*(N+Temperature)
  R<-conditional_effects(fit_brms, "temperature_std", resp = "aouumollstd", method = "predict", resolution = 1000)
  R5<-R$aouumollstd.aouumollstd_temperature_std %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           lower = lower__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           upper = upper__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")) %>%
    mutate(temperature = temperature_std*attr (SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = temperature, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = temperature,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = temperature, y = aou_umol_l)) +
    xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
    ylab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ggtitle("Model 5")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # R<-conditional_effects(fit_brms, "temperature_std:DayNight", resp = "aouumollstd", method = "predict", resolution = 1000)
  # R5<-R$`aouumollstd.aouumollstd_temperature_std:DayNight` %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")) %>%
  #   mutate(temperature = temperature_std*attr (SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = temperature, y = estimate, color = DayNight), lwd = 1)+
  #   geom_ribbon(aes(x = temperature,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = temperature, y = aou_umol_l, color = DayNight)) +
  #   xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  #   ylab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ggtitle("Model 5")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  # 
  # # Model 5b
  # # AOU ~ DayNight*(N+Temperature)
  # 
  R<-conditional_effects(fit_brms, "nox_u_m_std:DayNight", resp = "aouumollstd", method = "predict", resolution = 1000)
  R5b<-R$`aouumollstd.aouumollstd_nox_u_m_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           lower = lower__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           upper = upper__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")) %>%
    mutate(NN = nox_u_m_std*attr (SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NN, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = NN,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = nox_u_m, y = aou_umol_l, color = DayNight)) +
    xlab(expression(atop("Nitrate + Nitrite", paste("(",mu, "mol L"^-1,")"))))+
    ylab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ggtitle("Model 5")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 6
  # excess CO2 ~ DayNight*(SGD +AOU)
  
  # R<-conditional_effects(fit_brms, "rn_bq_m3_std:DayNight", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
  #   R6<-R$`excessco2umolkgstd.excessco2umolkgstd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
  #     mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #            lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #            upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
  #     mutate(radon = rn_bq_m3_std*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
  #     )%>%
  #     ggplot()+ # back trasform the log transformed data for better visual
  #     geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
  #     geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  #     geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = excess_co2_umol_kg, color = DayNight)) +
  #     xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  #     ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
  #     ggtitle("Model 6")+
  #     # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #     #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #     theme_minimal()+
  #     theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # R<-conditional_effects(fit_brms, "rn_bq_m3_std", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
  # R6<-R$excessco2umolkgstd.excessco2umolkgstd_rn_bq_m3_std %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
  #   mutate(radon = rn_bq_m3_std*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = radon, y = estimate), lwd = 1)+
  #   geom_ribbon(aes(x = radon,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = excess_co2_umol_kg)) +
  #   xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  #   ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
  #   ggtitle("Model 6")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 6b
  # excess CO2 ~ DayNight*(SGD +AOU)
  
  # R<-conditional_effects(fit_brms, "aou_umol_l_std:DayNight", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
  # R6b<-R$`excessco2umolkgstd.excessco2umolkgstd_aou_umol_l_std:DayNight` %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
  #   mutate(aou = aou_umol_l_std*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = aou, y = estimate, color = DayNight), lwd = 1)+
  #   geom_ribbon(aes(x = aou,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = excess_co2_umol_kg, color = DayNight)) +
  #   xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
  #   ggtitle("Model 6")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  # 
  # R<-conditional_effects(fit_brms, "aou_umol_l_std", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
  # R6b<-R$excessco2umolkgstd.excessco2umolkgstd_aou_umol_l_std %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
  #   mutate(aou = aou_umol_l_std*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = aou, y = estimate), lwd = 1)+
  #   geom_ribbon(aes(x = aou,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = excess_co2_umol_kg)) +
  #   xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
  #   ggtitle("Model 6")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  ## bring them all together in patchwork
  R<-(R2|R3)/(R4|R4b)/(R5|R5b)+
    #/(R6|R6b)+
    plot_layout(guides = "collect")+
    plot_annotation(tag_levels = "A")&
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 11))
  
  ggsave(here("Output",paste(site,"marginaleffects.pdf")),R, width = 12, height = 18, useDingbats = FALSE)
  
  var_name <- paste(fit_brms, site, sep="_") # Construct the name
  var_name <- str_replace_all(var_name, " ", "")   # remove white space
  return(fit_brms)
  
}
