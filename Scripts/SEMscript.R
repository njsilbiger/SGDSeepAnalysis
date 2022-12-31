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

Cdata <- Cdata%>%
  anti_join(remove_varari)%>% # remove the bad data
  anti_join(remove_cabral) %>%
  filter(Silicate_umolL <15, NEC.proxy<100)

# extract the model data
ModelData<-Cdata %>%
#  anti_join(remove_varari)%>%
#  anti_join(remove_cabral)%>%
  filter(Plate_Seep == "Plate") %>%
  select(Location, Tide, Day.Night = Day_Night, TimeBlock, Season, pH, Salinity, Silicate_umolL, NEP, NEP.proxy, NEC, NEC.proxy, NN_umolL, Phosphate_umolL, Ammonia_umolL, Temperature, VisibleHumidic_Like, MarineHumic_Like, Tryptophan_Like, Tyrosine_Like, TA) %>%
  mutate(SilicateumolL = log(Silicate_umolL), # log transform all the nutrient data
         NNumolL = log(NN_umolL),
         PhosphateumolL = log(Phosphate_umolL),
         AmmoniaumolL = log(Ammonia_umolL),
         Humics = VisibleHumidic_Like+MarineHumic_Like,
         Proteinaceous = Tryptophan_Like +Tyrosine_Like) %>% 
  mutate_at(vars("pH":"Proteinaceous"),  function(x) scale(x, center = TRUE)) 

### Create the models ####

SGDmod<-bf(Silicate_umolL~Salinity, family = "student")
NNmod<-bf(NNumolL ~SilicateumolL, family = "student") ### N and P are highly correlated 
Pmod<-bf(PhosphateumolL ~SilicateumolL, family = "student")
NH4mod<-bf(AmmoniaumolL ~SilicateumolL, family = "student") # N and NH4 are weakly correlated
NEPmod<-bf(NEP.proxy~(Day.Night*NNumolL) +Temperature, family = "student") ## N or P or NH4?
pHmod<-bf(pH~ SilicateumolL+NEP.proxy + (1|Day.Night), family = "student")
NECmod<-bf(NEC.proxy~pH+Temperature, family = "student")
Humicsmod<-bf(Humics~(Day.Night*NEP.proxy)+SilicateumolL, family = "student")
Protmod<-bf(Proteinaceous~NEC.proxy, family = "student")

# Function to run Bayesian SEM and make the posterior predictive checks and plot marginal effects
RunSEM<-function(site, season){
  
  
  fit_brms <- brm(#SGDmod+
                    NNmod+
                    NEPmod+
                    pHmod+
                    NECmod+
                    Humicsmod+
                    Protmod+
                    set_rescor(FALSE),
    data=ModelData[ModelData$Location == site & ModelData$Season == season,]
    ,cores=4, chains = 3)
  # calculate LOO (leave one out) diagnostics
  SGD_loo<-loo(fit_brms, reloo = TRUE) # looks good!
  
   # p1<-pp_check(fit_brms, resp="SilicateumolL") +
   #   scale_color_manual(values=c("red", "black"))+
   #   ggtitle("Silicate")
  
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

  #p1+ 
    p2+p3+p4+p5+p6+p7+
    #p6+
    plot_layout(guides = "collect") +
    plot_annotation(title = 'Posterior Predictive Checks', tag_levels = "A")
  ggsave(here("Output",paste(site,season, "Posteriorchecks.pdf")), width = 5, height = 5)
  
  # plot the results
  # Model 1
  # Silicate ~ Silicate
  # R<-conditional_effects(fit_brms, "Salinity", resp = "SilicateumolL",  resolution = 1000)
  # R1<-R$SilicateumolL.SilicateumolL_Salinity %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(ModelData$Silicate_umolL,"scaled:scale")+attr(ModelData$Silicate_umolL,"scaled:center"),
  #          lower = lower__*attr(ModelData$Silicate_umolL,"scaled:scale")+attr(ModelData$Silicate_umolL,"scaled:center"),
  #          upper = upper__*attr(ModelData$Silicate_umolL,"scaled:scale")+attr(ModelData$Silicate_umolL,"scaled:center")) %>%
  #   mutate(Salinity = Salinity*attr (ModelData$Salinity,"scaled:scale")+attr(ModelData$Salinity,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = Salinity, y = exp(estimate)), lwd = 1, color = 'grey')+
  #   geom_ribbon(aes(x = Salinity,ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.3, fill = "grey")+
  #   geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Salinity, y = Silicate_umolL)) +
  #   xlab("Salinity (psu)")+
  #   ylab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
  #   ggtitle("Model 1")+
  #   coord_trans(y="log")+
  #  # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  # #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  ## Model 2
  # NN ~ Silicate
  R<-conditional_effects(fit_brms, "SilicateumolL", resp = "NNumolL",  resolution = 1000)
  R2<-R$NNumolL.NNumolL_SilicateumolL %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NNumolL,"scaled:scale")+attr(ModelData$NNumolL,"scaled:center"),
           lower = lower__*attr(ModelData$NNumolL,"scaled:scale")+attr(ModelData$NNumolL,"scaled:center"),
           upper = upper__*attr(ModelData$NNumolL,"scaled:scale")+attr(ModelData$NNumolL,"scaled:center")) %>%
    mutate(SilicateumolL = SilicateumolL*attr (ModelData$SilicateumolL,"scaled:scale")+attr(ModelData$SilicateumolL,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(SilicateumolL), y = exp(estimate)), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = exp(SilicateumolL),ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.3, fill = "grey")+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = NN_umolL)) +
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
    ggtitle("Model 1")+
    coord_trans(x = "log", y="log")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  
  # Model 3
  # NEP ~ DayNight*(Temperature + logNN)
  R<-conditional_effects(fit_brms, "NNumolL:Day.Night", resp = "NEPproxy", resolution = 1000)
  R3<-R$`NEPproxy.NEPproxy_NNumolL:Day.Night` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           lower = lower__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           upper = upper__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")) %>%
    mutate(NNumolL = NNumolL*attr (ModelData$NNumolL,"scaled:scale")+attr(ModelData$NNumolL,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(NNumolL), y = estimate, lty = Day.Night), color = "grey", lwd = 1)+
    geom_ribbon(aes(x = exp(NNumolL),ymin=lower, ymax=upper, lty = Day.Night), fill = "grey", alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NN_umolL, y = NEP.proxy, color = Temperature, shape = Day_Night)) +
    ylab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
    ggtitle("Model 2")+
    coord_trans(x = "log")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    scale_color_gradient(low = "blue", high = "red", name = "Temperature")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
           lty = "none",
         shape = guide_legend(title = "Day/Night", title.position="top", title.hjust = 0.5)) 
 

   # NEP ~ DayNight*(Temperature + logNN)
  R<-conditional_effects(fit_brms, "Temperature", resp = "NEPproxy", resolution = 1000)
  R3a<-R$`NEPproxy.NEPproxy_Temperature` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           lower = lower__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           upper = upper__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")) %>%
    mutate(Temperature = Temperature*attr (ModelData$Temperature,"scaled:scale")+attr(ModelData$Temperature,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = Temperature, y = estimate), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = Temperature,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Temperature, y = NEP.proxy, color = NN_umolL)) +
    ylab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab(expression(atop("Temperature", paste("(",degree, "C)"))))+
    ggtitle("Model 2")+
    scale_color_gradient(name = "N+N", trans = "log", breaks =c(0.1,.2, .3, 0.5,1), low = "lightgreen", high = "darkgreen")+
      # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) 
  
 # model 4
  # pH ~ NEP + log(Silicate)
  R<-conditional_effects(fit_brms, "NEP.proxy", resp = "pH", resolution = 1000)
  R4<-R$pH.pH_NEP.proxy %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           lower = lower__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           upper = upper__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center")) %>%
    mutate(NEP.proxy = NEP.proxy*attr (ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NEP.proxy, y = estimate), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = NEP.proxy,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NEP.proxy, y = pH, color = Silicate_umolL)) +
    scale_color_gradient(name = "Silicate", trans = "log", breaks =c(0.1,1,2,3,5,10,15), low = "lightblue", high = "darkblue")+
    xlab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    ylab("pH")+
    ggtitle("Model 3")+
      
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
      guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) 
  
  # pH ~ NEP + log(Silicate)
  
R<-conditional_effects(fit_brms, "SilicateumolL", resp = "pH", resolution = 1000, re_formula = NULL)

  R4a<-R$`pH.pH_SilicateumolL` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           lower = lower__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           upper = upper__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center")) %>%
    mutate(SilicateumolL = SilicateumolL*attr (ModelData$SilicateumolL,"scaled:scale")+attr(ModelData$SilicateumolL,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(SilicateumolL), y = estimate), lwd = 1, color = "grey")+
    geom_ribbon(aes(x = exp(SilicateumolL),ymin=lower, ymax=upper), fill = "grey", linetype=1.5, alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = pH, color = NEP.proxy)) +
    scale_color_gradient2(low = "#D8B365",
                          mid = "gray88",
                          high = "#097969",
                          midpoint = 0)+
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab("pH")+
    ggtitle("Model 3")+
    coord_trans(x="log")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) 
  
  # Model 5
  # NEC ~ pH + Temperature

  R<-conditional_effects(fit_brms, "pH", resp = "NECproxy", resolution = 1000)
  R5<-R$NECproxy.NECproxy_pH %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center"),
           lower = lower__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center"),
           upper = upper__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center")) %>%
    mutate(pH = pH*attr (ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = pH, y = estimate), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = pH,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = pH, y = NEC.proxy, color = Temperature)) +
    ylab(expression(atop("NEC", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab("pH")+
    ggtitle("Model 4")+
    scale_color_gradient(low = "blue", high = "red", name = "Temperature")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) 
  
  # NEC ~ pH + Temperature
  
  R<-conditional_effects(fit_brms, "Temperature", resp = "NECproxy", resolution = 1000)
  R5a<-R$NECproxy.NECproxy_Temperature %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center"),
           lower = lower__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center"),
           upper = upper__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center")) %>%
    mutate(Temperature = Temperature*attr (ModelData$Temperature,"scaled:scale")+attr(ModelData$Temperature,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = Temperature, y = estimate), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = Temperature,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Temperature, y = NEC.proxy, color = pH)) +
    ylab(expression(atop("NEC", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab(expression(atop("Temperature", paste("(",degree, "C)"))))+
    scale_color_gradient(low = "peachpuff", high = "sienna", name = expression("pH"[t]))+
    ggtitle("Model 4")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) 
  
  # Model 6
  # Humics~NEP.proxy*Day.Night+Silicate_umolL
  
  R<-conditional_effects(fit_brms, "NEP.proxy:Day.Night", resp = "Humics", resolution = 1000)
  R6<-R$`Humics.Humics_NEP.proxy:Day.Night` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           lower = lower__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           upper = upper__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center")) %>%
    mutate(NEP.proxy = NEP.proxy*attr (ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NEP.proxy, y = estimate, lty = Day.Night), lwd = 1)+
    geom_ribbon(aes(x = NEP.proxy,ymin=lower, ymax=upper, lty = Day.Night),fill = "grey", alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NEP.proxy, y = VisibleHumidic_Like +MarineHumic_Like, color = Silicate_umolL, shape = Day_Night)) +
    scale_color_gradient(name = "Silicate", trans = "log", breaks =c(0.1,1,2,3,5,10,15), low = "lightblue", high = "darkblue")+
    xlab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    ylab(expression(atop("Humics","(RU)")))+
    ggtitle("Model 5")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
           lty = "none",
           shape = guide_legend(title = "Day/Night", title.position="top", title.hjust = 0.5)) 
  
  # Humics~NEP.proxy*Day.Night+Silicate_umolL
  R<-conditional_effects(fit_brms, "SilicateumolL", resp = "Humics", resolution = 1000)
  R6a<-R$Humics.Humics_SilicateumolL %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           lower = lower__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           upper = upper__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center")) %>%
    mutate(SilicateumolL = SilicateumolL*attr (ModelData$SilicateumolL,"scaled:scale")+attr(ModelData$SilicateumolL,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(SilicateumolL), y = estimate), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = exp(SilicateumolL),ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = VisibleHumidic_Like +MarineHumic_Like, color = NEP.proxy)) +
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab(expression(atop("Humics","(RU)")))+
    ggtitle("Model 5")+
    coord_trans(x="log")+
    scale_color_gradient2(low = "#D8B365",
                          mid = "gray88",
                          high = "#097969",
                          midpoint = 0)+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))  
 
  # Model 7
  # bf(Proteinaceous~NEC.proxy+Day.Night)
  R<-conditional_effects(fit_brms, "NEC.proxy", resp = "Proteinaceous", resolution = 1000)
  R7<-R$`Proteinaceous.Proteinaceous_NEC.proxy` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Proteinaceous,"scaled:scale")+attr(ModelData$Proteinaceous,"scaled:center"),
           lower = lower__*attr(ModelData$Proteinaceous,"scaled:scale")+attr(ModelData$Proteinaceous,"scaled:center"),
           upper = upper__*attr(ModelData$Proteinaceous,"scaled:scale")+attr(ModelData$Proteinaceous,"scaled:center")) %>%
    mutate(NEC.proxy = NEC.proxy*attr (ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NEC.proxy, y = estimate), color = "grey", lwd = 1)+
    geom_ribbon(aes(x = NEC.proxy,ymin=lower, ymax=upper),fill = "grey", linetype=1.5, alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NEC.proxy, y = Tryptophan_Like +Tyrosine_Like)) +
    xlab(expression(atop("NEC", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    ylab(expression(atop("Proteinaceous","(RU)")))+
    ggtitle("Model 7")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
    ## bring them all together in patchwork
  R<-(R2|R7)/(R3|R3a)/(R4|R4a)/(R5|R5a)/(R6|R6a)+
    plot_annotation(tag_levels = "A")+
    #/(R6|R6b)+
    plot_layout(guides = "collect")&
    guides(shape = guide_legend(order = 1, title.position = "top", title = "Day/Night",title.hjust = 0.5))&
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 11))
  
  ggsave(here("Output",paste(site,season,"marginaleffects.pdf")),R, width = 12, height = 18, useDingbats = FALSE)
  
  var_name <- paste(fit_brms, site,season, sep="_") # Construct the name
  var_name <- str_replace_all(var_name, " ", "")   # remove white space
  return(fit_brms)
  
}

# Run the SEM
V_Dry_fit<-RunSEM(site ="Varari", season = "Dry")
V_Wet_fit<-RunSEM(site ="Varari", season = "Wet")
C_Dry_fit<-RunSEM(site ="Cabral", season = "Dry")
C_Wet_fit<-RunSEM(site ="Cabral", season = "Wet")

## Get the posterior
## get the posterior
get_post<-function(site,season, fit){
  
  post <- as_draws_df(fit)
  
  Cof<-post %>% 
    select(starts_with("b"),-ends_with("Intercept")) %>%
    gather() %>% 
    group_by(key)%>%
    median_hdci()%>%
    mutate(sig = ifelse(sign(.lower)==sign(.upper),'yes','no'))%>%# if not significant make it grey
    separate(col = key,into = c("b", "dependent", "independent"),sep = "_")%>% #loose the b and bring the values back together
    mutate(independent  = recode(independent, Day = "Day/Night", `NEP.proxy:Day` = "NEP.proxy:Day/Night")) %>%# mutate(key = paste(dependent, independent))%>%
    mutate(dependent = factor(dependent, levels = c("NNumolL","NEPproxy","pH","NECproxy","Humics","Proteinaceous")))%>%
    mutate(Location = site, Season = season)
  
  return(Cof)
}

VDryCoF<-get_post(site = "Varari",season = "Dry", fit = V_Dry_fit)
VWerCoF<-get_post(site = "Varari",season = "Wet", fit = V_Wet_fit)
CDryCoF<-get_post(site = "Cabral",season = "Dry", fit = C_Dry_fit)
CWetCoF<-get_post(site = "Cabral",season = "Wet", fit = C_Wet_fit)

AllCoefs <- VDryCoF %>%
  bind_rows(VWerCoF,CDryCoF,CWetCoF)

# Make plots
IndivPlots<-function(site, season){
  AllCoefs%>%
    filter(Location == site, Season == season) %>%
    ggplot(aes(x = value, y = independent, shape = sig)) +  # note how we used `reorder()` to arrange the coefficients
    geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
    geom_point(size = 2)+
    geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
    #scale_alpha_manual(values = c(0.2,1))+
    scale_shape_manual(values = c(1,16))+
    #scale_color_brewer(palette = "Set2")+
    scale_color_manual(values = c("firebrick4"), name = " ")+
    labs(title = paste(site, season),
         x = NULL, y = NULL) +
    theme_bw() +
    guides(shape = "none")+
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
          axis.text.y  = element_text(hjust = 0),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          legend.text=element_text(size=14),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold")
    )+
    facet_grid(~dependent, scales = "free_y", space='free')
}

IndivPlots("Varari","Dry")
IndivPlots("Varari","Wet")
IndivPlots("Cabral","Wet")
IndivPlots("Cabral","Dry")

### Plot all together
#Make the plot
CoefPlot<-AllCoefs%>%
  # ggplot(aes(x = value, y = reorder(independent, value), shape = sig, color = Site)) +  # note how we used `reorder()` to arrange the coefficients
  ggplot(aes(x = value, y = independent, shape = sig, color = Location)) +  # note 
  geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
  geom_point(size = 2)+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
  #scale_alpha_manual(values = c(0.2,1))+
  scale_shape_manual(values = c(1,16))+
  scale_color_brewer(palette = "Set2")+
  #  scale_color_manual(values = c("firebrick4", "orange"), name = " ")+
  
  labs(
    x = NULL, y = NULL) +
  theme_bw() +
  guides(shape = FALSE)+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold")
  )+
  facet_grid(Season~dependent, scales = "free_y", space='free')
ggplot2::ggsave("Output/coefficientsAll.pdf", width = 10, height = 8, useDingbats = FALSE)


## Plot showing nitrate leading to net heterotrophy
Cdata %>% filter(Plate_Seep == "Plate") %>%
  group_by(Location, CowTagID, Season) %>%
  summarise(meanNEP = mean(NEP.proxy, na.rm = TRUE),
            sumNEP = sum(NEP.proxy),
            meanNN = mean(NN_umolL),
            meanSi = mean(Silicate_umolL),
            maxSi = max(Silicate_umolL, na.rm = TRUE)) %>%
  drop_na()%>%
  ggplot(aes(x = maxSi, y = sumNEP, color = Season, fill = Season))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", formula = "y~log(x)")+
  coord_trans(x = "log") +
  facet_wrap(~Location, scales = "free_x")+
  xlab(expression(atop("Max Silicate", paste("(",mu,"mol L"^-1,")"))))+
  ylab(expression(atop("Integrated NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
  scale_color_manual(values = c("#800000","#767676"))+
  scale_fill_manual(values = c("#800000","#767676"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"))
ggsave(here("Output","heterotrophic.pdf"),width = 7, height = 5)
