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
  mutate_at(vars("pH":"Proteinaceous"),  function(x) scale(x, center = TRUE)) 

### Create the models ####

SGDmod<-bf(Silicate_umolL~Salinity)
NNmod<-bf(NN_umolL ~Silicate_umolL) ### N and P are highly correlated 
Pmod<-bf(Phosphate_umolL ~Silicate_umolL)
NH4mod<-bf(Ammonia_umolL ~Silicate_umolL) # N and NH4 are weakly correlated
NEPmod<-bf(NEP.proxy~Day_Night*(NN_umolL) +Temperature) ## N or P or NH4?
pHmod<-bf(pH~ NEP.proxy+Silicate_umolL)
NECmod<-bf(NEC.proxy~pH+Temperature)
Humicsmod<-bf(Humics~NEP.proxy*Day_Night+Silicate_umolL)
Protmod<-bf(Proteinaceous~NEC.proxy+Day_Night)

# Function to run Bayesian SEM and make the posterior predictive checks and plot marginal effects
RunSEM<-function(site, season){
  
  
  fit_brms <- brm(SGDmod+
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

  p1+ p2+p3+p4+p5+p6+p7+
    #p6+
    plot_layout(guides = "collect") +
    plot_annotation(title = 'Posterior Predictive Checks', tag_levels = "A")
  ggsave(here("Output",paste(site,season, "Posteriorchecks.pdf")), width = 5, height = 5)
  
  # plot the results
  # Model 1
  # Silicate ~ Silicate
  R<-conditional_effects(fit_brms, "Salinity", resp = "SilicateumolL",  resolution = 1000)
  R1<-R$SilicateumolL.SilicateumolL_Salinity %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Silicate_umolL,"scaled:scale")+attr(ModelData$Silicate_umolL,"scaled:center"),
           lower = lower__*attr(ModelData$Silicate_umolL,"scaled:scale")+attr(ModelData$Silicate_umolL,"scaled:center"),
           upper = upper__*attr(ModelData$Silicate_umolL,"scaled:scale")+attr(ModelData$Silicate_umolL,"scaled:center")) %>%
    mutate(Salinity = Salinity*attr (ModelData$Salinity,"scaled:scale")+attr(ModelData$Salinity,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = Salinity, y = exp(estimate)), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = Salinity,ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.3, fill = "grey")+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Salinity, y = Silicate_umolL)) +
    xlab("Salinity (psu)")+
    ylab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ggtitle("Model 1")+
    coord_trans(y="log")+
   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  ## Model 2
  # NN ~ Silicate
  R<-conditional_effects(fit_brms, "Silicate_umolL", resp = "NNumolL",  resolution = 1000)
  R2<-R$NNumolL.NNumolL_Silicate_umolL %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NN_umolL,"scaled:scale")+attr(ModelData$NN_umolL,"scaled:center"),
           lower = lower__*attr(ModelData$NN_umolL,"scaled:scale")+attr(ModelData$NN_umolL,"scaled:center"),
           upper = upper__*attr(ModelData$NN_umolL,"scaled:scale")+attr(ModelData$NN_umolL,"scaled:center")) %>%
    mutate(Silicate_umolL = Silicate_umolL*attr (ModelData$Silicate_umolL,"scaled:scale")+attr(ModelData$Silicate_umolL,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(Silicate_umolL), y = exp(estimate)), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = exp(Silicate_umolL),ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.3, fill = "grey")+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = NN_umolL)) +
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
    ggtitle("Model 2")+
    coord_trans(x = "log", y="log")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  
  # Model 3
  # NEP ~ DayNight*(Temperature + logNN)
  R<-conditional_effects(fit_brms, "NN_umolL:Day_Night", resp = "NEPproxy", resolution = 1000)
  R3<-R$`NEPproxy.NEPproxy_NN_umolL:Day_Night` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           lower = lower__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           upper = upper__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")) %>%
    mutate(NN_umolL = NN_umolL*attr (ModelData$NN_umolL,"scaled:scale")+attr(ModelData$NN_umolL,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(NN_umolL), y = estimate, color = Day_Night), lwd = 1)+
    geom_ribbon(aes(x = exp(NN_umolL),ymin=lower, ymax=upper, fill = Day_Night), linetype=1.5, alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NN_umolL, y = NEP.proxy, color = Day_Night)) +
    ylab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
    ggtitle("Model 3")+
    coord_trans(x = "log")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

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
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Temperature, y = NEP.proxy)) +
    ylab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab(expression(atop("Temperature", paste("(",degree, "C)"))))+
    ggtitle("Model 3")+
      # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
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
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NEP.proxy, y = pH)) +
    xlab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    ylab("pH")+
    ggtitle("Model 4")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # pH ~ NEP + log(Silicate)
  R<-conditional_effects(fit_brms, "Silicate_umolL", resp = "pH", resolution = 1000)
  R4a<-R$pH.pH_Silicate_umolL %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           lower = lower__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           upper = upper__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center")) %>%
    mutate(Silicate_umolL = Silicate_umolL*attr (ModelData$Silicate_umolL,"scaled:scale")+attr(ModelData$Silicate_umolL,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(Silicate_umolL), y = estimate), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = exp(Silicate_umolL),ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = pH)) +
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab("pH")+
    ggtitle("Model 4")+
    coord_trans(x="log")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom') 
  
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
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = pH, y = NEC.proxy)) +
    ylab(expression(atop("NEC", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab("pH")+
    ggtitle("Model 5")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
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
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Temperature, y = NEC.proxy)) +
    ylab(expression(atop("NEC", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab(expression(atop("Temperature", paste("(",degree, "C)"))))+
    ggtitle("Model 5")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 6
  # Humics~NEP.proxy*Day_Night+Silicate_umolL
  
  R<-conditional_effects(fit_brms, "NEP.proxy:Day_Night", resp = "Humics", resolution = 1000)
  R6<-R$`Humics.Humics_NEP.proxy:Day_Night` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           lower = lower__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           upper = upper__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center")) %>%
    mutate(NEP.proxy = NEP.proxy*attr (ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NEP.proxy, y = estimate, color = Day_Night), lwd = 1)+
    geom_ribbon(aes(x = NEP.proxy,ymin=lower, ymax=upper, fill = Day_Night), linetype=1.5, alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NEP.proxy, y = VisibleHumidic_Like +MarineHumic_Like, color = Day_Night)) +
    xlab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    ylab(expression(atop("Humics","(RU)")))+
    ggtitle("Model 6")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Humics~NEP.proxy*Day_Night+Silicate_umolL
  R<-conditional_effects(fit_brms, "Silicate_umolL", resp = "Humics", resolution = 1000)
  R6a<-R$Humics.Humics_Silicate_umolL %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           lower = lower__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           upper = upper__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center")) %>%
    mutate(Silicate_umolL = Silicate_umolL*attr (ModelData$Silicate_umolL,"scaled:scale")+attr(ModelData$Silicate_umolL,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(Silicate_umolL), y = estimate), lwd = 1, color = 'grey')+
    geom_ribbon(aes(x = exp(Silicate_umolL),ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = VisibleHumidic_Like +MarineHumic_Like)) +
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab(expression(atop("Humics","(RU)")))+
    ggtitle("Model 7")+
    coord_trans(x="log")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom') 
 
  # Model 7
  # bf(Proteinaceous~NEC.proxy+Day_Night)
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
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NEC.proxy, y = Tryptophan_Like +Tyrosine_Like, color = Day_Night)) +
    xlab(expression(atop("NEC", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    ylab(expression(atop("Proteinaceous","(RU)")))+
    ggtitle("Model 7")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
    ## bring them all together in patchwork
  R<-(R1|R2)/(R3|R3a)/(R4|R4a)/(R5|R5a)/(R6|R6a)/(R7)+
    #/(R6|R6b)+
    plot_layout(guides = "collect")+
    plot_annotation(tag_levels = "A")&
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