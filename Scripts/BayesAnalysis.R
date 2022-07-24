### Bayesian SEM model showing the relationship between SGD and ecosystem functioning library(tidybayes)
### Created on 7/22/2022
### Created by Nyssa Silbiger

### Load Libraries ##########
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(brms)
library(posterior)
library(modelr)
library(rstan)
library(ggfortify)
library(lubridate)
library(ggforce)
library(viridis)
library(patchwork)
library(here)
library(ggridges)
library(seacarb)
library(ggtext)
library(lme4)
library(lmerTest)
library(broom)


##### Load data ###########
# load the 24 hour chemistry data #####################
Data_Dry<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC2.csv") %>%
  mutate(Season = "Dry") %>%
  mutate(Date = mdy(Date))
turbdata<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Nutrients/Turb_NC.csv") %>%
  mutate(Season = "Dry")
# wet season
Data_wet<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/CarbonateChemistry/pHProbe_Data_calculated_POcorrect.csv") %>%
  mutate(Season = "Wet") %>%
  rename(Temperature = TempInSitu) %>%
  left_join(read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/Sandwich_Locations_Final.csv"))

turb_wet<- read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/Nutrients/Turb_NC.csv") %>%
  mutate(Season = "Wet")


##### Clean Data ##########
remove2<-Data_Dry %>% filter(CowTagID=="V2", Tide =="Low", Day_Night=="Day", Date == ymd("2021-08-08"))
removelow<- Data_Dry %>% # remove the not real low tide
  filter(Date == ymd("2021-08-06") & Tide == "Low" & Plate_Seep == "Plate")


# bring both seasons together
Data <- Data_Dry %>%
  bind_rows(Data_wet)

### Calculate all carbonate parameters #############
#remove rows with NAs
Cdata<-Data %>%
  select(-c(Jamie_Plate_ID, Bottom_Plate_ID, Top_Plate_ID)) %>%
  mutate(Tide_Time = paste(Tide, Day_Night)) %>%
  drop_na(TA) # keep only complete cases


#calculate rest of carbonate params
CO2<-carb(flag=8, Cdata$pH, Cdata$TA/1000000, S=Cdata$Salinity, 
          T=Cdata$Temperature, Patm=1, P=0, Pt=Cdata$Phosphate_umolL/1000000, Sit=Cdata$Silicate_umolL/1000000, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")

#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propogation
er<-errors(flag=8, Cdata$pH, Cdata$TA/1000000, 
           S=Cdata$Salinity, T=Cdata$Temperature, 
           Patm=1, P=0,Pt=Cdata$Phosphate_umolL/1000000,
           Sit=Cdata$Silicate_umolL/1000000,evar1 = 0.01, evar2 = 5e-6) 

#average error for DIC based on pH and TA
mean(er$DIC*1000000)
sd(er$DIC*1000000)/sqrt(nrow(er))

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
CO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-CO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

Cdata[,c("CO2","HCO3","CO3","DIC","OmegaArag","OmegaCalcite","pCO2","fCO2")]<-
  CO2[,c("CO2","HCO3","CO3","DIC","OmegaAragonite","OmegaCalcite","pCO2","fCO2")]

# Create a TA mixing line the frequentist way--- think of maybe doing one for each season for better accuarcy
VarariMixModel<-lmer(TA~Salinity+(1|Season), data = Cdata %>% filter(Location == "Varari", Plate_Seep == "Seep"))
Vco<-coef(VarariMixModel)

CabralMixModel<-lmer(TA~Salinity+(1|Season), data = Cdata %>% filter(Location == "Cabral", Plate_Seep == "Seep"))
Cco<-coef(CabralMixModel)

# DIC mixing line
VarariMixModelDIC<-lmer(DIC~Salinity+(1|Season), data = Cdata %>% filter(Location == "Varari", Plate_Seep == "Seep"))
VcoDIC<-coef(VarariMixModelDIC)

CabralMixModelDIC<-lmer(DIC~Salinity+(1|Season), data = Cdata %>% filter(Location == "Cabral", Plate_Seep == "Seep"))
CcoDIC<-coef(CabralMixModelDIC)


## put them in the dataframe
Cdata <- Cdata %>% # add the predicted mixing line
  mutate(TA.mix = case_when(Location  == "Varari" & Season == "Dry" ~Salinity*Vco$Season[1,2]+Vco$Season[1,1],
                            Location  == "Varari" & Season == "Wet" ~Salinity*Vco$Season[2,2]+Vco$Season[2,1],
                            Location  == "Cabral" & Season == "Dry" ~Salinity*Cco$Season[1,2]+Cco$Season[1,1],
                            Location  == "Cabral" & Season == "Wet" ~Salinity*Cco$Season[2,2]+Cco$Season[2,1]
  )) %>%
  mutate(DIC.mix = case_when(Location  == "Varari" & Season == "Dry" ~Salinity*VcoDIC$Season[1,2]+VcoDIC$Season[1,1],
                             Location  == "Varari" & Season == "Wet" ~Salinity*VcoDIC$Season[2,2]+VcoDIC$Season[2,1],
                             Location  == "Cabral" & Season == "Dry" ~Salinity*CcoDIC$Season[1,2]+CcoDIC$Season[1,1],
                             Location  == "Cabral" & Season == "Wet" ~Salinity*CcoDIC$Season[2,2]+CcoDIC$Season[2,1]
  )) %>%
  mutate(
    TA.diff = TA.mix- TA, # TA expected by mixing alone - TA at plate
    DIC.diff = DIC.mix - DIC)
   
    #### Run Analysis #####

## Calculate mixing line the Bayesian Way #########

# extract just the Varari Data from the seep to create the mxing line
SeepDataVarari<-Cdata %>% 
  filter(Location == "Varari", Plate_Seep == "Seep") # Pull out just the seep data

# plate data
PlateDataVarari<-Cdata %>% 
  filter(Location == "Varari", Plate_Seep == "Plate") # Pull out just 

### run the model now where we calculate the difference between the TA at the plate and the expected TA due to mixing based on salinity

#####################3
# Data list
stan_data <- list(N=length(SeepDataVarari$TA), 
                  x = SeepDataVarari$Salinity,
                  y = SeepDataVarari$TA, 
                  season  = ifelse(SeepDataVarari$Season == "Dry",1,2),                   n_season = 2,
                  TAplate = PlateDataVarari$TA,
                  Salplate = PlateDataVarari$Salinity,
                  n_plate = length(PlateDataVarari$TA),
                  season_plate = ifelse(PlateDataVarari$Season == "Dry",1,2),
                  y_DIC =SeepDataVarari$DIC,
                  DICplate = PlateDataVarari$DIC
                  )

### Write the STAN Model 
write("// Stan model for simple linear regression

data {
 int < lower = 1 > N; // Sample size
 int < lower = 1 > n_season; // number of seasons
 vector[N] x; // Predictor
 vector[N] y; // Outcome
 vector[N] y_DIC; // Outcome
 int<lower = 0, upper = n_season> season [N]; // season prediction
 int < lower = 1> n_plate; // sample size for plates
 vector[n_plate] Salplate; // Salinity from the plate
 vector[n_plate] TAplate; // TA from the plate
 vector[n_plate] DICplate; // DIC from the plate
 int <lower = 0, upper = n_season> season_plate [n_plate]; // season prediction

}

parameters {
 vector[n_season] alpha; // vector of intercepts
 real beta; // Slope (regression coefficients)
 real < lower = 0 > sigma; // Error SD

// DIC params
 vector[n_season] alpha_DIC; // vector of intercepts
 real beta_DIC; // Slope (regression coefficients)
 real < lower = 0 > sigma_DIC; // Error SD
}

model {

// conditional mean
vector[N] mu;
vector[N] mu_DIC;

// linear combination
 mu = alpha[season] + x * beta;
 mu_DIC = alpha_DIC[season] + x * beta_DIC;

// likelihood function
 y ~ normal(mu, sigma);
 y_DIC ~ normal(mu_DIC, sigma_DIC);
}

generated quantities {
vector[N] y_rep;
vector[N] y_rep_DIC;
vector[n_plate] TA_pred;
vector[n_plate] TA_diff;
vector[n_plate] DIC_pred;
vector[n_plate] DIC_diff;

 for (i in 1:N) {
 // generate predictions
 real y_hat = alpha[season[i]] + x[i] * beta;
 real y_hat_DIC = alpha_DIC[season[i]] + x[i] * beta_DIC;
 
 // generate replication values
 y_rep[i] = normal_rng(y_hat, sigma);
 y_rep_DIC[i] = normal_rng(y_hat_DIC, sigma_DIC);
} // The posterior predictive distribution

// calculate predicted TA based on mixing from plate data
for (i in 1: n_plate){
 TA_pred[i] = alpha[season_plate[i]] + Salplate[i] * beta;
 TA_diff[i] = TA_pred[i] - TAplate[i];
 DIC_pred[i] = alpha_DIC[season_plate[i]] + Salplate[i] * beta_DIC;
 DIC_diff[i] = DIC_pred[i] - DICplate[i];


  }

}",

"stan_model2.stan")

stan_model2 <- "stan_model2.stan"

# run the model
fit2 <- stan(file = stan_model2, data = stan_data, warmup = 1000, iter = 3000, chains = 3, cores = 2, thin = 2,
             pars = c('beta', 'sigma', 'alpha', 'TA_diff', 'y_rep','beta_DIC','alpha_DIC','y_rep_DIC', 'DIC_diff'), seed = 11)

# assess the fit
posterior <- extract(fit2)
#stan_dens(fit)

## extract the TA_diff values
TA_diff<-summary(fit2, pars = 'TA_diff')$summary

## extract the DIC_diff values
DIC_diff<-summary(fit2, pars = 'DIC_diff')$summary

## look at the relationship between bayesian and frequentist TA diff values
bind_cols(PlateDataVarari, TA_diff) %>%
  ggplot()+
  geom_point(aes(x = TA.diff, y = mean,color = Season))+
  xlab("TA difference from lmer")+
  ylab("TA difference from STAN")

## look at the relationship between bayesian and frequentist TA diff values
bind_cols(PlateDataVarari, DIC_diff) %>%
  ggplot()+
  geom_point(aes(x = DIC.diff, y = mean,color = Season, fill = Season))+
  xlab("DIC difference from lmer")+
  ylab("DIC difference from STAN")

# Posterior predictive checks for TA
y_rep <- as.matrix(fit2, pars = "y_rep")
y_rep_sum <- data.frame(summary(fit2, pars = "y_rep")$summary)
ppc_dens_overlay(SeepDataVarari$TA, y_rep[1:200, ])

# Posterior predictive checks for DIC
y_rep_DIC <- as.matrix(fit2, pars = "y_rep_DIC")
y_rep_sum_DIC <- data.frame(summary(fit2, pars = "y_rep_DIC")$summary)
ppc_dens_overlay(SeepDataVarari$DIC, y_rep_DIC[1:200, ])


# Plot the mixing line with 95% Bayes prediction intervals
# TA
SeepDataVarari %>%
  mutate(mean_y_rep = y_rep_sum$mean,
         lower = y_rep_sum$X2.5.,
         upper = y_rep_sum$X97.5.)%>%
  ggplot(aes(color = Season, fill = Season)) + 
  geom_point(aes(Salinity, TA)) + 
  geom_line(aes(Salinity, mean_y_rep), size = 1.5) + 
  geom_ribbon(aes(Salinity, ymin = lower, ymax = upper), alpha = 0.35, lty =2)+
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()

# DIC
SeepDataVarari %>%
  mutate(mean_y_rep = y_rep_sum_DIC$mean,
         lower = y_rep_sum_DIC$X2.5.,
         upper = y_rep_sum_DIC$X97.5.)%>%
  ggplot(aes(color = Season, fill = Season)) + 
  geom_point(aes(Salinity, DIC)) + 
  geom_line(aes(Salinity, mean_y_rep), size = 1.5) + 
  geom_ribbon(aes(Salinity, ymin = lower, ymax = upper), alpha = 0.35, lty =2)+
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()



