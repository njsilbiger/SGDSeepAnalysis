### Bayesian  model showing the relationship between SGD and PR by species interactions for Jamie Kerlin
### Created on 7/27/2022
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

Data<-read_csv(here("Data","PR_turb_biogeochem_wide_2022_06_28.csv"))

## Just use a subset of the data to test the code

Data_sub<-Data %>%
  filter(Response_measurement == "GP") # Let's just look at GP

# extrat the data for the controls only
Data_sub_control<-Data %>%
  filter(Response_measurement == "GP", Treatment == "Monoculture")
#####################3
# Data list for GP ~ pHmin with residuals for competition treatment
stan_data <- list(N=length(Data_sub_control$Treatment), # the length of the control dat 
                  x = Data_sub_control$minpH, # use pH min for this example as the independent variable
                  x2 = Data_sub_control$Initial_PR, # Covariate
                  y = Data_sub_control$Final_PR, # Final pR is the y 
                  Treatment  = as.numeric(as.factor(Data_sub$Treatment)), # make the treatment numbers 1-4
                  n_treatment = 4, # number of treatments
                  PlateID = Data_sub$PlateID, # plate ID for the random intercepts
                  n_resid = length(Data_sub$PlateID), # the length of the data to calculate residuals from
                  AllPR = Data_sub$Final_PR, # PR from all treatments
                  AllInitial = Data_sub$Initial_PR,  # all initial 
                  AllpHmin = Data_sub$minpH, # All pH min
                  n_PlateID = length(unique(Data_sub$PlateID))# number of plates
                  )

### Write the STAN Model 
write("// Stan model for simple linear regression with residuals

data {
 int < lower = 1 > N; // Sample size for control regression
 int < lower = 1 > n_treatment; // number of treatments
 int < lower = 1 > n_PlateID; // number of plates
 vector[N] y; // Response
 vector[N] x; // Predictor
 vector[N] x2; // Covariate
 int < lower = 1 > n_resid; // Sample size for residual model
 int<lower = 1, upper = n_treatment> Treatment[n_resid]; // Treatment variable
 vector[n_resid] PlateID; // PlateID for random effect of plate
 vector[n_resid] AllPR; // All final PR
 vector[n_resid] AllInitial; // All initial PR
 vector[n_resid] AllpHmin; // All pHmin
 

}

parameters {
 real alpha; // Intercept
 real beta; // Slope (regression coefficients for x)
 real beta2; // Slope for x2
 real < lower = 0 > sigma; // Error SD
 real < lower = 0 > sigma_resid; // Error SD for residual model
 real beta_resid;
 //vector[n_PlateID] alpha_plate;
 real alpha_plate;

}

model {

// conditional mean
vector[N] mu;

// linear combination
 mu = alpha + x * beta + x2*beta2;

// likelihood function
 y ~ normal(mu, sigma);
 }

generated quantities {
vector[N] y_rep;
vector[n_resid] PR_expected;
vector[n_resid] PR_Resid;
vector[n_resid] y_resid;


 for (i in 1:N) {
 // generate predictions
 real y_hat = alpha + x[i] * beta + x2[i] *beta2; // Posterior predictions

 // generate replication values
 y_rep[i] = normal_rng(y_hat, sigma);
} // The posterior predictive distribution

// calculate predicted TA based on mixing from plate data
for (i in 1: n_resid){
 PR_expected[i] = alpha + AllpHmin[i] * beta + AllInitial[i]* beta2; // calculate the expected value
 PR_Resid[i] = AllPR[i] - PR_expected[i]; // calculate the residuals

 
}


 
}
  
model {
// conditional mean
 real mu_resid;


 // generate model for residuals ~ treatment with plateID as random
 
 mu_resid = alpha_plate + beta_resid*Treatment[i];
 
 // for (j in 1: n_PlateID){
 // mu_resid = alpha_plate[PlateID[j]] + beta_resid*PR_Resid[i];
 // }
 
 PR_Resid ~ normal(mu_resid, sigma_resid);

}



}",

"stan_modelJamie.stan")

# Run the model
stan_model <- "stan_modelJamie.stan"

# run the model --  increase the interations for better fit
fit2 <- stan(file = stan_model, data = stan_data, warmup = 1000, iter = 3000, chains = 3, cores = 2, thin = 2,
             seed = 11)

# assess the fit
posterior <- extract(fit2)

# Posterior predictive checks for pH min model
y_rep <- as.matrix(fit2, pars = "y_rep")
y_rep_sum <- data.frame(summary(fit2, pars = "y_rep")$summary)
ppc_dens_overlay(Data_sub_control$Final_PR, y_rep[1:200, ])

# Plot the fit line with 95% Bayes prediction intervals

Data_sub_control %>%
  mutate(mean_y_rep = y_rep_sum$mean,
         lower = y_rep_sum$X2.5.,
         upper = y_rep_sum$X97.5.)%>%
  ggplot() + 
  geom_point(aes(minpH, Final_PR)) + 
  geom_line(aes(minpH, mean_y_rep), size = 1.5) + 
  geom_ribbon(aes(minpH, ymin = lower, ymax = upper), alpha = 0.35, lty =2)+
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()


# make a boxplot of the residuals
y_resid <- as.matrix(fit2, pars = "y_resid")
y_resid_sum <- data.frame(summary(fit2, pars = "y_resid")$summary)

Data_sub_control %>%
  mutate(mean_y_resid = y_rep_sum$mean) %>%
  ggplot()+
  geom_boxplot(aes(x = Treatment, y =mean_y_resid ))