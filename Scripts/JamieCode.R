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
library(ggthemes)

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




stanvars <- stanvar(scode = "vector[N] y_rep;
 vector[n_resid] PR_expected;
 vector[n_resid] PR_Resid;

 for (i in 1:N) {
 // generate predictions
 real y_hat = alpha + x[i] * beta + x2[i] *beta2; // Posterior predictions

 // generate replication values
 y_rep[i] = normal_rng(y_hat, sigma);
} // The posterior predictive distribution





// calculate predicted PR based on mixing from plate data
for (i in 1: n_resid){

 // generate model for residuals ~ treatment with plateID as random
 
 PR_expected[i] = alpha + AllpHmin[i] * beta + AllInitial[i]* beta2; // calculate the expected value
 PR_Resid[i] = AllPR[i] - PR_expected[i]; // calculate the residuals",
          block = "genquant")

mod1<-brm(bf(y~x+x2), data = data.frame(x = Data_sub_control$minpH, # use pH min for this example as the independent variable
                                    x2 = Data_sub_control$Initial_PR, # Covariate
                                    y = Data_sub_control$Final_PR), iter = 10000, chains = 3,
          save_pars = save_pars(all = TRUE)) 



#posterio predictive checm
pp_check(mod1, nsamples = 50)

# model validation. Can use loo_compare for model selection to compare models 
loo(mod1,
    moment_match = TRUE)# this deals with problem observations

# Extract posterior values for each parameter
samples1 <- posterior_samples(mod1, "^b")
head(samples1)

ggplot(samples1, aes(x = b_x))+
  geom_density()

pred_draws<-mod1 %>% 
  epred_draws(newdata = expand_grid(x2 = seq(1,2, by = 0.05),
                                    x = seq(7.97, 8.02, by = 0.001)), 
              re_formula = NA)

ggplot(pred_draws, 
       aes(x = x, y = .epred)) +
  stat_lineribbon() +
  geom_point(data =Data_sub_control, aes(x = minpH, y =Final_PR ) )+
  scale_fill_brewer(palette = "Reds") +
  labs(x = "min pH", y = "GP",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")

# Calculate the residuals using the posterior distributions

predicteddata<-matrix(NA, nrow = length(Data_sub$minpH), ncol = nrow(samples1))
residualdata<-matrix(NA, nrow = length(Data_sub$minpH), ncol = nrow(samples1))


for (i in 1:length(Data_sub$minpH)){
  # This is the distribution of all the predicted values
  predicteddata[i,]<-samples1$b_Intercept + samples1$b_x*Data_sub$minpH[i] + samples1$b_x2*Data_sub$Initial_PR[i]
  # This is the distribution of all the residual values. The mean of each row is the mean residual for each coral
  residualdata[i,]<-Data_sub$Final_PR[i]-predicteddata[i,]
}


## calculate the mean and SD for each residual value so that I can carry the SD value to the next model as a prior on y

# calculate means by row
Data_combd<-bind_cols(Data_sub, #add it to the dataset
tibble(residualdata) %>% 
  rowwise() %>% 
  mutate(
    resid = mean(c_across()),
    resid_sd = sd(c_across())
  ) %>%
  select(resid, resid_sd)
) 
  


### brms model 2

dist_fm <- bf(
  resid ~ Treatment + (1|PlateID),
  sigma ~ resid_sd # let sigma vary by the SD of the residuals data
)

dist_brm <- brm(dist_fm, data = Data_combd)


pp_check(dist_brm)

pp_check(dist_brm, type = "scatter_avg_grouped", group = "Treatment") + 
    geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)

#plot(dist_brm)

p1<- conditional_effects(dist_brm) 

plot(p1, plot = FALSE)[[1]] + # plot the effects 
  geom_hline(yintercept = 0)

### Write the STAN Model 
write("// Stan model for simple linear regression with residuals

data {
 int < lower = 1 > N; // Sample size for control regression
 int < lower = 1 > n_treatment; // number of treatments
 // int < lower = 1 > n_PlateID; // number of plates
 vector[N] y; // Response
 vector[N] x; // Predictor
 vector[N] x2; // Covariate
 int < lower = 1 > n_resid; // Sample size for residual model
 // int<lower = 1, upper = n_treatment> Treatment[n_resid]; // Treatment variable
 // vector[n_resid] PlateID; // PlateID for random effect of plate
 vector[n_resid] AllPR; // All final PR
 vector[n_resid] AllInitial; // All initial PR
 vector[n_resid] AllpHmin; // All pHmin
 // vector[n_resid] Treatment; // Treatment
 

}

parameters {
 real alpha; // Intercept
 real beta; // Slope (regression coefficients for x)
 real beta2; // Slope for x2
 real < lower = 0 > sigma; // Error SD
 //  real < lower = 0 > sigma_resid; // Error SD for residual model
 // real beta_resid;
 //vector[n_PlateID] alpha_plate;
 // real alpha_plate;

}


model {

// conditional mean
vector[N] mu;
// vector[n_resid] mu_resid;

// linear combination
 mu = alpha + x * beta + x2*beta2;

// likelihood function
 y ~ normal(mu, sigma);
 
 }

generated quantities {
 vector[N] y_rep;
 vector[n_resid] PR_expected;
 vector[n_resid] PR_Resid;

 for (i in 1:N) {
 // generate predictions
 real y_hat = alpha + x[i] * beta + x2[i] *beta2; // Posterior predictions

 // generate replication values
 y_rep[i] = normal_rng(y_hat, sigma);
} // The posterior predictive distribution





// calculate predicted PR based on mixing from plate data
for (i in 1: n_resid){

 // generate model for residuals ~ treatment with plateID as random
 
// mu_resid = alpha_plate + beta_resid*Treatment[i];
 
 // for (j in 1: n_PlateID){
 // mu_resid = alpha_plate[PlateID[j]] + beta_resid*PR_Resid[i];
 // }
 
//PR_Resid ~ normal(mu_resid, sigma_resid);

 PR_expected[i] = alpha + AllpHmin[i] * beta + AllInitial[i]* beta2; // calculate the expected value
 PR_Resid[i] = AllPR[i] - PR_expected[i]; // calculate the residuals
}



}",

"stan_modelJamie.stan")

# Run the model
stan_model <- "stan_modelJamie.stan"

# run the model --  increase the interations for better fit
fit2 <- stan(file = stan_model, data = stan_data, warmup = 3000, iter = 6000, chains = 3, cores = 2, thin = 2,
             seed = 11)

# assess the fit
posterior <- extract(fit2)

posterior_array<-as.array(fit2)

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

# plot the distributions for each parameter
mcmc_dens(posterior_array, c("sigma","beta", "beta2"),
          facet_args = list(nrow = 2))

## extract the mean and variance for the variables
params<- as.matrix(fit2, pars = c("sigma","beta","beta2"))
params_sum <- data.frame(summary(fit2, pars = c("sigma","beta","beta2"))$summary)

## Extract the PR residuals
PR_Resid<- as.matrix(fit2, pars = c("PR_Resid"))
PR_Resid_sum <- data.frame(summary(fit2, pars = c("PR_Resid"))$summary)


### Make the second stan model for the treatments

stan_data2 <- list(y = PR_Resid_sum$mean, # PR residual from first model is the y 
                  Treatment  = as.numeric(as.factor(Data_sub$Treatment)), # make the treatment numbers 1-4
                  n_treatment = 4, # number of treatments
                  PlateID = Data_sub$PlateID, # plate ID for the random intercepts
                  n_resid = length(Data_sub$PlateID), # the length of the data to calculate residuals from
                  n_PlateID = length(unique(Data_sub$PlateID)), # number of plates
                  m1 = PR_Resid_sum$sd, # sigma for each y
                  n_SD = length(PR_Resid_sum$sd) # count for each SD replicate
                  
                  
)


stan_data2 <- list(y = PR_Resid_sum$mean, # PR residual from first model is the y 
                   Treatment  = Data_sub$Treatment, # make the treatment numbers 1-4
                   n_treatment = 4, # number of treatments
                   PlateID = Data_sub$PlateID, # plate ID for the random intercepts
                   m1 = PR_Resid_sum$sd # sigma for each y
)

dist_fm <- bf(
  y ~ Treatment + (1|PlateID),
  sigma ~ m1
)

dist_brm <- brm(dist_fm, data = stan_data2)

plot(dist_brm)

write("// Stan model for simple linear regression with residuals

data {
 int < lower = 1 > n_treatment; // number of treatments
 int < lower = 1 > n_PlateID; // number of plates
 int < lower = 1 > n_resid; // number of residual values
 vector[n_resid] y; // Response the PR resid from prior model
 int<lower = 1, upper = n_treatment> Treatment[n_resid]; // Treatment variable
 int<lower = 1, upper = n_PlateID> PlateID [n_resid]; // PlateID for random effect of plate
 vector[n_resid] m1; // SD for each y
 
 
  }

parameters {

 vector[n_treatment] beta; // Slope (regression coefficients for x)
 vector[n_PlateID] alpha_plate; // random intercept for each plate
 

}


model {

// conditional mean
 vector[n_resid] mu;
 vector[n_resid] sigma;
// real<lower = 0> sigma[n_resid];

// informative prior on sigma based on model 1

  
 sigma ~ lognormal(m1,1);

// linear combination
 mu = alpha_plate[n_PlateID] + beta[Treatment];

// likelihood function
 y ~ normal(mu, sigma);
 

 
 }

generated quantities {
   vector[n_resid] y_rep;
   vector[n_resid] sigma_hat;
   

  for (i in 1:n_resid) {
 // generate predictions
    
    real y_hat = alpha_plate[PlateID[i]] + beta[Treatment[i]]; // Posterior predictions
      // generate replication values
    
     sigma_hat[i] = lognormal_rng (m1[i],1);
     
     y_rep[i] = normal_rng(y_hat, sigma_hat[i]);
   } // The posterior predictive distribution


}",

"stan_modelJamie2.stan")


# Run the model
stan_model <- "stan_modelJamie2.stan"

# run the model --  increase the interations for better fit
fit3 <- stan(file = stan_model, data = stan_data2, warmup = 5000, iter = 10000, chains = 3, cores = 2, thin = 2,
             seed = 11)


posterior_array3<-as.array(fit3)

# Posterior predictive checks for pH min model
y_rep <- as.matrix(fit3, pars = "y_rep")
y_rep_sum <- data.frame(summary(fit3, pars = "y_rep")$summary)
ppc_dens_overlay(stan_data2$y, y_rep[1:200, ])


# plot the distributions for each parameter
mcmc_dens(posterior_array3, c("beta[1]","beta[2]","beta[3]", "beta[4]"),
          facet_args = list(nrow = 2))



