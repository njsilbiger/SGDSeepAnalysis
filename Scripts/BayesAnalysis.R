### Bayesian SEM model showing the relationship between SGD and ecosystem functioning library(tidybayes)


### Load Libraries ##########
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(brms)
library(posterior)
library(modelr)
library(rstan)

##### Load data ###########


##### Clean Data ##########


#### Run Analysis #####

## Calculate mixing line

# extract just the Varari Data
SeepDataVarari<-Cdata %>% 
  filter(Location == "Varari", Plate_Seep == "Seep") # Pull out just the seep data

# run the code for the mixing line
set.seed(11)


#TAmix<-br(TA~Salinity*Season)

mod_mixingline<-brm(TA~Salinity+(1|Season), data = SeepDataVarari)


#TA.diff ~ TA[sal = 35] - .prediction[Sal = 35]

# make it a latent variable

make_stancode(TA~Salinity+(1|Season), data = SeepDataVarari)


# Plot the prediction lines
SeepDataVarari %>%
  group_by(Season) %>%
  data_grid(Salinity = seq_range(Salinity, n = 51)) %>%
  add_predicted_draws(mod_mixingline) %>%
 # add_epred_draws(mod_mixingline) %>%
  ggplot(aes(x = Salinity, y = TA, color = Season)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
#  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = SeepDataVarari) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")

SeepDataVarari %>%
  group_by(Season) %>%
  data_grid(Salinity = seq_range(Salinity, n = 51)) %>%
  add_predicted_draws(mod_mixingline)%>%
  head(10)


SeepDataVarari %>%
  group_by(Season) %>%
  data_grid(Salinity = seq_range(Salinity, n = 51)) %>%
  add_predicted_draws(mod_mixingline) %>%
  # add_epred_draws(mod_mixingline) %>%
  ggplot(aes(x = Salinity, y = TA, color = Season)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  #  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = SeepDataVarari) +
  # geom_point(data = Cdata %>% 
  #              filter(Location == "Varari", Plate_Seep == "Plate"), aes(shape = Season)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")


## Same but with silicate

# run the code for the mixing line
set.seed(11)
mod_mixingline2<-brm(TA~Silicate_umolL*Season, data = SeepDataVarari)



# Plot the prediction lines
SeepDataVarari %>%
  group_by(Season) %>%
  data_grid(Silicate_umolL = seq_range(Silicate_umolL, n = 51)) %>%
  add_predicted_draws(mod_mixingline2) %>%
  # add_epred_draws(mod_mixingline) %>%
  ggplot(aes(x = Silicate_umolL, y = TA, color = Season)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  #  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = SeepDataVarari) +
  geom_point(data = Cdata %>% 
               filter(Location == "Varari", Plate_Seep == "Plate"), shape = 5) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")
 
 # ylim(2200, 2300)

### stan data
stan_data <- list(N=length(SeepDataVarari$TA), x = SeepDataVarari$Salinity,
                  y = SeepDataVarari$TA, season  = SeepDataVarari$Season, n_season = 2)

write("// Stan model for simple linear regression

data {
 int < lower = 1 > N; // Sample size
 int n_season; // number of seasons
 vector[N] x; // Predictor
 vector[N] y; // Outcome
 int<lower = 0, upper = n_season> season [N]; // season prediction
}

parameters {
 vector[n_season] alpha; // vector of intercepts
 real beta; // Slope (regression coefficients)
 real < lower = 0 > sigma; // Error SD
}

model {

// conditional mean
vector[N] mu;

// linear combination
 mu = alpha[season] + x * beta;

// likelihood function
 y ~ normal(mu, sigma);
}

generated quantities {
vector[N] y_rep;

 for (i in 1:N) {
 // generate predictions
 real y_hat = alpha[season[i]] + x[i] * beta;
 
 // generate replication values
 y_rep[i] = normal_rng(y_hat, sigma);
} // The posterior predictive distribution

}",

"stan_model1.stan")

stan_model1 <- "stan_model1.stan"

fit <- stan(file = stan_model1, data = stan_data, warmup = 500, iter = 1000, chains = 4, cores = 2, thin = 1)

posterior <- extract(fit)
stan_dens(fit)
plot(fit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon")


y_rep <- as.matrix(fit, pars = "y_rep")

ppc_dens_overlay(SeepDataVarari$TA, y_rep[1:200, ])
