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
# Bayes model
mod1<-brm(bf(y~x+x2), data = data.frame(x = Data_sub_control$minpH, # use pH min for this example as the independent variable
                                    x2 = Data_sub_control$Initial_PR, # Covariate
                                    y = Data_sub_control$Final_PR), iter = 10000, chains = 3,
          save_pars = save_pars(all = TRUE)) 



#posterior predictive chec
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
# set up the model
dist_fm <- bf(
  resid ~ Treatment + (1|PlateID),
  sigma ~ resid_sd # let sigma vary by the SD of the residuals data
)

# run the bayes model
dist_brm <- brm(dist_fm, data = Data_combd)

# posterior predictive checks
pp_check(dist_brm)

pp_check(dist_brm, type = "scatter_avg_grouped", group = "Treatment") + 
    geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)

#plot(dist_brm)

p1<- conditional_effects(dist_brm) 

plot(p1, plot = FALSE)[[1]] + # plot the effects 
  geom_hline(yintercept = 0)

