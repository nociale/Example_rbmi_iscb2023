# Clear environment
rm(list = ls())

# Load packages
library(rbmi)
library(dplyr)


# Please check out the rbmi vignettes
vignette("quickstart",package="rbmi")  # Simple example of the rbmi functions
vignette("advanced",package="rbmi")    # Advanced use of the rbmi functions
vignette("stat_specs", package="rbmi") # Statistical Specifications


######################## SIMULATE DATA
set.seed(27082023)

# Simulate data using simulate_data from rbmi (`?simulate_data` for info).
# Similar simulation parameters to `?get_example_data`

# Input parameters to simulate_data
n <- 100 # sample size per arm
time <- c(0, 4, 8, 12) # visits

# Assume 10 point yearly increase in the control arm.
# Assume treatment effect starts after 4 months with 50% relative reduction
muC <- c(50, 53.33333, 56.66667, 60) # outcome mean of control arm
muT <- c(50, 53.33333, 55, 56.66667) # outcome mean of intervention arm

# Covariance matrix (same for control and treatment arm)
sd_error <- 2.5
covRE <- rbind(c(25, 6.25), c(6.25, 25))
Sigma <- cbind(1, time/12) %*% covRE %*% rbind(1, time/12) + 
  diag(sd_error^2, nrow = length(time))

# Parameters to simulate treatment discontinuation due to SDCR reasons
probDisc_C <- 0.04
probDisc_T <- 0.06
or_outcome <- 1.1
prob_dropout <- 0.5

# Parameter to simulate treatment discontinuation due to NSDCR reasons
prob_ice2 <- 0.03

# Set simulation parameters of the control arm using `set_simul_pars`
parsC <- set_simul_pars(
  mu = muC,
  sigma = Sigma,
  n = n,
  prob_ice1 = probDisc_C, 
  or_outcome_ice1 = or_outcome, 
  prob_post_ice1_dropout = prob_dropout,
  prob_ice2 = prob_ice2
)

# Set simulation parameters of the treatment arm
parsT <- parsC
parsT$mu <- muT
parsT$prob_ice1 <- probDisc_T
post_ice_traj <- "CIR"

# Simulate data
dat <- simulate_data(
  pars_c = parsC,
  pars_t = parsT,
  post_ice1_traj = post_ice_traj
)

# Let's summarise the ICEs
dat %>% 
  group_by(visit) %>%
  summarise(
    frequency_ice1 = sum(ind_ice1 == 1)/n(),
    frequency_dropout_ice1 = sum(dropout_ice1 == 1)/n(),
    frequency_ice2 = sum(ind_ice2 == 1)/n()
  )





######################## PREPROCESSING

# Let's have a look at the first of the 4 core functions of rbmi: `draws()`
?draws


########### data

# Pre-process data: remove baseline from outcome variable (since we 
# model the change from baseline)
dat <- dat %>% 
  filter(visit != 0) %>%
  mutate(
    chg = outcome - outcome_bl,
    visit = factor(visit, levels = unique(visit))
  )

# Be careful: If you want to implement
# a hypothetical strategy which prescribes to remove
# data collected after the ICE, set to NA
# such data!
all(is.na(dat$chg[dat$ind_ice2 == 1]))
# ok, in our case we don't have data collected after ICE2


########### data_ice

# Specify CIR imputation for data after ICE1
dat_ice1 <- dat %>% 
  filter(ind_ice1 == 1) %>% 
  group_by(id) %>%
  slice(1) %>%
  mutate(strategy = "CIR") %>%
  select(id, visit, strategy)

# Specify MAR imputation for data after ICE2.
# This is not necessary since rbmi automatically imputes under MAR 
# if nothing else is specified
dat_ice2 <- dat %>% 
  filter(ind_ice2 == 1) %>% 
  group_by(id) %>%
  slice(1) %>%
  mutate(strategy = "MAR") %>%
  select(id, visit, strategy)

dat_ice <- rbind(
  dat_ice1,
  dat_ice2
)


########### vars

# Use the function `set_vars`
vars <- set_vars(
  subjid = "id",
  visit = "visit",
  outcome = "chg",
  group = "group",
  covariates = c("outcome_bl*visit", "group*visit"),
  strategy = "strategy"
)


########### method

# Two methods: Bayesian MI and condmean + jackknife

method_bayesian <- method_bayes(
  burn_in = 200,
  burn_between = 50,
  n_samples = 150
)

method_cm_jackknife <- method_condmean(
  type = "jackknife"
)




######################## RUN ESTIMATORS

# Bayesian MI

draws_obj <- draws(
  data = dat,
  data_ice = dat_ice,
  vars = vars,
  method = method_bayesian
)

# Set references (for reference-based imputation)
references <- c(
  "Control" = "Control",
  "Intervention" = "Control"
)

impute_obj <- impute(
  draws = draws_obj,
  references = references
)

# Analysis model: ANCOVA.
# The rbmi function `ancova` can be used as analysis function
# Additional arguments passed to `analyse` are here `vars` and 
# `visits` (used to specify at which visits we want to perform
# the analysis. For `vars` we use again `set_vars`.
vars_ancova <- vars
vars_ancova$covariates <- c("outcome_bl", "group")

an_obj <- analyse(
  impute_obj,
  fun = ancova,
  vars = vars_ancova,
  visits = "3" # last visit
)

# Pool analysis results using Rubin's rules
pool_obj_bayes <- pool(
  an_obj
)
pool_obj_bayes


# Conditional mean imputation + jackknife
# (Same code as before except for the `method` argument)

draws_obj <- draws(
  data = dat,
  data_ice = dat_ice,
  vars = vars,
  method = method_cm_jackknife
)

references <- c(
  "Control" = "Control",
  "Intervention" = "Control"
)

impute_obj <- impute(
  draws = draws_obj,
  references = references
)

vars_ancova <- vars
vars_ancova$covariates <- c("outcome_bl", "group")

an_obj <- analyse(
  impute_obj,
  vars = vars_ancova,
  visits = "3"
)

pool_obj_cm <- pool(
  an_obj
)
pool_obj_cm
pool_obj_bayes

# note the difference in SE estimation due to the difference between
# information-anchored inference (targeted by Rubin's rules) 
# and frequentist inference (targeted by conditional mean imputation +
# jackknife) under reference-based assumption





#################### DELTA-ADJUSTMENT 

# Let's assume we are interested in whether the results would be 
# significant if the imputed values in the treatment arm are 2 points
# worse than what is imputed in the primary analysis.

# The `delta` argument of `analyse()` allows users to modify the 
# outcome variable prior to the analysis. 
# To do this, the user needs to provide a data.frame 
# which contains columns for the subject and visit 
# (to identify the observation to be adjusted) plus an additional 
# column called delta which specifies the value which will be added 
# to the outcomes prior to the analysis.
# The `delta_template()` function supports the user in creating
# this data.frame
dat_delta <- delta_template(imputations = impute_obj) %>%
  mutate(delta = (is_missing & group == "Intervention") * 2)

ana_delta <- analyse(
  impute_obj,
  delta = dat_delta,
  vars = vars_ancova,
  visits = "3"
)
pool(ana_delta)
pool_obj_cm

# This can be used to implement a "tipping point analysis"