# load packages
library(rbmi)



######################## SIMULATE DATA
set.seed(27082023)

# simulate data using simulate_data from rbmi (?simulate_data for info).
# same simulation parameters as ?get_example_data 
# except that we simulate also ICE2 (probability of 2% after each visit)

# input parameters to simulate_data
n <- 100 # sample size per arm
time <- c(0, 4, 8, 12) # visits

# assume 10 point yearly increase in the control arm
# assume treatment effect starts after 2 months with 50% relative reduction
muC <- c(50, 53.33333, 56.66667, 60) # outcome mean of control arm
muT <- c(50, 53.33333, 55, 56.66667) # outcome mean of intervention arm

# covariance matrix (same for control and treatment arm)
sd_error <- 2.5
covRE <- rbind(c(25, 6.25), c(6.25, 25))
Sigma <- cbind(1, time/12) %*% covRE %*% rbind(1, time/12) + 
  diag(sd_error^2, nrow = length(time))

# parameters to simulate treatment discontinuation due to SDCR reasons
probDisc_C <- 0.02
probDisc_T <- 0.03
or_outcome <- 1.1
prob_dropout <- 0.5

# parameter to simulate treatment discontinuation due to NSDCR reasons
prob_ice2 <- 0.02

# set simulation parameters for the control arm using set_simul_pars
parsC <- set_simul_pars(
  mu = muC,
  sigma = Sigma,
  n = n,
  prob_ice1 = probDisc_C, 
  or_outcome_ice1 = or_outcome, 
  prob_post_ice1_dropout = prob_dropout,
  prob_ice2 = prob_ice2
)

# set simulation parameters of the treatment arm
parsT <- parsC
parsT$mu <- muT
parsT$prob_ice1 <- probDisc_T
post_ice_traj <- "CIR"

# simulate data
dat <- simulate_data(
  pars_c = parsC,
  pars_t = parsT,
  post_ice1_traj = post_ice_traj
)

# let's summarise the ICEs
dat %>% 
  group_by(visit) %>%
  summarise(
    frequency_ice1 = sum(ind_ice1 == 1)/n(),
    frequency_dropout_ice1 = sum(dropout_ice1 == 1)/n(),
    frequency_ice2 = sum(ind_ice2 == 1)/n()
  )

######################## PREPROCESSING

# pre-process data: remove baseline from outcome variable (since we 
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
# such data!!
all(is.na(dat$chg[dat$ind_ice2 == 1]))
# ok, in our case we don't have data collected after ICE2

# let's have a look at the first of the 4 core functions of rbmi: "draws()"
?draws

########### data_ice

# specify CIR imputation for data after ICE1
dat_ice1 <- dat %>% 
  filter(ind_ice1 == 1) %>% 
  group_by(id) %>%
  slice(1) %>%
  mutate(strategy = "CIR") %>%
  select(id, visit, strategy)

# specify MAR imputation for data after ICE2
# not necessary since rbmi imputes under MAR if nothing else
# is specified
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

# use the function set_vars 
vars <- set_vars(
  subjid = "id",
  visit = "visit",
  outcome = "chg",
  group = "group",
  covariates = c("outcome_bl*visit", "group*visit"),
  strategy = "strategy"
)


########### method

# two methods: Bayesian MI and condmean + jackknife

method_bayesian <- method_bayes(
  burn_in = 200,
  burn_between = 50,
  n_samples = 100
)

method_cm_jackknife <- method_condmean(
  type = "jackknife"
)


######################## RUN ESTIMATORS

# Byesian MI

draws_obj <- draws(
  data = dat,
  data_ice = dat_ice,
  vars = vars,
  method = method_bayesian
)

references <- c(
  "Control" = "Control",
  "Intervention" = "Control"
)

impute_obj <- impute(
  draws = draws_obj,
  references = references
)

# analysis model: ANCOVA
# set_vars for ANCOVA
vars_ancova <- vars
vars_ancova$covariates <- c("outcome_bl", "group")

an_obj <- analyse(
  impute_obj,
  vars = vars_ancova,
  visits = "6"
)

pool_obj_bayes <- pool(
  an_obj
)
pool_obj_bayes


# conditional mean imputation + jackknife
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

# analysis model: ANCOVA
# set_vars for ANCOVA
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




