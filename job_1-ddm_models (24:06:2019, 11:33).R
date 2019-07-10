# Install packages
# install.packages("tidyverse")
# install.packages("tidybayes")
# install.packages("brms")
# install.packages("modelr")
# install.packages("rstan")
# install.packages("rtdists")
# install.packages("RWiener")
# install.packages("here")
#
# Load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(modelr)
library(rstan)
library(rtdists)
library(RWiener)

#------------------------------------------------------------
# CONSTANTS
ADAPT_DELTA <- 0.99
ITER <- 2000
WARMUP <- floor(ITER/2)
#------------------------------------------------------------

setwd("/home/ubelix/psy/ellis/projects/LuSchumacher")
data <- read_csv("data/FINAL-Dataset-[21.6].csv")


#------------------------------------------------------------
# Remove ID's that did the task not accordingly
data <- data %>%
  filter(ID!=8 & ID!=12 & ID!=14)

data %<>%
  select(ID, block, condition,
         intensity, intF, sayhi, sayhigh,
         correct, rt, handedness, leftButton, order) %>%
  mutate(ID = as_factor(ID),
         instruction = fct_recode(as_factor(condition),
                                  speed = "speeded",
                                  accuracy = "accuracy"),
         order = as.factor(order),
         leftButton = as.factor(leftButton),
         handedness = as_factor(handedness),
         intensity = 100*intensity,
         intF = as_factor(intF),
         sayhigh = fct_relevel(as_factor(sayhigh),
                               "high")) %>%
  drop_na()

# Intensity as factor
data$intensity <- as.factor(data$intensity)


# Remove outlier according to EWMA plus min 200ms
data <- data %>%
  filter(rt>0.3 & rt<4.216) # New cutoff


#------------------------------------------------------------
# Function for initial values
initfun <- function() {
  list(
    b = rnorm(tmp_dat$K),
    b_bs = runif(tmp_dat$K_bs, 1, 2),
    b_ndt = runif(tmp_dat$K_ndt, 0.1, 0.15),
    b_bias = rnorm(tmp_dat$K_bias, 0.5, 0.1),
    sd_1 = runif(tmp_dat$M_1, 0.5, 1),
    z_1 = matrix(rnorm(tmp_dat$M_1*tmp_dat$N_1, 0, 0.01),
                 tmp_dat$M_1, tmp_dat$N_1),
    L_1 = diag(tmp_dat$M_1)
  )
}

# DDM helper function
predict_ddm <- function(fit) {
  pred <- fit$data %>%
    select(ID, instruction, intensity) %>%
    add_predicted_draws(model = fit,
                        negative_rt = TRUE,
                        n = 500) %>%
    mutate(decision = ifelse(.prediction > 0, 1, 0),
           rt = abs(.prediction))
  pred
}

#------------------------------------------------------------
# MODEL 1: Full model (all predictors)
#------------------------------------------------------------
formula <- bf(rt | dec(sayhi) ~ 0 + instruction:intensity + (0 + instruction:intensity | s | ID),
              bs ~ 0 + instruction + (0 + instruction | s | ID),
              ndt ~ 0 + instruction + (0 + instruction | s | ID),
              bias ~ 0 + instruction + (0 + instruction | s | ID))

# Defining priors
prior <- prior(normal(0, 2), class = b) +
  prior(normal(0.2, 1), class = b, dpar = bs) +
  prior(normal(0, 1), class = b, dpar = bias) +
  prior(normal(0.2, 0.1), class = b, dpar = ndt) + # same as Singmann
  prior(student_t(3, 0, 10), class = sd, group = ID)

tmp_dat <- make_standata(formula,
                         data = data,
                         prior = prior,
                         family = wiener(link = "identity",
                                         link_bs = "log",
                                         link_bias = "logit",
                                         link_ndt = "identity"))

final_ddm_fit_1 <- brm(formula,
                       family = wiener(link = "identity",
                                       link_bs = "log",
                                       link_bias = "logit",
                                       link_ndt = "identity"),
                       control = list(max_treedepth = 15,
                                      adapt_delta = ADAPT_DELTA),
                       warmup = WARMUP,
                       iter = ITER,
                       inits = initfun,
                       prior = prior,
                       cores = parallel::detectCores(),
                       data = data,
                       file = here::here("models/final_ddm_fit_1.rds"))


pred_final_ddm_1 <- predict_ddm(final_ddm_fit_1)
saveRDS(pred_final_ddm_1, file = "models/pred_final_ddm_1.rds")


#-- MODEL 2: Saturated model (no predictor for bias) ----
#------------------------------------------------------------
formula <- bf(rt | dec(sayhi) ~ 0 + instruction:intensity + (0 + instruction:intensity | s | ID),
              bs ~ 0 + instruction + (0 + instruction | s | ID),
              ndt ~ 0 + instruction + (0 + instruction | s | ID),
              bias ~ 1 + (1 | ID))

prior <- prior(normal(0, 2), class = b) +
  prior(normal(0.2, 0.5), class = b, dpar = bs) +
  prior(normal(0.2, 0.1), class = b, dpar = ndt) + # same as Singmann
  prior(student_t(3, 0, 10), class = sd, group = ID)

tmp_dat <- make_standata(formula,
                         data = data,
                         prior = prior,
                         family = wiener(link = "identity",
                                         link_bs = "log",
                                         link_bias = "logit",
                                         link_ndt = "identity"))

final_ddm_fit_2 <- brm(formula,
                       family = wiener(link = "identity",
                                       link_bs = "log",
                                       link_bias = "logit",
                                       link_ndt = "identity"),
                       control = list(max_treedepth = 15,
                                      adapt_delta = ADAPT_DELTA),
                       warmup = WARMUP,
                       iter = ITER,
                       inits = initfun,
                       prior = prior,
                       cores = parallel::detectCores(),
                       data = data,
                       file = here::here("models/final_ddm_fit_2.rds"))


pred_final_ddm_2 <- predict_ddm(final_ddm_fit_2)
saveRDS(pred_final_ddm_2, file = "models/pred_final_ddm_2.rds")


#------------------------------------------------------------
# MODEL 3: Saturated model (no predictor for bias) + (instruction does not predict drift rate)
#------------------------------------------------------------
formula <- bf(rt | dec(sayhi) ~ 0 + intensity + (0 + intensity | s | ID),
              bs ~ 0 + instruction + (0 + instruction | s | ID),
              ndt ~ 0 + instruction + (0 + instruction | s | ID),
              bias ~ 1 + (1 | ID))

prior <- prior(normal(0, 2), class = b) +
  prior(normal(0.2, 1), class = b, dpar = bs) +
  prior(normal(0.2, 0.1), class = b, dpar = ndt) + # same as Singmann
  prior(student_t(3, 0, 10), class = sd, group = ID)

tmp_dat <- make_standata(formula,
                         data = data,
                         prior = prior,
                         family = wiener(link = "identity",
                                         link_bs = "log",
                                         link_bias = "logit",
                                         link_ndt = "identity"))

final_ddm_fit_3 <- brm(formula,
                       family = wiener(link = "identity",
                                       link_bs = "log",
                                       link_bias = "logit",
                                       link_ndt = "identity"),
                       control = list(max_treedepth = 15,
                                      adapt_delta = ADAPT_DELTA),
                       warmup = WARMUP,
                       iter = ITER,
                       inits = initfun,
                       prior = prior,
                       cores = parallel::detectCores(),
                       data = data,
                       file = here::here("models/final_ddm_fit_3.rds"))


pred_final_ddm_3 <- predict_ddm(final_ddm_fit_3)
saveRDS(pred_final_ddm_3, file = "models/pred_final_ddm_3.rds")


#------------------------------------------------------------
# MODEL 4: Saturated model (no predictor for bias) + (instruction does not predict drift rate) + (no predictor for boundary seperation)
#------------------------------------------------------------
formula <- bf(rt | dec(sayhi) ~ 0 + intensity + (0 + intensity | s | ID),
              bs ~ 1 + (1 | ID),
              ndt ~ 0 + instruction + (0 + instruction | s | ID),
              bias ~ 1 + (1 | ID))


prior <- prior(normal(0, 2), class = b) +
  prior(normal(0.2, 0.1), class = b, dpar = ndt) + # same as Singmann
  prior(student_t(3, 0, 10), class = sd, group = ID)


tmp_dat <- make_standata(formula,
                         data = data,
                         prior = prior,
                         family = wiener(link = "identity",
                                         link_bs = "log",
                                         link_bias = "logit",
                                         link_ndt = "identity"))

final_ddm_fit_4 <- brm(formula,
                       family = wiener(link = "identity",
                                       link_bs = "log",
                                       link_bias = "logit",
                                       link_ndt = "identity"),
                       control = list(max_treedepth = 15,
                                      adapt_delta = ADAPT_DELTA),
                       warmup = WARMUP,
                       iter = ITER,
                       inits = initfun,
                       prior = prior,
                       cores = parallel::detectCores(),
                       data = data,
                       file = here::here("models/final_ddm_fit_4.rds"))


pred_final_ddm_4 <- predict_ddm(final_ddm_fit_4)
saveRDS(pred_final_ddm_4, file = "models/pred_final_ddm_4.rds")


#------------------------------------------------------------
# MODEL 5: Saturated model (no predictor for bias) + (instruction does not predict drift rate) + (no predictor for non decision time)
#------------------------------------------------------------
formula <- bf(rt | dec(sayhi) ~ 0 + intensity + (0 + intensity | s | ID),
              bs ~ 0 + instruction + (0 + instruction | s | ID),
              ndt ~ 1 + (1 | ID),
              bias ~ 1 + (1 | ID))

prior <- prior(normal(0, 2), class = b) +
  prior(normal(0.2, 1), class = b, dpar = bs) +
  prior(student_t(3, 0, 10), class = sd, group = ID)


tmp_dat <- make_standata(formula,
                         data = data,
                         prior = prior,
                         family = wiener(link = "identity",
                                         link_bs = "log",
                                         link_bias = "logit",
                                         link_ndt = "identity"))

final_ddm_fit_5 <- brm(formula,
                       family = wiener(link = "identity",
                                       link_bs = "log",
                                       link_bias = "logit",
                                       link_ndt = "identity"),
                       control = list(max_treedepth = 15,
                                      adapt_delta = ADAPT_DELTA),
                       warmup = WARMUP,
                       iter = ITER,
                       inits = initfun,
                       prior = prior,
                       cores = parallel::detectCores(),
                       data = data,
                       file = here::here("models/final_ddm_fit_5.rds"))


pred_final_ddm_5 <- predict_ddm(final_ddm_fit_5)
saveRDS(pred_final_ddm_5, file = "models/pred_final_ddm_5.rds")

#------------------------------------------------------------
# MODEL 6: Saturated model (instruction does not predict anything)
#------------------------------------------------------------
formula <- bf(rt | dec(sayhi) ~ 0 + intensity + (0 + intensity | s | ID),
              bs ~ 1 + (1 | ID),
              ndt ~ 1 + (1 | ID),
              bias ~ 1 + (1 | ID))

prior <- prior(normal(0, 2), class = b) +
  prior(student_t(3, 0, 10), class = sd, group = ID)

tmp_dat <- make_standata(formula,
                         data = data,
                         prior = prior,
                         family = wiener(link = "identity",
                                         link_bs = "log",
                                         link_bias = "logit",
                                         link_ndt = "identity"))

final_ddm_fit_6 <- brm(formula,
                       family = wiener(link = "identity",
                                       link_bs = "log",
                                       link_bias = "logit",
                                       link_ndt = "identity"),
                       control = list(max_treedepth = 15,
                                      adapt_delta = ADAPT_DELTA),
                       warmup = WARMUP,
                       iter = ITER,
                       inits = initfun,
                       prior = prior,
                       cores = parallel::detectCores(),
                       data = data,
                       file = here::here("models/final_ddm_fit_6.rds"))


pred_final_ddm_6 <- predict_ddm(final_ddm_fit_6)
saveRDS(pred_final_ddm_6, file = "models/pred_final_ddm_6.rds")


#------------------------------------------------------------
# Model comparison
#------------------------------------------------------------
model_comp <- brms::loo(final_ddm_fit_1,
                        final_ddm_fit_2,
                        final_ddm_fit_3,
                        final_ddm_fit_4,
                        final_ddm_fit_5,
                        final_ddm_fit_6)

saveRDS(model_comp, file = "models/model_comp.rds")
