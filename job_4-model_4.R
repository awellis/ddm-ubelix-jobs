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
  mutate(ID = as.factor(ID),
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
# initfun <- function() {
#   list(
#     b = rnorm(tmp_dat$K),
#     b_bs = runif(tmp_dat$K_bs, 0.9, 1.5),
#     b_ndt = runif(tmp_dat$K_ndt, -3.1, -2.9),
#     b_bias = rnorm(tmp_dat$K_bias, 0.0, 0.1),
#     sd_1 = runif(tmp_dat$M_1, 0.5, 1),
#     z_1 = matrix(rnorm(tmp_dat$M_1*tmp_dat$N_1, 0, 0.01),
#                  tmp_dat$M_1, tmp_dat$N_1),
#     L_1 = diag(tmp_dat$M_1)
#   )
# }

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
# MODEL 4: Saturated model (no predictor for bias) + (instruction does not predict drift rate) + (no predictor for boundary seperation)
#------------------------------------------------------------
formula <- bf(rt | dec(sayhi) ~ 0 + intensity + (0 + intensity | s | ID),
              bs ~ 1 + (1 | s | ID),
              ndt ~ 1 + instruction + (1 + instruction | s | ID),
              bias ~ 1 + (1 | s | ID))

prior <- prior(normal(0, 2), class = b) +
         prior(normal(0, 1), class = Intercept, dpar = bs) +
         prior(normal(0, 5), class = Intercept, dpar = ndt) +
         prior(normal(0, 1), class = b, dpar = ndt) +
         prior(normal(0, 1), class = Intercept, dpar = bias) +
         prior(student_t(3, 0, 1), class = sd, group = ID) +
         prior(lkj(2), class = cor, group = ID)

# tmp_dat <- make_standata(formula,
#                          data = data,
#                          prior = prior,
#                          family = wiener(link = "identity",
#                                          link_bs = "log",
#                                          link_bias = "logit",
#                                          link_ndt = "log"))


initfun <- function() {
  list(temp_ndt_Intercept = -3)
}

final_ddm_fit_4 <- brm(formula,
                       family = wiener(link = "identity",
                                       link_bs = "log",
                                       link_bias = "logit",
                                       link_ndt = "log"),
                       control = list(max_treedepth = 15,
                                      adapt_delta = ADAPT_DELTA),
                       warmup = WARMUP,
                       iter = ITER,
                       inits = initfun,
                       init_r = 0.05,
                       prior = prior,
                       cores = parallel::detectCores(),
                       data = data) %>%
                    add_loo()


pred_final_ddm_4 <- predict_ddm(final_ddm_fit_4)
saveRDS(final_ddm_fit_4, file = "models/fit_4_ddm_final.rds")
saveRDS(pred_final_ddm_4, file = "models/pred_4_ddm_final.rds")
