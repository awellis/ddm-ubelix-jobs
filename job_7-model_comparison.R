library(tidyverse)
library(brms)
library(tidybayes)


final_ddm_fit_1 <- readRDS("models/final_ddm_fit_1.rds.rds")
final_ddm_fit_2 <- readRDS("models/final_ddm_fit_2.rds.rds")
final_ddm_fit_3 <- readRDS("models/final_ddm_fit_3.rds.rds")
final_ddm_fit_4 <- readRDS("models/final_ddm_fit_4.rds.rds")
final_ddm_fit_5 <- readRDS("models/final_ddm_fit_5.rds")
final_ddm_fit_6 <- readRDS("models/final_ddm_fit_6.rds")

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
