# some preliminary analysis to do

# install

#install.packages("pacman")
#pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)

#devtools::install_github("daniel1noble/orchaRd", force = TRUE)

library(tidyverse)
library(here)
library(lme4)
library(orchaRd)
library(gptstudio)
library(metafor)
library(patchwork)
library(alluvial)

# install.packages("pak")
#pak::pak("MichelNivard/gptstudio")

#dat_full <- read.csv(here("data/dat_04_04_2023.csv"))
dat_full <- read.csv(here("data/dat_23_06_2023.csv"))

# function for calculating variance
Vd_func <- function(d, n1, n2, design, r = 0.5){
  # independent design
  if(design == "among"){
    var <- (n1 + n2) / (n1*n2) + d^2 / (2 * (n1 + n2 - 2)) # variance
  } else { # dependent design
    var <- 2*(1-r) / n1 + d^2 / (2*(n1 - 1)) # variance
  }
  var # return variance
}

# getting Hedges' g - get small size bias corrected effect size
dat_full$SMD <- dat_full$d / (1 - 3/(4 * (dat_full$NTreat + dat_full$Ncontrol) - 9))

# flipping d 
dat_full$SMD <- dat_full$d*dat_full$Direction*dat_full$PredictedDirection


# calucating Vd
dat_full$Vd <- with(dat_full, pmap_dbl(list(SMD, NTreat, Ncontrol, Design), Vd_func))


# extra useful function
# function for getting mean and sd from median, quartiles and sample size
# get_mean_sd <- function(median, q1, q3, n){
#   sd <- (q3 - q1) / (2 * (qnorm((0.75 * n - 0.125) / (n + 0.25)))) # sd
#   mean <- (median + q1 + q3)/3 # mean
#   c(mean, sd)
# }


# observation id
dat_full$Obs_ID <- 1:nrow(dat_full)

# filtering very large variance and also very small sample size
dat_int <- dat_full %>% filter(Vd < 10 & Ncontrol > 2 & NTreat > 2)

dim(dat_full)
dim(dat_int)


# sorting out modality stuff
# creat - 1,2,3 modality - also easier classification A, O, V (AOV = L) 

dat_int %>% mutate(Treat_mod = case_when(Treatment == "A" ~ "A",
                                          Treatment == "AV" ~ "AV",
                                          Treatment == "AVG" ~ "AV",
                                          Treatment == "AVM" ~ "AV",
                                          Treatment == "L" ~ "AVO",
                                          Treatment == "O" ~ "O",
                                          Treatment == "OV" ~ "OV",
                                          Treatment == "V" ~ "V",
                                          Treatment == "VG" ~ "V",
                                          Treatment == "VM" ~ "V",
                                          Treatment == "VP" ~ "V"),
                    # into how many
                    Treat_No = case_when(Treatment == "A" ~ 1,
                                         Treatment == "AV" ~ 2,
                                         Treatment == "AVG" ~ 2,
                                         Treatment == "AVM" ~ 2,
                                         Treatment == "L" ~ 3,
                                         Treatment == "O" ~ 1,
                                         Treatment == "OV" ~ 2,
                                         Treatment == "V" ~ 1,
                                         Treatment == "VG" ~ 1,
                                         Treatment == "VM" ~ 1,
                                         Treatment == "VP" ~ 1),
                    # des it have some add-ons
                    Add_on = case_when(Treatment == "A" ~ "No",
                                         Treatment == "AV" ~ "No",
                                         Treatment == "AVG" ~ "Yes",
                                         Treatment == "AVM" ~ "Yes",
                                         Treatment == "L" ~ "No",
                                         Treatment == "O" ~ "No",
                                         Treatment == "OV" ~ "No",
                                         Treatment == "V" ~ "No",
                                         Treatment == "VG" ~ "Yes",
                                         Treatment == "VM" ~ "Yes",
                                         Treatment == "VP" ~ "Yes"),

                      ) -> dat

# creating data just for A, V, and AV 
dat_short <- dat %>% filter(Treat_mod == "A" | Treat_mod == "V" | Treat_mod == "AV")


# some data exploration
# dat %>% group_by(Treat_mod) %>% summarise(n = n())
# using a table to see overalps between Treat_mod and Type
tab <- table(dat_short$Treat_mod, dat_short$Type)

# visualise this table using alluvial plot?
# https://cran.r-project.org/web/packages/alluvial/vignettes/alluvial.html

dat_short %>% group_by(Treat_mod, Type) %>%
  summarise(n = n()) -> tab1
alluvial(tab1[,1:2], freq = tab1$n)

# main meta-analysis

mod0 <- rma.mv(yi = SMD, 
       V = Vd, 
       random = list(~1|FocalSpL ,
                     ~1 | RecNo,
                     ~1 | Obs_ID), 
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat)

summary(mod0)

orchard_plot(mod0,
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)


# meta-regression
## Treatment
mod1 <- rma.mv(yi = SMD, 
               V = Vd, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo, 
                             ~1 | Obs_ID), 
               mod = ~ Treat_mod - 1, 
               test = "t",
               method = "REML", 
               sparse = TRUE,
               data = dat)

summary(mod1)

orchard_plot(mod1, 
             mod = "Treat_mod",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)


mod1a <- rma.mv(yi = SMD, 
                   V = Vd, 
                   random = list(~1|FocalSpL , 
                                 ~1 | RecNo, 
                                 ~1 | Obs_ID), 
                   mod = ~ Treat_mod, 
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat_short)

summary(mod1a)

mod1b <- rma.mv(yi = SMD, 
                V = Vd, 
                random = list(~1|FocalSpL , 
                              ~1 | RecNo, 
                              ~1 | Obs_ID), 
                mod = ~ relevel(factor(Treat_mod), ref = "AV"), 
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dat_short)

summary(mod1b)


orchard_plot(mod1a, 
             mod = "Treat_mod",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)


mod1c <- rma.mv(yi = SMD, 
                V = Vd, 
                random = list(~1|FocalSpL , 
                              ~1 | RecNo, 
                              ~ Treat_mod | Obs_ID), 
                mod = ~ Treat_mod, 
                test = "t",
                struct = "DIAG",
                method = "REML", 
                sparse = TRUE,
                data = dat_short)

summary(mod1c)

anova(mod1a, mod1c)

orchard_plot(mod1c, 
             mod = "Treat_mod",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)


# Type of responses
mod2 <- rma.mv(yi = SMD, 
               V = Vd, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo, 
                             ~1 | Obs_ID), 
               mod = ~ Type - 1, 
               test = "t",
               method = "REML", 
               sparse = TRUE,
               data = dat)

summary(mod2)

orchard_plot(mod2, 
             mod = "Type",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)

# heteroscadasticity model

mod2b <- rma.mv(yi = SMD, 
               V = Vd, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo, 
                             ~Type | Obs_ID), 
               mod = ~ Type - 1, 
               struct = "DIAG",
               test = "t",
               method = "REML", 
               sparse = TRUE,
               data = dat)

summary(mod2b)

orchard_plot(mod2b, 
             mod = "Type",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)

# Category of responses

mod3 <- rma.mv(yi = SMD, 
               V = Vd, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo, 
                             ~1 | Obs_ID), 
               mod = ~ Category - 1, 
               test = "t",
               method = "REML", 
               sparse = TRUE,
               data = dat)

summary(mod3)

orchard_plot(mod3, 
             mod = "Category",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             angle = 45,
             branch.size = 3)

# testing the number of stimuli


mod4 <- rma.mv(yi = SMD, 
               V = Vd, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo, 
                             ~1 | Obs_ID), 
               mod = ~ Treat_No, 
               test = "t",
               method = "REML", 
               sparse = TRUE,
               data = dat)

summary(mod4)

bubble_plot(mod4,
             mod = "Treat_No",
             group = "RecNo",
             xlab = "The number of simuli")


############
# Not run
############
## putting results separtely

dat_ant <- dat[dat$Category == "AntiPredator" , ]
dat_cost <- dat[dat$Category == "ParentalCare" | 
                dat$Category == "Intake"|
                dat$Category == "ConstlyBehaviours", ]

nrow(dat_ant)
nrow(dat_cost)
## Treatment
mod1a <- rma.mv(yi = SMD, 
               V = Vd, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo, 
                             ~1 | Obs_ID), 
               mod = ~ Treat_mod - 1, 
               test = "t",
               method = "REML", 
               sparse = TRUE,
               data = dat_ant)

summary(mod1a)

orchard_plot(mod1a, 
             mod = "Treat_mod",
             group = "RecNo", 
             data = dat_ant, 
             xlab = "SDM: Ant-Predatory Behaviour")

## Treatment
mod1b <- rma.mv(yi = SMD, 
               V = Vd, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mod = ~ Treat_mod - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat_cost)

summary(mod1b)

orchard_plot(mod1b, 
             mod = "Treat_mod",
             group = "RecNo",
             #data = dat_cost,
             xlab = "SMD: Care, Intake and Costly")


