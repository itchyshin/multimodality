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

# install.packages("pak")
#pak::pak("MichelNivard/gptstudio")

dat_full <- read.csv(here("data/dat_04_04_2023.csv"))

# function for calculating variance
Vd <- function(d, n1, n2, design, r = 0.5){
  # independent design
  if(design == "among"){
    var <- (n1 + n2) / (n1*n2) + d^2 / (2*(n1 + n2 -2))
  } else { # dependent design
    var <- 2*(1-r) / n1 + d^2 / (2*(n1 - 1))
  }
  var
}

# this does not work 
# this is because "Error in if (design == "among") { : the condition has length > 1"
#dat_full$Vd <- with(dat, Vd(d, NTreat, Ncontrol, design = Design))


# TODO - some data checking needs to be a bit later

# this works
dat_full$Vd <- with(dat_full, mapply(Vd_func, d, 
                                n1 = NTreat, 
                                n2 = Ncontrol, 
                                design = Design))

# observation id
dat_full$Obs_ID <- 1:nrow(dat_full)

# flipping d 

dat_full$SMD <- dat_full$d*dat_full$Direction

# filtering very large variance and also very small sample size
dat_full %>% filter(Vd < 10 & Ncontrol > 3 & NTreat > 3) -> dat_int

dim(dat_full)
dim(dat_int)

# TODO 
# sorting out modality stuff
# creat - 1,2,3 modality - also easier classification A, O, V (AOV = L)
# summary(as.factor(dat_full$Treatment))
#    A  AV AVG AVM   L   O  OV   V  VG  VM  VP 
# 344  78   6  16   9  24  21 161  22  16   4 
# expectation is

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
             data = dat, 
             xlab = "Standardised mean differnece (SMD)")


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
             data = dat, 
             xlab = "Standardised mean differnece (SMD)")

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
             data = dat, 
             xlab = "Standardised mean differnece (SMD)")

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
             data = dat, 
             xlab = "Standardised mean differnece (SMD)",
             angle = 45)

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
             data = dat_cost, 
             xlab = "SMD: Care, Intake and Costly")


