# some preliminary analysis to do

# install

#install.packages("pacman")
#pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)

#devtools::install_github("daniel1noble/orchaRd", force = TRUE)

library(tidyverse)
library(here)
library(lme4)
library(orchaRd)
#library(gptstudio)
library(metafor)
library(patchwork)
library(alluvial)
library(ggalluvial)
library(easyalluvial)
library(ape)
library(clubSandwich)
library(emmeans)
library(MuMIn)
# making metafor talk to MuMIn
eval(metafor:::.MuMIn)
# install.packages("pak")
#pak::pak("MichelNivard/gptstudio")

#dat_full <- read.csv(here("data/dat_04_04_2023.csv"))
dat_full <- read.csv(here("data/dat_28_06_2023.csv"))

# add phylogenetic tree - only topologies
# TODO? - we could get better tree from birdtree.org
# we can do 50 different trees as in 
# https://academic.oup.com/sysbio/article/68/4/632/5267840
tree_top <- read.tree(here("R/birds_MA.tre"))

# tree with branch lengths
tree <- compute.brlen(tree_top)
plot(tree)
# turning it into a correlation matrix
cor_tree <- vcv(tree,corr=T)

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
dat_full$Phylo <- gsub(" ", "_", dat_full$FocalSpL)

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

# for add-on, we only need V and AV
dat_short_add <- dat %>% filter(Treat_mod == "AV" | Treat_mod == "V")


# some data exploration
# dat %>% group_by(Treat_mod) %>% summarise(n = n())
# using a table to see overalps between Treat_mod and Type
tab <- table(dat_short$Treat_mod, dat_short$Type)

# visualise this table using alluvial plot?
# https://cran.r-project.org/web/packages/alluvial/vignettes/alluvial.html

# Treatment vs Type
dat_short %>% group_by(Treat_mod, Type) %>%
  summarise(n = n()) -> tab1
alluvial(tab1[,1:2], freq = tab1$n)

# using ggaruvial
ggplot(tab1,
       aes(y = n,
           axis1 = Treat_mod,
           axis2 = Type)) +
  geom_alluvium(aes(fill = Treat_mod)) +
  geom_stratum(alpha = 0.5) +
  geom_text(stat = "stratum", size = 6, aes(label = after_stat(stratum))) +
  theme(legend.position = "none") +
  theme(legend.position = "none",
        axis.text.x = element_blank()) + # remove x-axis labels
  ylab("Frequency") + 
  xlab("Treatment modality and trait type")

# Type vs Sex

dat %>% group_by(Sex, Type) %>%
  summarise(n = n()) -> tab2
alluvial(tab2[,1:2], freq = tab2$n)

# using ggaruvial
ggplot(tab2,
       aes(y = n,
           axis1 = Type,
           axis2 = Sex)) +
  geom_alluvium(aes(fill = Type)) +
  geom_stratum(alpha = 0.5) +
  geom_text(stat = "stratum", size = 6, aes(label = after_stat(stratum))) +
  theme(legend.position = "none") +
  theme(legend.position = "none",
        axis.text.x = element_blank()) + # remove x-axis labels
  ylab("Frequency") + 
  xlab("Trait type and sex")


# TOOD - use easyalluvial to make one fig
# other ones - easyalluvial
#https://www.r-bloggers.com/2018/10/data-exploration-with-alluvial-plots-an-introduction-to-easyalluvial/
# using easyalluvial and alluvial_wide
#factor_cols <- dat %>% select_if(is.factor) %>% names()

#alluvial_wide(dat_short, , max_variables = 5
#                , fill_by = 'first_variable' ) %>%
#  add_marginal_histograms(dat_short)

# exploratory analysis
# check each columns for missing values and other stuff
# dat %>% map_df(~sum(is.na(.)))

dat <- dat %>%
  mutate_if(is.character, as.factor)

summary(dat)

#####################
# main meta-analysis
######################

# VCV matrix

VCV <- vcalc(vi = dat$Vd,
             cluster = dat$SubjectID,
             rho = 0.5)

mod0 <- rma.mv(yi = SMD,
       V = VCV, 
       random = list(~1 | Phylo,
                     ~1 | FocalSpL,
                     ~1 | RecNo,
                     ~1 | SubjectID, # incoprated as VCV
                     ~1 | Obs_ID),
       R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat)

summary(mod0)

# TODO - think about whether we add this or not
robust(mod0, cluster = dat$SubjectID)

round(i2_ml(mod0), 2)

# reducted model

mod0r <- rma.mv(yi = SMD,
       V = VCV, 
       random = list(#~1 | Phylo,
                     ~1 | FocalSpL,
                     ~1 | RecNo,
                     #~1 | SubjectID, # incoprated as VCV
                     ~1 | Obs_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat)

summary(mod0r)

round(i2_ml(mod0r), 2)

orchard_plot(mod0,
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)



#################
# meta-regression
#################

## Treatment - A, V, AV etc 
#dat$Phylo <- as.character(dat$Phylo)
#match(dat$Phylo, colnames(cor_tree))
#match(colnames(cor_tree), dat$Phylo)

mod1 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1 | FocalSpL,
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mod = ~ Treat_mod - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)

summary(mod1)

round(r2_ml(mod1)*100, 2)

orchard_plot(mod1, 
             mod = "Treat_mod",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)

###############################################
###############################################
# focusing on A, V, and AV

VCVs <- vcalc(vi = dat_short$Vd,
             cluster = dat_short$SubjectID,
             rho = 0.5)


mod1a <- rma.mv(yi = SMD,
                   V = VCVs,
                   random = list(~1|FocalSpL,
                                 ~1 | RecNo,
                                 ~1 | Obs_ID),
                   mod = ~ Treat_mod, 
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat_short)

# to look at A vs AV and A vs V
# TODO - make these hetero model too
summary(mod1a)

mod1b <- rma.mv(yi = SMD, 
                V = VCVs, 
                random = list(~1|FocalSpL , 
                              ~1 | RecNo, 
                              ~1 | Obs_ID), 
                mod = ~ relevel(factor(Treat_mod), ref = "AV"), 
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dat_short)

# to look at AV vs V and AV vs A (already done)
summary(mod1b)


orchard_plot(mod1a, 
             mod = "Treat_mod",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)

# modeling heteroscedasticity
mod1c <- rma.mv(yi = SMD, 
                V = VCVs, 
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


# the effect of additions
# this is a part of sensitivity analysis

VCVs2 <- vcalc(vi = dat_short_add$Vd,
             cluster = dat_short_add$SubjectID,
             rho = 0.5)

mod5 <- rma.mv(yi = SMD, 
                V = VCVs2, 
                random = list(~1|FocalSpL , 
                              ~1 | RecNo, 
                              ~ Treat_mod | Obs_ID), 
                mod = ~ Treat_mod*Add_on - 1,
                test = "t",
                struct = "DIAG",
                method = "REML", 
                sparse = TRUE,
                data = dat_short_add)

summary(mod5)

################################################
################################################

# Type of responses
mod2 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1 | FocalSpL,
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mod = ~ Type - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)

summary(mod2)

mod2c <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL,
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mod = ~ Type,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)

summary(mod2c)

mod2d <- rma.mv(yi = SMD,
               V = VCV,
               random = list(~1|FocalSpL,
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mod = ~ relevel(Type, ref = "LifeHistory"),
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)

summary(mod2d)

orchard_plot(mod2,
             mod = "Type",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)

round(r2_ml(mod2)*100, 2)

# heteroscadasticity model
mod2b <- rma.mv(yi = SMD, 
               V = VCV, 
               mod = ~ Type - 1, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo, 
                             ~Type | Obs_ID), 
               struct = "DIAG",
               test = "t",
               method = "REML", 
               sparse = TRUE,
               data = dat)

summary(mod2b)


# make other hetero
mod2e <- rma.mv(yi = SMD,
               V = VCV,
               mod = ~ Type, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo, 
                             ~Type | Obs_ID), 
               struct = "DIAG",
               test = "t",
               method = "REML", 
               sparse = TRUE,
               data = dat)

summary(mod2e)

mod2f <- rma.mv(yi = SMD,
               V = VCV,
               mod = ~ relevel(Type, ref = "LifeHistory"), 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo, 
                             ~Type | Obs_ID), 
               struct = "DIAG",
               test = "t",
               method = "REML", 
               sparse = TRUE,
               data = dat)

summary(mod2f)

# heteroscadasticity model better than the homoscedasticity model
# note LifeHistory has lowest variation but this may be expected? 
# as it is less flexiable (e.g. the number of eggs?)
anova(mod2, mod2b)

orchard_plot(mod2b, 
             mod = "Type",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)

########################
# Category of responses

mod3 <- rma.mv(yi = SMD, 
               V = VCV, 
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
               V = VCV, 
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
             xlab = "The number of simuli",
             g = TRUE)

# Predactor guild
# quite heterogeneous
# TODO this could be in random effects - think abou thtis a bit later
mod6 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ PredGuild - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)

summary(mod6)

round(r2_ml(mod6)*100, 2)

orchard_plot(mod6, 
             mod = "PredGuild",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# Setting

mod7 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ Setting - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)

summary(mod7)

round(r2_ml(mod7)*100, 2)

orchard_plot(mod7, 
             mod = "Setting",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# Season

mod8 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ Season - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)

mod8b <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ Season,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)

summary(mod8b)

round(r2_ml(mod8)*100, 2)

orchard_plot(mod8,
             mod = "Season",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# Design
mod9 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ Design - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)
summary(mod9)


mod9b <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ Design,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)

summary(mod9b)

round(r2_ml(mod9)*100, 2)

orchard_plot(mod9,
             mod = "Design",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# Response period
mod10 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ ResponsePeriod - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)
summary(mod10)

round(r2_ml(mod10)*100, 2)

orchard_plot(mod10,
             mod = "ResponsePeriod",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# control type
mod11 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ ControlType - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)
summary(mod11)

round(r2_ml(mod11)*100, 2)

orchard_plot(mod11,
             mod = "ControlType",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# sex
# TODO - this could be interesting
# what is in males and females
mod12 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ Sex - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)
summary(mod12)

round(r2_ml(mod12)*100, 2)

orchard_plot(mod12,
             mod = "Sex",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# shoter data for just males and females
# hetero but no sex effect here
dat_sex <- dat %>% filter(Sex != "both")

VCV3 <- vcalc(vi = dat_sex$Vd,
             cluster = dat_sex$SubjectID,
             rho = 0.5)

mod12a <- rma.mv(yi = SMD, 
                 V = VCV3, 
                 mod = ~ Sex, 
                 random = list(~1|FocalSpL , 
                               ~1 | RecNo, 
                               ~1 | Obs_ID), 
                 #struct = "DIAG",
                 test = "t",
                 method = "REML", 
                 sparse = TRUE,
                 data = dat_sex)

mod12b <- rma.mv(yi = SMD, 
                V = VCV3, 
                mod = ~ Sex, 
                random = list(~1|FocalSpL , 
                              ~1 | RecNo, 
                              ~Sex | Obs_ID), 
                struct = "DIAG",
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dat_sex)

summary(mod12b)

anova(mod12a, mod12b)

orchard_plot(mod12b,
             mod = "Sex",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# age
mod13 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ Age - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)
summary(mod13)

round(r2_ml(mod13)*100, 2)

orchard_plot(mod13,
             mod = "Age",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# type of predator

dat$PredTo[dat$PredTo == ""] <- NA
mod14 <- rma.mv(yi = SMD, 
               V = VCV, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ PredTo - 1,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)
summary(mod14)

round(r2_ml(mod14)*100, 2)

orchard_plot(mod14,
             mod = "PredTo",
             group = "RecNo",
             xlab = "Standardised mean differnece (SMD)")

# treatment duration

dat$ln_duration <- log(dat$duration_days)

mod15 <- rma.mv(yi = SMD,
               V = VCV,
               random = list(~1|FocalSpL,
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ ln_duration,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)
summary(mod15)

mod16 <- rma.mv(yi = SMD,
               V = VCV,
               random = list(~1|FocalSpL,
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ ln_duration*Type,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat)
summary(mod16)

round(r2_ml(mod15)*100, 2)

bubble_plot(mod15,
             mod = "ln_duration",
             group = "RecNo",
             xlab = "log(duration in days)",
             g = TRUE) +
    geom_point(data = dat,
    aes(x = ln_duration, y = SMD,
    color = Type,
    fill = Type,
    size = 1/sqrt(Vd)), alpha = 0.6) +
    scale_color_discrete() + #+ # how to put the legend for colour
    guides(color = "legend")

#p + geom_point(aes(colour = Type))
#scale_colour_manual(values = c("red", "blue", "green"))

#######################
# Mulit-variable models
#######################

dat_short$sln_duration <- scale(log(dat_short$duration_days))

mod_full <- rma.mv(yi = SMD, 
               V = VCVs, 
               random = list(~1|FocalSpL , 
                             ~1 | RecNo,
                             ~1 | Obs_ID),
               mods =  ~ #Design +
                         sln_duration*Type +
                         sln_duration*Treat_mod +
                         Sex,
               test = "t",
               method = "REML",
               sparse = TRUE,
               data = dat_short)
summary(mod_full)

round(r2_ml(mod_full)*100, 2)

orchard_plot(mod_full,
             mod = "Type",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)

orchard_plot(mod_full,
             mod = "Treat_mod",
             group = "RecNo", 
             xlab = "Standardised mean differnece (SMD)",
             branch.size = 3)

int_type <- mod_results(mod_full, mod = "sln_duration", group = "RecNo", weights = "prop",
                                   by = "Type")

bubble_plot(int_type, group = "RecNo", mod = "sln_duration", xlab = "ln(duration in days)",
                     legend.pos = "top.left", condition.nrow = 3)

int_trt <- mod_results(mod_full, mod = "sln_duration", group = "RecNo", weights = "prop",
                        by = "Treat_mod")

bubble_plot(int_trt, group = "RecNo", mod = "sln_duration", xlab = "ln(duration in days)",
            legend.pos = "top.left", condition.nrow = 3)

# mulit-model selection
candidates <- dredge(mod_full, trace = 2)

# displays delta AICc <2
candidates_aic2 <- subset(candidates, delta < 5) 
# model averaging
mr_averaged_aic2 <- summary(model.avg(candidates, delta < 5)) 

# relative importance of each predictor for all the models
importance <- sw(candidates)


######################
# Publiction bias test
######################


# funnel plots
funnel(mod0r, 
       yaxis="seinv",
       type = "rstudent")

# Egger

dat$effectN <- (dat$Ncontrol * dat$NTreat) / (dat$Ncontrol + dat$NTreat)
dat$sqeffectN <- sqrt(dat$effectN)

mod0e <- rma.mv(yi = SMD,
                V = VCV,
                mods = ~ sqeffectN,
                random = list(#~1 | Phylo,
                  ~1 | FocalSpL,
                  ~1 | RecNo,
                  #~1 | SubjectID, # incoprated as VCV
                  ~1 | Obs_ID),
                #R = list(Phylo = cor_tree),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dat)

summary(mod0e)

bubble_plot(mod0e,
            mod = "sqeffectN",
            group = "RecNo",
            xlab = "Effective N",
            g = TRUE)
# decline effect
mod0d <- rma.mv(yi = SMD,
                V = VCV,
                mods = ~ Year,
                random = list(#~1 | Phylo,
                  ~1 | FocalSpL,
                  ~1 | RecNo,
                  #~1 | SubjectID, # incoprated as VCV
                  ~1 | Obs_ID),
                #R = list(Phylo = cor_tree),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dat)

summary(mod0d)

bubble_plot(mod0d,
            mod = "Year",
            group = "RecNo",
            xlab = "Publication year",
            g = TRUE)


# full model
dat_short$effectN <- (dat_short$Ncontrol * dat_short$NTreat) / (dat_short$Ncontrol + dat_short$NTreat)
dat_short$sqeffectN <- sqrt(dat_short$effectN)

mod_fulle <- rma.mv(yi = SMD, 
                   V = VCVs, 
                   random = list(~1|FocalSpL , 
                                 ~1 | RecNo,
                                 ~1 | Obs_ID),
                   mods =  ~ sqeffectN +
                     Year +
                     sln_duration*Type +
                     sln_duration*Treat_mod +
                     Sex,
                   test = "t",
                   method = "REML",
                   sparse = TRUE,
                   data = dat_short)
summary(mod_fulle)

dat_fulle <- qdrg(object = mod_fulle, 
                  data = dat_short)
# marginalized overall mean at vi = 0 and year.c = 0; also weights = "prop" or "cells" average things over proportionally. if not specified, all groups (levels) get the same weights
res_fulle1 <- emmeans(dat_fulle, 
                     specs = ~ sqeffectN,
                     df = mod_fulle$ddf, 
                     weights = "prop")

#####################################################
#####################################################
#####################################################
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

