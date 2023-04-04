# some preliminary analysis to do

# install

#install.packages("pacman")
#pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)

#devtools::install_github("daniel1noble/orchaRd", force = TRUE)

library(tidyverse)
library(here)
library(lme4)
library(orchaRd)

dat <- read.csv(here("data/dat_04_04_2023.csv"))

# function for calculating variance
Vd < function(d, n1, n2, design = c("within","between"), r = 0.5){}



# filtering some crazy values
dat %>% 





