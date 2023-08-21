# sorting out 100 trees for birds

library(ape)
library(tidyverse)
library(here)
library(readr)

dat <- read.csv(here("data/dat_19_07_2023.csv"))
dat$FocalSpL_corrected <- dat$FocalSpL

# tree data
tree_all <- read.tree(here("phy.tre"))
tree_50 <- head(tree_all, 50)
#tree <- tree_all[[1]]
#saveRDS(tree, here("misc", "tree.RDS"))
tree <- readRDS(here("tree", "tree.RDS"))
#write.tree(tree, here("misc", "tree.tre"))


#fix species in teh data set to match the tree
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Troglodytes aedon musculus", "Troglodytes aedon")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Troglodytes aedon aedon", "Troglodytes aedon")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Malothrus mater", "Molothrus ater")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Tetrastes sewesrzowi", "Tetrastes sewerzowi")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Puffinus pacificu", "Puffinus pacificus")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Aphelocoma coerulescents", "Aphelocoma coerulescens")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Thamnomanes schisogynus", "Thamnomanes schistogynus")

dat$FocalSpL_corrected <- gsub(" ","_", dat$FocalSpL_corrected) #get rid of the underscores

# how many unique species
length(unique(dat$FocalSpL_corrected))

# check overlap and differences with taxa list
intersect(unique(dat$FocalSpL_corrected), tree$tip.label) # 81 overlap
setdiff(unique(dat$FocalSpL_corrected), tree$tip.label) # 7 in data, not in tree - mostly synonyms

# note that Parus_minor should be Parus_major but we use the tree for Parus_moticolus
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Chloris_chloris", "Carduelis_chloris")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Cyanistes_caeruleus", "Parus_caeruleus")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Setophaga_petechia", "Dendroica_petechia")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Haemorhous_mexicanus", "Carpodacus_mexicanus")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Poecile_gambeli", "Parus_gambeli")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Poecile_atricapillus", "Parus_atricapillus")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Parus_minor", "Parus_monticolus")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Tetrastes_sewerzowi", "Bonasa_sewerzowi")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Setophaga_citrina", "Wilsonia_citrina")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Poecile_carolinensis", "Parus_carolinensis")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Oreothlypis_celata", "Vermivora_celata")
dat$FocalSpL_corrected <- replace(dat$FocalSpL_corrected, dat$FocalSpL_corrected == "Poecile_palustris", "Parus_palustris")
# checking again
intersect(unique(dat$FocalSpL_corrected), tree$tip.label) # 81 overlap
setdiff(unique(dat$FocalSpL_corrected), tree$tip.label) # 7 in data, not in tree - mostly synonyms

length(unique(dat$FocalSpL_corrected)) # 85 species

take_out <- setdiff(tree$tip.label, unique(dat$FocalSpL_corrected))

tree_50 <- map(tree_50, ~drop.tip(.x, take_out))

saveRDS(tree_50, here("tree", "tree_50.RDS"))

write.csv(dat, file = "data/dat_19_07_2023_spp.csv", row.names = F)
