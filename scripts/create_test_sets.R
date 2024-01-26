# The plan for this is to create fake data sets for each data type that have a
# known number of engrafted features in each of three categories (placebo,
# fmt-res, fmt-nores) and then run the functions and make sure the numbers come
# out right.

library(tidyverse)
source('permutation_functions.R')

#### ASVs ####

# Import 16s mapfile and create named vectors of patient IDs/groups for permuting
mapfile <- read.csv("../data/UCFMT1_16Smap.txt")
tx_pv = pat_to_vect(mapfile, 'Fig_lab', 'Treatment')
rs_pv = pat_to_vect(filter(mapfile, Treatment == 'FMT'), 'Fig_lab', 
                          'Remission')

# Import the 16s ASV data
marker_lvl_16s <- read.csv("../data/16s_lvl_RA.csv")

# Import donor B 16s ASVs and create a vector
strains16s <- read.csv("../data/jcsz_16s_core_DonB.csv")
strains16s %>%
  pull() -> strains16s

# Set the cutoffs for 16s presence and absence
cutoff_pres_16S = 0.0001
cutoff_abs_16S = 0

dim(marker_lvl_16s)
marker_lvl_16s[1:10,1:10]

#### Species ####

#### Strains ####

#### MAGs ####

#### Genes ####