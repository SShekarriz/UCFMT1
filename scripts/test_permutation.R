#!/usr/bin/env Rscript

##### Setup #####
print('Setup')

# Import packages, set working directory, source functions
set.seed(4)
library(tidyverse)
theme_set(theme_bw())
source('./permutation_functions.R')

# Create the p-value storage matrix
pvals = matrix(nrow = 10, ncol = 3)
colnames(pvals) = c('PEngraftedFeatures', 'PEngraftmentEvents', 
                     'PWeightedEvents')
rownames(pvals) = c('FMTvPl_16s','ResvNoRes_16s',
                    'FMTvPl_species','ResvNoRes_species',
                    'FMTvPl_strain','ResvNoRes_strain',
                    'FMTvPl_mags','ResvNoRes_mags',
                    'FMTvPl_genes','ResvNoRes_genes')


##### 16S Setup #####
print('16S Setup')

### Import the data ###

# Import 16s mapfile and create named vectors of patient IDs/groups for permuting
mapfile <- read.csv("../test_data/test_mapfile.csv")
tx_pv = pat_to_vect(mapfile, 'Fig_lab', 'Treatment')
rs_pv = pat_to_vect(filter(mapfile, Treatment == 'FMT'), 'Fig_lab',
                          'Remission')

# Import the 16s ASV data
marker_lvl_16s <- read.csv("../test_data/test_asv_data1.csv")

# Import donor B 16s ASVs and create a vector
strains16s = marker_lvl_16s$Marker

# Set the cutoffs for 16s presence and absence
cutoff_pres_16S = 0.0001
cutoff_abs_16S = 0

# mrkr_t = t(marker_lvl_16s)
# mrkr_t
# mrkr_t = mrkr_t[-1,]
# mrkr_t = matrix(as.numeric(mrkr_t), ncol = 7)
# 
# colnames(mrkr_t) = marker_lvl_16s$Marker
# rownames(mrkr_t) = colnames(marker_lvl_16s)[-1]
# 
# mrkr_t = (as.data.frame(mrkr_t)
#           %>% rownames_to_column('Study_ID'))
# mrkr_t = (mrkr_t
#           %>% left_join(mapfile))
# mrkr_t %>% filter(Treatment == 'FMT')

### Get the engraftment matrix ###

# convert the 16s data to long format
long_mark_16s = mark_to_long(marker_lvl_16s, strains16s)
# Create a matrix of engraftment counts per patient
engr_16s = get_engraft(long_mark_16s, mapfile, cutoff_abs_16S, cutoff_pres_16S,
                       Treatment %in% c('FMT', 'Placebo'))

engr_16s

#### 16S FMT vs Placebo ####
print('16S FMT vs Placebo')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

#### Fig 1 
cts_16s_plt = count_engraft(engr_16s, tx_pv, 'Treatment')
cts_16s_plt

#### 16S Plot Res vs. NoRes ####

#### Fig 1 
cts_16s_plt = count_engraft(engr_16s, rs_pv, 'Remission')
cts_16s_plt

#### Species Setup ####
print('Species Setup')

### Import the data

# Read in the species-level markers
marker_lvl_sp <- read.csv("../test_data/test_species_data1.csv")

# get a vector of donor B species
strains_sp = strains16s

# Set the cutoff for "presence" with species. This is used for both patient
# samples and donor B
cutoff_pres_sp = 0.005

# Set the species absence cutoff
cutoff_abs_sp = 0

### Get the wide matrix features engrafted (1) or not (0) in patients
long_mark_sp = mark_to_long(marker_lvl_sp, strains_sp)
engr_sp = get_engraft(long_mark_sp, mapfile, cutoff_abs_sp, cutoff_pres_sp,
                       Treatment %in% c('FMT', 'Placebo'))
engr_sp


#### Species FMT vs Placebo ####
print('Species FMT vs Placebo')

### Calculate the observed test statistics
print('Calculate the observed test statistics')

# Get the engraftment counts and calculate the test statistics
cts_sp = count_engraft(engr_sp, tx_pv, 'Treatment')
cts_sp

#### Species Res vs NoRes ####
print('Species Res vs NoRes')

### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_sp_rem = count_engraft(engr_sp, rs_pv, txrm = 'Remission')
cts_sp_rem

#### Strains Setup ####
print('Strains Setup')

### Import the data

# Import the strain-level marker data
marker_lvl_st <- read.csv("../test_data/test_strain_data1.csv")

# Set the presence cutoff used in both donor B and patient samples
cutoff_pres_st = 0.0001

strains_st = strains16s

# Set the strain absence cutoff
cutoff_abs_st = 0

### Get the engraftment matrix

long_mark_st = mark_to_long(marker_lvl_st, strains_st)
engr_st = get_engraft(long_mark_st, mapfile, cutoff_abs_st, cutoff_pres_st,
                       Treatment %in% c('FMT', 'Placebo'))
engr_st


#### Strains FMT vs Placebo ####
print('Strains FMT vs Placebo')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_st = count_engraft(engr_st, tx_pv, 'Treatment')
cts_st

#### Strains Res vs NoRes ####
print('Strains Res vs NoRes')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_st_rem = count_engraft(engr_st, rs_pv, txrm = 'Remission')
cts_st_rem

#### MAGs Setup ####
print('MAGs Setup')

### Import the data

marker_lvl_mg = read.csv('../test_data/test_mag_data1.csv')

# Set the presence and absence cutoffs for patients. DonorB doesn't need these
# because all the MAGs come from donor B
cutoff_pres_mg = 90
cutoff_abs_mg = 25


### Get the engraftment matrix

long_mark_mg = mark_to_long(marker_lvl_mg)
engr_mg = get_engraft(long_mark_mg, mapfile, cutoff_abs_mg, cutoff_pres_mg,
                       Treatment %in% c('FMT', 'Placebo'))
engr_mg


#### MAGs FMT vs Placebo ####
print('MAGs FMT vs Placebo')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_mg = count_engraft(engr_mg, tx_pv, 'Treatment')
cts_mg


#### MAGs Res vs NoRes ####
print('MAGs Res vs NoRes')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_mg_rem = count_engraft(engr_mg, rs_pv, txrm = 'Remission')
cts_mg_rem

#### Genes Setup ####
print('Genes Setup')

### Import the data

marker_lvl_ge <- read.csv("../test_data/test_genes_data1.csv")

# Set the presence and absence cutoffs for patients. These aren't used for donor
# B because all the genes came from the donor B assemblies.
cutoff_pres_ge = 90
cutoff_abs_ge = 25

### Get the engraftment matrix

long_mark_ge = mark_to_long(marker_lvl_ge)
engr_ge = get_engraft(long_mark_ge, mapfile, cutoff_abs_ge, cutoff_pres_ge,
                       Treatment %in% c('FMT', 'Placebo'))
engr_ge


#### Genes FMT vs Placebo ####
print('Genes FMT vs Placebo')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_ge = count_engraft(engr_ge, tx_pv, 'Treatment')
cts_ge

#### Genes Res vs NoRes ####
print('Genes Res vs NoRes')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_ge_rem = count_engraft(engr_ge, rs_pv, txrm = 'Remission')
cts_ge_rem

