# Data Processing ####

# This script is to do all data wrangling, organization, etc. to create the data
# structures that will be used by run_permutation.R, make_paperFigs.R, and any
# other analysis or visualizing in the UCFMT1 metagenomics project.

# Must be run from the top-level project directory. Can be run interactively or 
# using Rscript

## Setup ####

library(tidyverse)
source('./scripts/permutation_functions.R')
source('./scripts/paperfig_functions.R')

## 16S ####

### Setup ####

# Import 16s mapfile and create named vectors of patient IDs/groups for 
# permuting
mapfile_16s <- read.csv("./data/UCFMT1_16Smap.txt")
tx_pv_16s = pat_to_vect(mapfile_16s, 'Fig_lab', 'Treatment')
rs_pv_16s = pat_to_vect(filter(mapfile_16s, Treatment == 'FMT'), 'Fig_lab', 
                        'Remission')

### Save Setup ####

setup_16s = c('mapfile_16s', 'tx_pv_16s', 'rs_pv_16s')
if (!dir.exists('./processed_data/')){
  dir.create('./processed_data/')
}
save(list = setup_16s, file = './processed_data/setup_16s.rds')

### Identify and Count Engraftment Events ####

# Import the 16s ASV data
marker_lvl_16s <- read.csv("./data/16s_lvl_RA.csv")

# Import donor B 16s ASVs and create a vector
strains16s <- read.csv("./data/jcsz_16s_core_DonB.csv")
strains16s %>%
    pull() -> strains16s

# Set the cutoffs for 16s presence and absence
cutoff_pres_16S = 0.00013
cutoff_abs_16S = 0


### Get the engraftment matrix ####

# convert the 16s data to long format
long_mark_16s = mark_to_long(marker_lvl_16s, strains16s)
# Create a matrix of engraftment counts per patient
engr_16s = get_engraft(long_mark_16s, mapfile_16s, cutoff_abs_16S, 
                       cutoff_pres_16S, Treatment %in% c('FMT', 'Placebo'))

# Count engraftments per group
cts_16s_plt = count_engraft(engr_16s, tx_pv_16s, 'Treatment')
cts_16s_plt_rem = count_engraft(engr_16s, rs_pv_16s, 'Remission')

### Process Data for Visualization ####

#### Benchmarking and profiling. ####

# Using the same core DonorB ASVs as the engraftment uses: strains16s
donorB_16s_samples <- read.csv("./data/16s_lvl_RA_donorB_samples.csv")
donorB_16s_samples %>%
    gather(sample, abundance, -Marker) %>%
    # any marker with abundance > 0
    # markers that were present in at least one donor B sample
    filter(!is.na(abundance)) %>%
    filter(as.numeric(abundance) > 0) -> donorB_16s_samples

# profile the feature types:
profile_16s = get_profile_pat(long_mark_16s, mapfile_16s, cutoff_abs_16S, 
                              cutoff_pres_16S,
                              Treatment %in% c('FMT', 'Placebo'))
profile_16s_don = get_profile_don(long_mark_16s, mapfile_16s, cutoff_abs_16S,
                                  cutoff_pres_16S,
                                  Treatment %in% c('FMT', 'Placebo'))

### Save 16S Objects ####

eng_16s = c('engr_16s', 'cts_16s_plt', 'cts_16s_plt_rem')
            
save(list = eng_16s, file = './processed_data/eng_16s.rds')

fig_16s = c('donorB_16s_samples','profile_16s', 'profile_16s_don','cts_16s_plt',
            'cts_16s_plt_rem')
save(list = fig_16s, file = './processed_data/fig_16s.rds')

oth_16s = c('marker_lvl_16s', 'strains16s', 'long_mark_16s','cutoff_abs_16S', 
            'cutoff_pres_16S', 'long_mark_16s')
save(list = unique(c(eng_16s, fig_16s, oth_16s)),
     file = './processed_data/all_objs_16s.rds')

## Setup Metagenomics ####

# This mapfile is used for every non-16S feature type
mapfile_mgm <- read.csv("./data/UCFMT1_METAmap.txt", sep = "\t")

# These patient vectors are used for every non-16S feature type
tx_pv_mgm = pat_to_vect(mapfile_mgm, 'Fig_lab', 'Treatment')
rs_pv_mgm = pat_to_vect(filter(mapfile_mgm, Treatment == 'FMT'), 'Fig_lab', 
                        'Remission')

### Save Metagenomics Objects ####

mgm_setup = c('mapfile_mgm', 'tx_pv_mgm','rs_pv_mgm')
save(list = mgm_setup, file = './processed_data/mgm_setup.rds')

## Species ####

### Run Engraftment Functions ####

# Read in the species-level markers
marker_lvl_sp <- read.csv("./data/species_lvl_RA.csv")

# get a vector of donor B species
donorB_species_all = (marker_lvl_sp
                      %>% gather(sample, abundance, -Marker) 
                      #filtering donor B markers
                      %>% filter(str_detect(sample, "DonorB_D_")) 
                      # any marker with abundance > 0
                      # markers that were present in at least one donor B sample
                      %>% filter(!is.na(abundance)) 
                      %>% filter(as.numeric(abundance) > 0) 
                      %>% mutate(sample= gsub("DonorB_D_", "", sample)))

# Set the cutoffs and get the donorB markers
cutoff_pres_sp = 0.0073
cutoff_abs_sp = 0

# Filter to only include donor B species that meet that cutoff
donorB_species = (donorB_species_all 
                  %>% filter(as.numeric(abundance) > cutoff_pres_sp))

# Turn it into a vector
donorB_species %>%
    select(Marker) %>% distinct() %>%
    pull() -> strains_sp


### Get the wide matrix features engrafted (1) or not (0) in patients
long_mark_sp = mark_to_long(marker_lvl_sp, strains_sp)
engr_sp = get_engraft(long_mark_sp, mapfile_mgm, cutoff_abs_sp, cutoff_pres_sp,
                      Treatment %in% c('FMT', 'Placebo'))
cts_sp = count_engraft(engr_sp, tx_pv_mgm, 'Treatment')
cts_sp_rem = count_engraft(engr_sp, rs_pv_mgm, txrm = 'Remission')

### Process Data for Visualization ####

# profile the feature types:
profile_sp = get_profile_pat(long_mark_sp, mapfile_mgm, cutoff_abs_sp, 
                             cutoff_pres_sp, Treatment %in% c('FMT', 'Placebo'))
profile_sp_don = get_profile_don(long_mark_sp, mapfile_mgm, cutoff_abs_sp,
                                 cutoff_pres_sp, 
                                 Treatment %in% c('FMT', 'Placebo'))

### Save Species Engraftment ####

eng_sp = c('cts_sp', 'engr_sp', 'cts_sp_rem')
save(list = eng_sp, file = './processed_data/eng_sp.rds')
fig_sp = c('donorB_species_all', 'profile_sp',
           'profile_sp_don', 'cts_sp', 'cts_sp_rem')
save(list = fig_sp, file = './processed_data/fig_sp.rds')
oth_sp = c('marker_lvl_sp', 'cutoff_abs_sp', 'strains_sp', 'long_mark_sp', 
           'cutoff_pres_sp', 'donorB_species')
save(list = unique(c(eng_sp, fig_sp, oth_sp)), 
     file = './processed_data/all_objs_sp.rds')


## Strains ####

### Run the Engraftment Functions ####

# Import the strain-level marker data
marker_lvl_st <- read.csv("./data/marker_lvl_RA.csv")

# Get the strains that are found in donor B
donorB_strains_all = (marker_lvl_st 
                      %>% gather(sample, abundance, -sp_marker, -st_marker, -Marker) 
                      %>% filter(str_detect(sample, "DonorB_D_"))
                      #filtering donor B markers
                      %>% filter(!is.na(abundance)) 
                      %>% filter(as.numeric(abundance) > 0) 
                      %>% mutate(sample= gsub("DonorB_D_", "", sample)))

# Set the presence cutoff used in both donor B and patient samples
cutoff_pres_st = 0.00013
cutoff_abs_st = 0

# Filter donor B to only the strains that meet the cutoff
donorB_strains = (donorB_strains_all
                  %>% filter(as.numeric(abundance) > cutoff_pres_st))

# Turn it into a vector
donorB_strains %>%
    select(Marker) %>% distinct() %>%
    pull() -> strains_st

### Get the engraftment matrix

long_mark_st = mark_to_long(select(marker_lvl_st, -sp_marker, -st_marker), 
                            strains_st)
engr_st = get_engraft(long_mark_st, mapfile_mgm, cutoff_abs_st, cutoff_pres_st,
                      Treatment %in% c('FMT', 'Placebo'))
cts_st = count_engraft(engr_st, tx_pv_mgm, 'Treatment')
cts_st_rem = count_engraft(engr_st, rs_pv_mgm, 'Remission')

### Process Data for Visualization ####

# profile the feature types:
profile_st = get_profile_pat(long_mark_st, mapfile_mgm, cutoff_abs_st, cutoff_pres_st,
                             Treatment %in% c('FMT', 'Placebo'))
profile_st_don = get_profile_don(long_mark_st, mapfile_mgm, cutoff_abs_st,
                                 cutoff_pres_st, Treatment %in% c('FMT', 'Placebo'))

### Save Strain Engraftment ####

eng_st = c('cts_st', 'engr_st', 'cts_st_rem')
save(list = eng_st, file = './processed_data/eng_st.rds')
fig_st = c('donorB_strains_all', 'profile_st', 'cts_st', 'cts_st_rem',
           'profile_st_don')
save(list = fig_st, file = './processed_data/fig_st.rds')
oth_st = c('marker_lvl_st', 'strains_st', 'cutoff_pres_st', 'cutoff_abs_st', 
           'long_mark_st', 'donorB_strains')
save(list = unique(c(eng_st, fig_st, oth_st)),
     file = './processed_data/all_objs_st.rds')		   

## MAG ####

### Run Engraftment Functions ####

load('./data/MAGs.RData')
marker_lvl_mg_all = read.csv('./data/marker_lvl_mg.csv')
marker_lvl_mg = (marker_lvl_mg_all 
                 %>% filter(str_detect(sample, "PMCL")) 
                 %>% select(sample, bin_id, bin_coverage) 
                 %>% spread(sample, bin_coverage) 
                 %>% rename(Marker= bin_id))
# Set the presence and absence cutoffs for patients. DonorB doesn't need these
# because all the MAGs come from donor B
cutoff_pres_mg = 75
cutoff_abs_mg = 20

### Get the engraftment matrix ####

long_mark_mg = mark_to_long(marker_lvl_mg)
engr_mg = get_engraft(long_mark_mg, mapfile_mgm, cutoff_abs_mg, cutoff_pres_mg,
                      Treatment %in% c('FMT', 'Placebo'))
cts_mg = count_engraft(engr_mg, tx_pv_mgm, 'Treatment')
cts_mg_rem = count_engraft(engr_mg, rs_pv_mgm, 'Remission')

### Process Data for Visualization ####

# profile the feature types:
profile_mg = get_profile_pat(long_mark_mg, mapfile_mgm, cutoff_abs_mg, 
                             cutoff_pres_mg, Treatment %in% c('FMT', 'Placebo'))
profile_mg_don = get_profile_don(long_mark_mg, mapfile_mgm, cutoff_abs_mg, 
                                 cutoff_pres_mg,
                                 Treatment %in% c('FMT', 'Placebo'))

### Save MAG Objects ####

eng_mg = c('cts_mg', 'engr_mg', 'cts_mg_rem')
save(list = eng_mg, file = './processed_data/eng_mg.rds')
fig_mg = c('marker_lvl_mg_all', 'profile_mg', 'profile_mg_don', 'cts_mg',
           'cts_mg_rem', 'MAGs', 'engr_mg')
save(list = fig_mg, file = './processed_data/fig_mg.rds')
oft_mg = c('marker_lvl_mg', 'cutoff_pres_mg', 
           'cutoff_abs_mg', 'long_mark_mg')
save(list = unique(c(eng_mg, fig_mg, oft_mg)), 
     file = './processed_data/all_objs_mg.rds')

## Genes ####

### Run Engraftment Functions ####

#### Import the data ####

marker_lvl_ge <- read.csv("./data/genes_lvl.lrg")

# Set the presence and absence cutoffs for patients. These aren't used for donor
# B because all the genes came from the donor B assemblies.
cutoff_pres_ge = 90
cutoff_abs_ge = 25

### Get the engraftment matrix ####

long_mark_ge = mark_to_long(marker_lvl_ge)
engr_ge = get_engraft(long_mark_ge, mapfile_mgm, cutoff_abs_ge, cutoff_pres_ge,
                      Treatment %in% c('FMT', 'Placebo'))
cts_ge = count_engraft(engr_ge, tx_pv_mgm, 'Treatment')
cts_ge_rem = count_engraft(engr_ge, rs_pv_mgm, 'Remission')
```

### Process Data for Visualization ####

# profile the feature types:
profile_ge = get_profile_pat(long_mark_ge, mapfile_mgm, cutoff_abs_ge, cutoff_pres_ge,
                             Treatment %in% c('FMT', 'Placebo'))
profile_ge_don = get_profile_don(long_mark_ge, mapfile_mgm, cutoff_abs_ge, cutoff_pres_ge,
                                 Treatment %in% c('FMT','Placebo'))

### Save Gene Engraftment ####

eng_ge = c('cts_ge', 'engr_ge', 'cts_ge_rem')
save(list = eng_ge, file = './processed_data/eng_ge.rds')
fig_ge = c('profile_ge', 'profile_ge_don', 'cts_ge', 'cts_ge_rem', 'engr_ge')
save(list = fig_ge, file = './processed_data/fig_ge.rds')
oth_ge = c('marker_lvl_ge', 'cutoff_pres_ge', 'cutoff_abs_ge', 'long_mark_ge')
save(list = unique(c(eng_ge, fig_ge, oth_ge)), 
     file = './processed_data/all_objs_ge.rds')
