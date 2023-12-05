#!/usr/bin/env Rscript

##### Setup #####
print('Setup')

# Import packages, set working directory, source functions
library(tidyverse)
# setwd('~/Projects/Favours/Sharok/UCFMT1/scripts/')
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
mapfile <- read.csv("../data/UCFMT1_16Smap.txt")
tx_pv = pat_to_vect(mapfile, 'Fig_lab', 'Treatment')
rs_pv = pat_to_vect(filter(mapfile, Treatment == 'FMT'), 'Fig_lab', 
                          'Remission')

# Import the 16s ASV data
marker_lvl_16s <- read.csv("../data/16s_lvl_RA.csv")

# Import donor B 16s ASVs and create a vector
strains16s <- read.csv("../data/16s_core_DonB.csv")
strains16s %>%
  pull() -> strains16s

# Set the cutoffs for 16s presence and absence
cutoff_pres_16S = 0.0001
cutoff_abs_16S = 0


### Get the engraftment matrix ###

# convert the 16s data to long format
long_mark_16s = mark_to_long(marker_lvl_16s, strains16s)
# Create a matrix of engraftment counts per patient
engr_16s = get_engraft(long_mark_16s, mapfile, cutoff_abs_16S, cutoff_pres_16S,
                       Treatment %in% c('FMT', 'Placebo'))


#### 16S FMT vs Placebo ####
print('16S FMT vs Placebo')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

# The uneven numbers in the two groups make our test invalid. Subsample to fix
# that

# Find out how many of each group there are. It's 20
table(tx_pv)
nsamp = min(table(tx_pv))

# Create an array to store 100 count matrices
cts_16s_all = array(dim = c(1105, 2, 100), 
                    dimnames = list(rownames(engr_16s),
                                    c('tx_tot','bl_tot'),
                                    NULL))
# Resample (without permuting) from the Placebo group 100 times to get 100 count
# matrices
set.seed(4)
for (i in 1:100){
    # Subsample to 20
    tx_pv_sub = tx_pv[tx_pv == 'FMT']
    tx_pv_sub = c(tx_pv_sub, sample(tx_pv[tx_pv == 'Placebo'], nsamp))
    cts_16s_all[,,i] = count_engraft(engr_16s, tx_pv_sub, 'Treatment')
}

# Take the average of the 100 matrices to get the observed count
cts_16s = apply(cts_16s_all, c(1,2), mean)
cts_16s_sd = apply(cts_16s_all, c(1,2), sd)
# Calculate the observed value of the three test statistics from the averaged
# counts
obs_16s = get_stats(cts_16s)

#### Permute 
print('Permute')

# We are doing 2000 permutations in total. This includes 1999 permuted values
# plus our observed value in the null distribution to which we then compare the
# observed value.
nperm = 2000
perms_16s = do_permute(engr_16s, cts_16s, obs_16s, tx_pv, nperm, 
                       subsample = nsamp, txrm = 'Treatment')
perms_16s[['cts_sd_array']][,,1] = cts_16s_sd
# Save the permuted data
save(perms_16s, file = '../permut_data/intermed/perms_16s.RData')

# Calculate the p-values
pvals_16s = get_pvals(perms_16s[['stat_mat']])

# x2 because we're doing a two-tailed test
(pvals_16s = 2*pvals_16s)

# add the p-value to the matrix for plotting
pvals['FMTvPl_16s',] = round(pvals_16s, 5)

#### 16S Res vs NoRes ####
print('16S Res vs NoRes')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

# Subsample 100 times to get means
cts_16s_rem_all = array(dim = c(1105, 2, 100), 
                    dimnames = list(rownames(engr_16s),
                                    c('tx_tot','bl_tot'),
                                    NULL))

# Find out how many of each group there are. It's 6
table(rs_pv)
nsamp = min(table(rs_pv))

# Resample (without permuting) from the Placebo group 100 times to get 100 count
# matrices
for (i in 1:100){
    # Subsample to 20
    rs_pv_sub = rs_pv[rs_pv == 'Res']
    rs_pv_sub = c(rs_pv_sub, sample(rs_pv[rs_pv == 'NoRes'], nsamp))
    cts_16s_rem_all[,,i] = count_engraft(engr_16s, rs_pv_sub, 'Remission')
}

# Take the average of the 100 matrices to get the observed count
cts_16s_rem = apply(cts_16s_rem_all, c(1,2), mean)
cts_16s_rem_sd = apply(cts_16s_rem_all, c(1,2), sd)
# Calculate the observed value of the three test statistics from the averaged
# counts
obs_16s_rem = get_stats(cts_16s_rem)


#### Permute
print('Permute')

# We are doing 2000 permutations in total. This includes 1999 permuted values
# plus our observed value in the null distribution to which we then compare the
# observed value.
nperm = 2000
perms_16s_rem = do_permute(engr_16s, cts_16s_rem, obs_16s_rem, rs_pv, nperm,
                           subsample = nsamp, txrm = 'Remission')
perms_16s_rem[['cts_sd_array']][,,1] = cts_16s_rem_sd
# Save the permuted data
save(perms_16s_rem, file = '../permut_data/intermed/perms_16s_rem.RData')

# Calculate, double, and store the p-values
pvals_16s_rem = get_pvals(perms_16s_rem[['stat_mat']])
pvals_16s_rem
(pvals_16s_rem = 2 * pvals_16s_rem)
pvals['ResvNoRes_16s',] = round(pvals_16s_rem, 5)


#### Setup Metagenomics ####
print('Setup Metagenomics')

# This mapfile is used for every non-16S feature type
mapfile <- read.csv("../data/UCFMT1_METAmap.txt", sep = "\t")

# These patient vectors are used for every non-16S feature type
tx_pv = pat_to_vect(mapfile, 'Fig_lab', 'Treatment')
rs_pv = pat_to_vect(filter(mapfile, Treatment == 'FMT'), 'Fig_lab', 
                          'Remission')


#### Species Setup ####
print('Species Setup')

### Import the data

# Read in the species-level markers
marker_lvl_sp <- read.csv("../data/species_lvl_RA.csv")

# get a vector of donor B species
donorB_species = (marker_lvl_sp
                  %>% gather(sample, abundance, -Marker) 
                  #filtering donor B markers
                  %>% filter(str_detect(sample, "DonorB_D_")) 
                  # any marker with abundance > 0
                  # markers that were present in at least one donor B sample
                  %>% filter(!is.na(abundance)) 
                  %>% filter(as.numeric(abundance) > 0) 
                  %>% mutate(sample= gsub("DonorB_D_", "", sample)))

# Set the cutoff for "presence" with species. This is used for both patient
# samples and donor B
cutoff_pres_sp = 0.005

# Filter to only include donor B species that meet that cutoff
donorB_species = (donorB_species 
                  %>% filter(as.numeric(abundance) > cutoff_pres_sp))

# Turn it into a vector
donorB_species %>%
  select(Marker) %>% distinct() %>%
  pull() -> strains_sp

# Set the species absence cutoff
cutoff_abs_sp = 0

### Get the wide matrix features engrafted (1) or not (0) in patients
long_mark_sp = mark_to_long(marker_lvl_sp, strains_sp)
engr_sp = get_engraft(long_mark_sp, mapfile, cutoff_abs_sp, cutoff_pres_sp,
                       Treatment %in% c('FMT', 'Placebo'))


#### Species FMT vs Placebo ####
print('Species FMT vs Placebo')

### Calculate the observed test statistics
print('Calculate the observed test statistics')

# Get the engraftment counts and calculate the test statistics
cts_sp = count_engraft(engr_sp, tx_pv, 'Treatment')
obs_sp = get_stats(cts_sp)


### Permute 
print('Permute')

# Permute 1999 times and include the observed value in the null distribution
nperm = 2000
perms_sp = do_permute(engr_sp, cts_sp, obs_sp, tx_pv, nperm, 
                      txrm = 'Treatment')

# Save the permuted data so we don't have to do this all again
save(perms_sp, file = '../permut_data/intermed/perms_sp.RData')

# Calculate the p-values, double them, and store them
pvals_sp = get_pvals(perms_sp[['stat_mat']])
(pvals_sp = 2*pvals_sp)
pvals['FMTvPl_species',] = round(pvals_sp, 5)


#### Species Res vs NoRes ####
print('Species Res vs NoRes')

### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_sp_rem = count_engraft(engr_sp, rs_pv, txrm = 'Remission')
obs_sp_rem = get_stats(cts_sp_rem)

#### Permute
print('Permute')

# Permute 1999 times and include the observed value in the null distribution
nperm = 2000
perms_sp_rem = do_permute(engr_sp, cts_sp_rem, obs_sp_rem, rs_pv, nperm,
                          txrm = 'Remission')

# Save the permuted values
save(perms_sp_rem, file = '../permut_data/intermed/perms_sp_rem.RData')

# Calculate the p-values, double them, and store them
pvals_sp_rem = get_pvals(perms_sp_rem[['stat_mat']])
(pvals_sp_rem = 2*pvals_sp_rem)
pvals['ResvNoRes_species',] = round(pvals_sp_rem, 5)

#### Strains Setup ####
print('Strains Setup')

### Import the data

# Import the strain-level marker data
marker_lvl_st <- read.csv("../data/marker_lvl_RA.csv")

# Get the strains that are found in donor B
donorB_strains = (marker_lvl_st 
                  %>% gather(sample, abundance, -sp_marker, -st_marker, -Marker) 
                  %>% filter(str_detect(sample, "DonorB_D_"))
                  #filtering donor B markers
                  %>% filter(!is.na(abundance)) 
                  %>% filter(as.numeric(abundance) > 0) 
                  %>% mutate(sample= gsub("DonorB_D_", "", sample)))

# Set the presence cutoff used in both donor B and patient samples
cutoff_pres_st = 0.0001

# Filter donor B to only the strains that meet the cutoff
donorB_strains = (donorB_strains
                  %>% filter(as.numeric(abundance) > cutoff_pres_st))

# Turn it into a vector
donorB_strains %>%
  select(Marker) %>% distinct() %>%
  pull() -> strains_st

# Set the strain absence cutoff
cutoff_abs_st = 0

### Get the engraftment matrix

long_mark_st = mark_to_long(select(marker_lvl_st, -sp_marker, -st_marker), 
                            strains_st)
engr_st = get_engraft(long_mark_st, mapfile, cutoff_abs_st, cutoff_pres_st,
                       Treatment %in% c('FMT', 'Placebo'))


#### Strains FMT vs Placebo ####
print('Strains FMT vs Placebo')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_st = count_engraft(engr_st, tx_pv, 'Treatment')
obs_st = get_stats(cts_st)

#### Permute 
print('Permute')

# 1999 permutations plus the observed value go in the null distribution
nperm = 2000
perms_st = do_permute(engr_st, cts_st, obs_st, tx_pv, nperm, 
                      txrm = 'Treatment')

# Save the permuted values
save(perms_st, file = '../permut_data/intermed/perms_st.RData')

# Calculate, double, and store the p-values
pvals_st = get_pvals(perms_st[['stat_mat']])
(pvals_st = 2*pvals_st)
pvals['FMTvPl_strain',] = round(pvals_st, 5)

#### Strains Res vs NoRes ####
print('Strains Res vs NoRes')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_st_rem = count_engraft(engr_st, rs_pv, txrm = 'Remission')
obs_st_rem = get_stats(cts_st_rem)

#### Permute
print('Permute')

nperm = 2000
perms_st_rem = do_permute(engr_st, cts_st_rem, obs_st_rem, rs_pv, nperm,
                          txrm = 'Remission')

# Save the permuted values
save(perms_st_rem, file = '../permut_data/intermed/perms_st_rem.RData')

# calculate, double, and store the p-values
pvals_st_rem = get_pvals(perms_st_rem[['stat_mat']])
(pvals_st_rem = 2*pvals_st_rem)
pvals['ResvNoRes_strain',] = round(pvals_st_rem, 5)

#### MAGs Setup ####
print('MAGs Setup')

### Import the data

# Import the MAG and bin data
c_inBin <- read.csv("../data/Assembly/contigs_inBins.txt", sep = "\t", header = F)
c_inMAG <- read.csv("../data/Assembly/contigs_inMAGs.txt", sep = "\t", header = F)
colnames(c_inBin) <- c("contig", "bin_id")
colnames(c_inMAG) <- c("contig", "bin_id")

# Get a vector of the mags
MAGs = (c_inMAG 
        %>% select(bin_id) 
        %>% distinct() 
        %>% pull())


## read in the MAG coverage
heads_cover <- c("sample.id",
                 "contig", "startpos", "endpos", 
                 "numreads", "covbases", "coverage", "meandepth",
                 "meanbaseq", "meanmapq",
                 "ReadSampling", "Quality")

path="../data/perfect_mapping/DBD"
patt="_perfect_.*.cover"
marker_lvl_mg = (data.frame(sample.id = paste(dir(path,
                                                  pattern = patt))) 
             # read each file from the directory (current path) and then unnest
                 %>% mutate(file_contents = map(sample.id, ~ 
                                                    read.csv(file.path(path, .),
                                                         sep = "\t", header = F,
                                                           comment.char = "#"))) 
                 %>% unnest() 
                 %>%  mutate(sample.id = gsub("_perfect_BDB.cover", 
                                              "", sample.id)))
colnames(marker_lvl_mg) <- heads_cover

# Join the coverage info with the bin info, filter to only keep mags, calculate
# normalized coverage
marker_lvl_mg = (marker_lvl_mg 
                 %>% left_join(c_inBin) 
                 %>% mutate(bin_q= case_when(bin_id %in% MAGs ~ paste("MAG"),
                                             TRUE ~ paste("Bin"))) 
                 %>% filter(bin_q == "MAG") 
                 %>% select(-bin_q) 
                 %>% dplyr::rename(sample=sample.id) 
                 %>% group_by(sample, bin_id) 
                 %>% summarise(bin_len= sum(endpos),
                               bin_covbases=sum(covbases),
                               bin_numreads=sum(numreads),
                               bin_coverage=bin_covbases/bin_len*100))

# Keep only the patient samples, reshape the numbers, and name the Marker column
marker_lvl_mg = (marker_lvl_mg 
                 %>% filter(str_detect(sample, "PMCL")) 
                 %>% select(sample, bin_id, bin_coverage) 
                 %>% spread(sample, bin_coverage) 
                 %>% rename(Marker= bin_id))

# Set the presence and absence cutoffs for patients. DonorB doesn't need these
# because all the MAGs come from donor B
cutoff_pres_mg = 90
cutoff_abs_mg = 25


### Get the engraftment matrix

long_mark_mg = mark_to_long(marker_lvl_mg)
engr_mg = get_engraft(long_mark_mg, mapfile, cutoff_abs_mg, cutoff_pres_mg,
                       Treatment %in% c('FMT', 'Placebo'))


#### MAGs FMT vs Placebo ####
print('MAGs FMT vs Placebo')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_mg = count_engraft(engr_mg, tx_pv, 'Treatment')
obs_mg = get_stats(cts_mg)

#### Permute 
print('Permute')

nperm = 2000
perms_mg = do_permute(engr_mg, cts_mg, obs_mg, tx_pv, nperm,
                      txrm = 'Treatment')
# save the permuted data
save(perms_mg, file = '../permut_data/intermed/perms_mg.RData')

# calculate, double, and store the p-values
pvals_mg = get_pvals(perms_mg[['stat_mat']])
(pvals_mg = 2*pvals_mg)
pvals['FMTvPl_mags',] = round(pvals_mg, 5)

#### MAGs Res vs NoRes ####
print('MAGs Res vs NoRes')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_mg_rem = count_engraft(engr_mg, rs_pv, txrm = 'Remission')
obs_mg_rem = get_stats(cts_mg_rem)

#### Permute
print('Permute')

nperm = 2000
perms_mg_rem = do_permute(engr_mg, cts_mg_rem, obs_mg_rem, rs_pv, nperm,
                          txrm = 'Remission')
# save the permutations
save(perms_mg_rem, file = '../permut_data/intermed/perms_mg_rem.RData')

# calculate, double, and store the p-values
pvals_mg_rem = get_pvals(perms_mg_rem[['stat_mat']])
(pvals_mg_rem = 2*pvals_mg_rem)
pvals['ResvNoRes_mags',] = round(pvals_mg_rem, 5)


#### Genes Setup ####
print('Genes Setup')

### Import the data

marker_lvl_ge <- read.csv("../data/genes_lvl.csv")

# Set the presence and absence cutoffs for patients. These aren't used for donor
# B because all the genes came from the donor B assemblies.
cutoff_pres_ge = 90
cutoff_abs_ge = 5

### Get the engraftment matrix

long_mark_ge = mark_to_long(marker_lvl_ge)
engr_ge = get_engraft(long_mark_ge, mapfile, cutoff_abs_ge, cutoff_pres_ge,
                       Treatment %in% c('FMT', 'Placebo'))


#### Genes FMT vs Placebo ####
print('Genes FMT vs Placebo')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_ge = count_engraft(engr_ge, tx_pv, 'Treatment')
obs_ge = get_stats(cts_ge)

#### Permute 
print('Permute')

# This uses tens of GB of RAM and should be done on alpsr only
nperm = 2000
perms_ge = do_permute(engr_ge, cts_ge, obs_ge, tx_pv, nperm,
                      txrm = 'Treatment')

# Save the permutations
save(perms_ge, file = '../permut_data/intermed/perms_ge.RData')

# Calculate, double, and store the p-values
pvals_ge = get_pvals(perms_ge[['stat_mat']])
(pvals_ge = 2*pvals_ge)
pvals['FMTvPl_genes',] = round(pvals_ge, 5)

#### Genes Res vs NoRes ####
print('Genes Res vs NoRes')

#### Calculate the observed test statistics
print('Calculate the observed test statistics')

cts_ge_rem = count_engraft(engr_ge, rs_pv, txrm = 'Remission')
obs_ge_rem = get_stats(cts_ge_rem)

#### Permute
print('Permute')

# This uses tens of GB of RAM and should be done on alpsr only
nperm = 20
perms_ge_rem = do_permute(engr_ge, cts_ge_rem, obs_st_rem, rs_pv, nperm,
                          txrm = 'Remission')
# save the permutations
save(perms_ge_rem, file = '../permut_data/intermed/perms_ge_rem.RData')
# calculated, double, and store the p-values
pvals_ge_rem = get_pvals(perms_ge_rem[['stat_mat']])
(pvals_ge_rem = 2*pvals_ge_rem)
pvals['ResvNoRes_genes',] = round(pvals_ge_rem, 5)

# Write the p-value table
print('Write the p-value table')

write.csv(pvals, file = '../results/permuation_pvals.csv', 
          row.names = TRUE, col.names = NA)
