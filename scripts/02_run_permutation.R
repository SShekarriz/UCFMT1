#!/usr/bin/env Rscript

# Run Permutations ####

# This must be run from inside the top-level directory. I recommend running it 
# using Rscript, because it can take a long time.

## Setup #####
print('Setup')

# Import packages, set working directory, source functions
set.seed(4)
library(tidyverse)
theme_set(theme_bw())
source('./scripts/permutation_functions.R')

# Create the p-value storage matrix
pvals = matrix(nrow = 10, ncol = 3)
colnames(pvals) = c('PEngraftedFeatures', 'PEngraftmentEvents', 
                     'PWeightedEvents')
rownames(pvals) = c('FMTvPl_16s','ResvNoRes_16s',
                    'FMTvPl_species','ResvNoRes_species',
                    'FMTvPl_strain','ResvNoRes_strain',
                    'FMTvPl_mags','ResvNoRes_mags',
                    'FMTvPl_genes','ResvNoRes_genes')


## 16S Setup #####
print('16S Setup')
load('./processed_data/setup_16s.rds')

### Import the data ####

load('./processed_data/eng_16s.rds')

### 16S FMT vs Placebo ####
print('16S FMT vs Placebo')

# Calculate the observed test statistics
print('Calculate the observed test statistics')

# The uneven numbers in the two groups make our test invalid. Subsample to fix
# that

# Find out how many of each group there are. It's 20
table(tx_pv_16s)
nsamp = min(table(tx_pv_16s))

# Create an array to store 100 count matrices
cts_16s_all = array(dim = c(nrow(engr_16s), 2, 100),
                    dimnames = list(rownames(engr_16s),
                                    c('tx_tot','bl_tot'),
                                    NULL))
# Create a matrix to store the 100 test statistics
obs_16s_all = matrix(nrow = 100, ncol = 3)

# Resample (without permuting) from the Placebo group 100 times to get 100 count
# matrices
for (i in 1:100){
    # Subsample to 20
    tx_pv_sub = tx_pv_16s[tx_pv_16s == 'FMT']
    tx_pv_sub = c(tx_pv_sub, sample(tx_pv_16s[tx_pv_16s == 'Placebo'], nsamp))
    cts_16s_all[,,i] = count_engraft(engr_16s, tx_pv_sub, 'Treatment')
    obs_16s_all[i,] = get_stats(cts_16s_all[,,i])
}

# Calculate the means of the 100 subsampled test statistics
obs_16s = colMeans(obs_16s_all)

cts_16s = apply(cts_16s_all, c(1,2), mean)

#### Plot ####

# I think the only sensible way to plot this is to plot the full data and just
# talk about what we did for the statistics in the methods section. The counting
# is done in prep_data.R

# Fig 1 
f1_16s = plot_fig_13(cts_16s_plt, 'Treatment')
ggsave('./plots/f1_16s.pdf',f1_16s, height = 5, width = 5.5)
f1_16s
       

# Fig 2 
f2_16s = plot_fig_2(cts_16s_plt)
ggsave('./plots/f2_16s.pdf',f2_16s, height = 5, width = 5.5)
f2_16s

#### Permute ####
print('Permute')

# We are doing 2000 permutations in total. This includes 1999 permuted values
# plus our observed value in the null distribution to which we then compare the
# observed value.
nperm = 2000
perms_16s = do_permute(engr_16s, cts_16s, obs_16s, tx_pv_16s, nperm,
                       subsample = nsamp, txrm = 'Treatment')
head(perms_16s[['stat_mat']])
# Save the permuted data
save(perms_16s, file = './permut_data/intermed/perms_16s.RData')
# load('./permut_data/intermed/perms_16s.RData')

# Calculate the p-values
(pvals_16s = get_pvals(perms_16s[['stat_mat']]))

# add the p-value to the matrix for plotting
pvals['FMTvPl_16s',] = round(pvals_16s, 5)

pl_pv_16s = plot_permutation(perms_16s[['stat_mat']])
pl_pv_16s
ggsave('./plots/pl_pv_16s.pdf',pl_pv_16s, height = 10, width = 5.5)

### 16S Res vs. NoRes ####

print('16S Res vs NoRes')

#### Plot ####

# I think the only sensible way to plot this is to plot the full data and just
# talk about what we did for the statistics in the methods section. The counts
# is done in prep_data.R.

# Fig 1 
f3_16s = plot_fig_13(cts_16s_plt_rem, 'Remission')
f3_16s
ggsave('./plots/f3_16s.pdf',f3_16s, height = 5, width = 5.5)

#### Subsample ####

# Calculate the observed test statistics
print('Calculate the observed test statistics')

# Subsample 100 times to get means
cts_16s_rem_all = array(dim = c(nrow(engr_16s), 2, 100), 
                    dimnames = list(rownames(engr_16s),
                                    c('tx_tot','bl_tot'),
                                    NULL))
# Create a matrix to store the 100 test statistics
obs_16s_rem_all = matrix(nrow = 100, ncol = 3)

# Find out how many of each group there are. It's 6
table(rs_pv_16s)
nsamp = min(table(rs_pv_16s))


# Resample (without permuting) from the Placebo group 100 times to get 100 count
# matrices
for (i in 1:100){
    # Subsample to 20
    rs_pv_sub = rs_pv_16s[rs_pv_16s == 'Res']
    rs_pv_sub = c(rs_pv_sub, sample(rs_pv_16s[rs_pv_16s == 'NoRes'], nsamp))
    cts_16s_rem_all[,,i] = count_engraft(engr_16s, rs_pv_sub, 'Remission')
    obs_16s_rem_all[i,] = get_stats(cts_16s_rem_all[,,i])
}

# Take the average of the 100 matrices to get the observed count
cts_16s_rem = apply(cts_16s_rem_all, c(1,2), mean)
# Calculate the observed value of the three test statistics from the averaged
# counts
obs_16s_rem = colMeans(obs_16s_rem_all)

#### Permute ####
print('Permute')

# We are doing 2000 permutations in total. This includes 1999 permuted values
# plus our observed value in the null distribution to which we then compare the
# observed value.
nperm = 2000
perms_16s_rem = do_permute(engr_16s, cts_16s_rem, obs_16s_rem, rs_pv_16s, nperm,
                           subsample = nsamp, txrm = 'Remission')
# Save the permuted data
save(perms_16s_rem, file = './permut_data/intermed/perms_16s_rem.RData')
# load('./permut_data/intermed/perms_16s_rem.RData')
head(perms_16s_rem[['stat_mat']])

# Calculate and store the p-values
(pvals_16s_rem = get_pvals(perms_16s_rem[['stat_mat']]))
pvals['ResvNoRes_16s',] = round(pvals_16s_rem, 5)

pl_pv_16s_rem = plot_permutation(perms_16s_rem[['stat_mat']])
pl_pv_16s_rem
ggsave('./plots/pl_pv_16s_rem.pdf',pl_pv_16s_rem, height = 10, width = 5.5)


## Setup Metagenomics ####

print('Setup Metagenomics')
load('./processed_data/mgm_setup.rds')

## Species ####

### Setup ####

print('Species Setup')

### Import the data ####
load('./processed_data/eng_sp.rds')

### Species FMT vs Placebo ####
print('Species FMT vs Placebo')

# Calculate the observed test statistics
print('Calculate the observed test statistics')

# Get the engraftment counts and calculate the test statistics
obs_sp = get_stats(cts_sp)


#### Plot ####

# Fig 1
f1_sp = plot_fig_13(cts_sp, 'Treatment')
f1_sp
ggsave('./plots/f1_sp.pdf',f1_sp, height = 5, width = 5.5)

# Fig 2 
f2_sp = plot_fig_2(cts_sp)
f2_sp
ggsave('./plots/f2_sp.pdf',f2_sp, height = 5, width = 5.5)

#### Permute ####
print('Permute')

# Permute 1999 times and include the observed value in the null distribution
nperm = 2000
perms_sp = do_permute(engr_sp, cts_sp, obs_sp, tx_pv_mgm, nperm, 
                      txrm = 'Treatment')
head(perms_sp[['stat_mat']])

# Save the permuted data so we don't have to do this all again
save(perms_sp, file = './permut_data/intermed/perms_sp.RData')
#  load('./permut_data/intermed/perms_sp.RData')

# Calculate the p-values, double them, and store them
(pvals_sp = get_pvals(perms_sp[['stat_mat']]))
pvals['FMTvPl_species',] = round(pvals_sp, 5)

# Look at the distribution
pl_pv_sp = plot_permutation(perms_sp[['stat_mat']])
pl_pv_sp
ggsave('./plots/pl_pv_sp.pdf',pl_pv_sp, height = 10, width = 5.5)


### Species Res vs NoRes ####
print('Species Res vs NoRes')

# Calculate the observed test statistics
print('Calculate the observed test statistics')

obs_sp_rem = get_stats(cts_sp_rem)

#### Plot ####
f3_sp = plot_fig_13(cts_sp_rem, 'Remission')
f3_sp
ggsave('./plots/f3_sp.pdf',f3_sp, height = 5, width = 5.5)

#### Permute ####
print('Permute')

# Permute 1999 times and include the observed value in the null distribution
nperm = 2000
perms_sp_rem = do_permute(engr_sp, cts_sp_rem, obs_sp_rem, rs_pv_mgm, nperm,
                          txrm = 'Remission')
head(perms_sp_rem[['stat_mat']])

# Save the permuted values
save(perms_sp_rem, file = './permut_data/intermed/perms_sp_rem.RData')
# load('./permut_data/intermed/perms_sp_rem.RData')

# Calculate the p-values, double them, and store them
(pvals_sp_rem = get_pvals(perms_sp_rem[['stat_mat']]))
pvals['ResvNoRes_species',] = round(pvals_sp_rem, 5)

# Look at the distribution
pl_pv_sp_rem = plot_permutation(perms_sp_rem[['stat_mat']])
pl_pv_sp_rem
ggsave('./plots/pl_pv_sp_rem.pdf',pl_pv_sp_rem, height = 10, width = 5.5)

## Strains ####

### Setup ####
print('Strains Setup')

# Import the data
load('./processed_data/eng_st.rds')

### Strains FMT vs Placebo ####
print('Strains FMT vs Placebo')

# Calculate the observed test statistics
print('Calculate the observed test statistics')

obs_st = get_stats(cts_st)

#### Plot ####

f1_st = plot_fig_13(cts_st, 'Treatment')
f1_st
ggsave('./plots/f1_st.pdf',f1_st, height = 5, width = 5.5)

f2_st = plot_fig_2(cts_st)
f2_st
ggsave('./plots/f2_st.pdf',f2_st, height = 5, width = 5.5)

#### Permute ####
print('Permute')

# 1999 permutations plus the observed value go in the null distribution
nperm = 2000
perms_st = do_permute(engr_st, cts_st, obs_st, tx_pv_mgm, nperm, 
                      txrm = 'Treatment')
head(perms_st[['stat_mat']])

# Save the permuted values
save(perms_st, file = './permut_data/intermed/perms_st.RData')
# load('./permut_data/intermed/perms_st.RData')

# Calculate, double, and store the p-values
(pvals_st = get_pvals(perms_st[['stat_mat']]))
pvals['FMTvPl_strain',] = round(pvals_st, 5)

# Look at the distribution
pl_pv_st = plot_permutation(perms_st[['stat_mat']])
pl_pv_st
ggsave('./plots/pl_pv_st.pdf',pl_pv_st, height = 10, width = 5.5)

### Strains Res vs NoRes ####
print('Strains Res vs NoRes')

### Calculate the observed test statistics
print('Calculate the observed test statistics')

obs_st_rem = get_stats(cts_st_rem)

#### Plot ####

f3_st = plot_fig_13(cts_st_rem, 'Remission')
f3_st
ggsave('./plots/f3_st.pdf',f3_st, height = 5, width = 5.5)

#### Permute ####
print('Permute')

nperm = 2000
perms_st_rem = do_permute(engr_st, cts_st_rem, obs_st_rem, rs_pv_mgm, nperm,
                          txrm = 'Remission')
head(perms_st_rem[['stat_mat']])

# Save the permuted values
save(perms_st_rem, file = './permut_data/intermed/perms_st_rem.RData')
# load('./permut_data/intermed/perms_st_rem.RData')

# calculate, double, and store the p-values
(pvals_st_rem = get_pvals(perms_st_rem[['stat_mat']]))
pvals['ResvNoRes_strain',] = round(pvals_st_rem, 5)

# Look at the distribution
pl_pv_st_rem = plot_permutation(perms_st_rem[['stat_mat']])
pl_pv_st_rem
ggsave('./plots/pl_pv_st_rem.pdf',pl_pv_st_rem, height = 10, width = 5.5)

## MAGs ####

### Setup ####
print('MAGs Setup')

# Import the data

load('./processed_data/eng_mg.rds')

### MAGs FMT vs Placebo ####
print('MAGs FMT vs Placebo')

# Calculate the observed test statistics
print('Calculate the observed test statistics')

obs_mg = get_stats(cts_mg)

#### Plot ####

f1_mg = plot_fig_13(cts_mg, 'Treatment')
f1_mg
ggsave('./plots/f1_mg.pdf',f1_mg, height = 5, width = 5.5)

f2_mg = plot_fig_2(cts_mg)
f2_mg
ggsave('./plots/f2_mg.pdf',f2_mg, height = 5, width = 5.5)

#### Permute ####

print('Permute')

nperm = 2000
perms_mg = do_permute(engr_mg, cts_mg, obs_mg, tx_pv_mgm, nperm,
                      txrm = 'Treatment')
head(perms_mg[['stat_mat']])
# save the permuted data
save(perms_mg, file = './permut_data/intermed/perms_mg.RData')
# load('./permut_data/intermed/perms_mg.RData')

# calculate, double, and store the p-values
(pvals_mg = get_pvals(perms_mg[['stat_mat']]))
pvals['FMTvPl_mags',] = round(pvals_mg, 5)

# Look at the distribution
pl_pv_mg = plot_permutation(perms_mg[['stat_mat']])
pl_pv_mg
ggsave('./plots/pl_pv_mg.pdf',pl_pv_mg, height = 10, width = 5.5)

### MAGs Res vs NoRes ####
print('MAGs Res vs NoRes')

# Calculate the observed test statistics
print('Calculate the observed test statistics')

obs_mg_rem = get_stats(cts_mg_rem)

#### Plot ####

f3_mg = plot_fig_13(cts_mg_rem,'Remission')
f3_mg
ggsave('./plots/f3_mg.pdf',f3_mg, height = 5, width = 5.5)

#### Permute ####
print('Permute')

nperm = 2000
perms_mg_rem = do_permute(engr_mg, cts_mg_rem, obs_mg_rem, rs_pv_mgm, nperm,
                          txrm = 'Remission')
head(perms_mg_rem[['stat_mat']])
# save the permutations
save(perms_mg_rem, file = './permut_data/intermed/perms_mg_rem.RData')
# load('./permut_data/intermed/perms_mg_rem.RData')

# calculate, double, and store the p-values
(pvals_mg_rem = get_pvals(perms_mg_rem[['stat_mat']]))
pvals['ResvNoRes_mags',] = round(pvals_mg_rem, 5)

# Look at the distribution
pl_pv_mg_rem = plot_permutation(perms_mg_rem[['stat_mat']])
pl_pv_mg_rem
ggsave('./plots/pl_pv_mg_rem.pdf',pl_pv_mg_rem, height = 10, width = 5.5)


## Genes ####

### Setup ####
print('Genes Setup')

load('./processed_data/eng_ge.rds')

### Genes FMT vs Placebo ####
print('Genes FMT vs Placebo')

# Calculate the observed test statistics
print('Calculate the observed test statistics')

obs_ge = get_stats(cts_ge)

#### Plot ####

f1_ge = plot_fig_13(cts_ge, 'Treatment')
f1_ge
ggsave('./plots/f1_ge.pdf',f1_ge, height = 5, width = 5.5)

f2_ge = plot_fig_2(cts_ge)
f2_ge
ggsave('./plots/f2_ge.pdf',f2_ge, height = 5, width = 5.5)

#### Permute ####
print('Permute')

# This uses tens of GB of RAM and should only be done on a server that can 
# handle it
nperm = 2000
perms_ge = do_permute(engr_ge, cts_ge, obs_ge, tx_pv_mgm, nperm,
                      txrm = 'Treatment')
head(perms_ge[['stat_mat']])

# Save the permutations
save(perms_ge, file = './permut_data/intermed/perms_ge.RData')
# load('./permut_data/intermed/perms_ge.RData')

# Calculate, double, and store the p-values
(pvals_ge = get_pvals(perms_ge[['stat_mat']]))
pvals['FMTvPl_genes',] = round(pvals_ge, 5)

# Look at the distribution
pl_pv_ge = plot_permutation(perms_ge[['stat_mat']])
pl_pv_ge
ggsave('./plots/pl_pv_ge.pdf',pl_pv_ge, height = 10, width = 5.5)

### Genes Res vs NoRes ####
print('Genes Res vs NoRes')

# Calculate the observed test statistics
print('Calculate the observed test statistics')

obs_ge_rem = get_stats(cts_ge_rem)

#### Plot ####

f3_ge = plot_fig_13(cts_ge_rem, 'Remission')
f3_ge
ggsave('./plots/f3_ge.pdf',f3_ge, height = 5, width = 5.5)

#### Permute ####
print('Permute')

# This uses tens of GB of RAM and should be done on alpsr only
nperm = 2000
perms_ge_rem = do_permute(engr_ge, cts_ge_rem, obs_st_rem, rs_pv_mgm, nperm,
                          txrm = 'Remission')
head(perms_ge_rem[['stat_mat']])
# save the permutations
save(perms_ge_rem, file = './permut_data/intermed/perms_ge_rem.RData')
# load('./permut_data/intermed/perms_ge_rem.RData')
# calculated, double, and store the p-values
(pvals_ge_rem = get_pvals(perms_ge_rem[['stat_mat']]))
pvals['ResvNoRes_genes',] = round(pvals_ge_rem, 5)

# Look at the distribution
pl_pv_ge_rem = plot_permutation(perms_ge_rem[['stat_mat']])
pl_pv_ge_rem
ggsave('./plots/pl_pv_ge_rem.pdf',pl_pv_ge_rem, height = 10, width = 5.5)

# Write the p-value table
print('Write the p-value table')

write.csv(pvals, file = './results/permuation_pvals.csv', 
          row.names = TRUE, col.names = NA)
