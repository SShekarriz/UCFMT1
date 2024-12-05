#!/usr/bin/env Rscript

##### Setup #####
print('Setup')

# Import packages, set working directory, source functions
set.seed(4)
library(tidyverse)
theme_set(theme_bw())
# setwd('~/Projects/Favours/Sharok/UCFMT1/scripts/')
source('./permutation_functions.R')
source('./sensitivity_functions.R')

# Create the p-value storage matrix
pvals = matrix(nrow = 10, ncol = 3)
colnames(pvals) = c('PEngraftedFeatures', 'PEngraftmentEvents', 
                     'PWeightedEvents')
rownames(pvals) = c('FMTvPl_16s','ResvNoRes_16s',
                    'FMTvPl_species','ResvNoRes_species',
                    'FMTvPl_strain','ResvNoRes_strain',
                    'FMTvPl_mags','ResvNoRes_mags',
                    'FMTvPl_genes','ResvNoRes_genes')


#### Setup Metagenomics ####
print('Setup Metagenomics')

load('../processed_data/mgm_setup.rds')

#### Genes Setup ####
print('Genes Setup')

load('../processed_data/eng_ge.rds')
load('../sensitivity/intermed/engr_arr.RData')
load('../sensitivity/intermed/cts_arrs.RData')

#### Genes FMT vs Placebo ####
print('Genes FMT vs Placebo')

#### Calculate the observed test statistics
print('Calculate the observed test statistics for each permutation')

# make a place to put them
n = dim(arr_cts)[1]
obs_sn = matrix(nrow = n, ncol = 3)
for (i in 1:n){
    obs_sn[1,] = get_stats(arr_cts[i,,])
}
rownames(obs_sn) = dimnames(arr_cts)$replicate
obs_ge = get_stats(cts_ge)

#### Genes Plot FMT vs Placebo ####
plt_all = sens_fig_13_df(arr_cts, 'Treatment')
sens_plt = sens_plot_fig_13(plt_all, 'Treatment')
sens_plt
ggsave('../sensitivity/intermed/f1_sens_plt.png', sens_plt,
       height = 10, width = 12)

f2_df = sens_fig_2_df(arr_cts)
f2_ge = sens_plot_fig2(f2_df)
f2_ge
ggsave('../sensitivity/f2_sens_plt.png',f2_ge, height = 10, width = 12)

#### Permute 
print('Permute')

# This uses tens of GB of RAM and should be done on alpsr only

# Iterate through the replicates to do the permutations
nperm = 2000
perm_lst = rep_perms(engr_arr, arr_cts, nperm, tx_pv_mgm, 'Treatment')
save(perm_lst, file = '../sensitivity/intermed/perm_lst.RData')
perm_df = sens_perm_df(perm_lst[['stat_arr']])
perm_fx = sens_perm_plot(perm_df, 'fx')
perm_xfx = sens_perm_plot(perm_df, 'xfx')
perm_x2fx = sens_perm_plot(perm_df, 'x2fx')

perm_fx
perm_xfx
perm_x2fx

ggsave('../sensitivity/intermed/perm_fx.png', perm_fx, 
       height = 10, width = 12)
ggsave('../sensitivity/intermed/perm_xfx.png', perm_xfx, 
       height = 10, width = 12)
ggsave('../sensitivity/intermed/perm_x2fx.png', perm_x2fx, 
       height = 10, width = 12)


# PICK UP HERE


# # Look at the distribution
# pl_pv_ge = plot_permutation(perms_ge[['stat_mat']])
# pl_pv_ge
# ggsave('../plots/pl_pv_ge.pdf',pl_pv_ge, height = 10, width = 5.5)
# 
# #### Genes Res vs NoRes ####
# print('Genes Res vs NoRes')
# 
# #### Calculate the observed test statistics
# print('Calculate the observed test statistics')
# 
# obs_ge_rem = get_stats(cts_ge_rem)
# 
# #### Genes Plot Res vs NoRes
# f3_ge = plot_fig_13(cts_ge_rem, 'Remission')
# f3_ge
# ggsave('../plots/f3_ge.pdf',f3_ge, height = 5, width = 5.5)
# 
# #### Permute
# print('Permute')
# 
# # This uses tens of GB of RAM and should be done on alpsr only
# nperm = 2000
# perms_ge_rem = do_permute(engr_ge, cts_ge_rem, obs_st_rem, rs_pv_mgm, nperm,
#                           txrm = 'Remission')
# head(perms_ge_rem[['stat_mat']])
# # save the permutations
# save(perms_ge_rem, file = '../permut_data/intermed/perms_ge_rem.RData')
# # load('../permut_data/intermed/perms_ge_rem.RData')
# # calculated, double, and store the p-values
# (pvals_ge_rem = get_pvals(perms_ge_rem[['stat_mat']]))
# pvals['ResvNoRes_genes',] = round(pvals_ge_rem, 5)
# 
# # Look at the distribution
# pl_pv_ge_rem = plot_permutation(perms_ge_rem[['stat_mat']])
# pl_pv_ge_rem
# ggsave('../plots/pl_pv_ge_rem.pdf',pl_pv_ge_rem, height = 10, width = 5.5)
# 
# # Write the p-value table
# print('Write the p-value table')
# 
# write.csv(pvals, file = '../results/permuation_pvals.csv', 
#           row.names = TRUE, col.names = NA)
