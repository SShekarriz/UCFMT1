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
print('Calculate the observed test statistics for each replicate')

# make a place to put them
n = dim(arr_cts)[1]
obs_sn = matrix(nrow = n, ncol = 3)
for (i in 1:n){
    obs_sn[i,] = get_stats(arr_cts[i,,])
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
# nperm = 2000
# perm_lst = rep_perms(engr_arr, arr_cts, nperm, tx_pv_mgm, 'Treatment')
# save(perm_lst, file = '../sensitivity/intermed/perm_lst.RData')
# perm_df = sens_perm_df(perm_lst[['stat_arr']])
# perm_fx = sens_perm_plot(perm_df, 'fx')
# perm_xfx = sens_perm_plot(perm_df, 'xfx')
# perm_x2fx = sens_perm_plot(perm_df, 'x2fx')
# 
# perm_fx
# perm_xfx
# perm_x2fx
# 
# ggsave('../sensitivity/intermed/perm_fx.png', perm_fx, 
#        height = 10, width = 12)
# ggsave('../sensitivity/intermed/perm_xfx.png', perm_xfx, 
#        height = 10, width = 12)
# ggsave('../sensitivity/intermed/perm_x2fx.png', perm_x2fx, 
#        height = 10, width = 12)

load('../sensitivity/intermed/perm_lst.RData')

# Get the observed p-values

load('../permut_data/intermed/perms_ge.RData')
obs_pvals_ge = get_pvals(perms_ge[['stat_mat']])

# Arrange the data for plotting

pval_arr = perm_lst[['pval_arr']]
colnames(pval_arr) = names(obs_pvals_ge)
head(pval_arr)

reps = dimnames(perm_lst[[1]])$replicate
cutoffs = get_reps(reps)
pval_df = (pval_arr
           %>% data.frame()
           %>% mutate(Replicate = reps,
                      AbsCutoff = as.numeric(cutoffs[['abs']]),
                      PresCutoff = as.numeric(cutoffs[['pres']]))
           %>% pivot_longer(c(-Replicate, -AbsCutoff, -PresCutoff), 
                            names_to  = 'Statistic', 
                            values_to = 'P')
           %>% mutate(Observed = case_when(AbsCutoff == 25 & PresCutoff == 90 ~
                                               'Used',
                                           TRUE ~ 'Sens')))
head(pval_df)

sens_pval = ggplot(pval_df, aes(P)) +
    geom_histogram(bins = 4) +
    facet_wrap(~Statistic) +
    geom_linerange(aes(ymin = 0, ymax = 50),
               data = filter(pval_df, Observed == 'Used'),
               colour = 'red') +
    xlab('Permutation p-value') +
    scale_x_continuous(labels = calc_lab)
sens_pval

ggsave('../sensitivity/sens_pval.png', sens_pval)

#### Genes Res vs NoRes ####
print('Genes Res vs NoRes')

#### Calculate the observed test statistics
print('Calculate the observed test statistics for each replicate')
 
n = dim(arr_cts_rem)[1]
obs_sn_rm = matrix(nrow = n, ncol = 3)

for (i in 1:n){
	obs_sn_rm[i,] = get_stats(arr_cts_rem[i,,])
}

rownames(obs_sn_rm) = dimnames(arr_cts_rem)$replicate
obs_ge_rem = get_stats(cts_ge_rem)

#### Genes Plot Res vs NoRes
plt_all_rm = sens_fig_13_df(arr_cts_rem, 'Remission')
sens_plt_rm = sens_plot_fig_13(plt_all_rm, 'Remission')
sens_plt_rm
ggsave('../sensitivity/intermed/f1_sens_rm_plt.png', sens_plt_rm,
		height = 10, width = 10)

#### Permute
print('Permute')
 
# Iterate through the replicates to do the permutations
nperm = 2000
perm_lst_rm = rep_perms(engr_arr, arr_cts_rem, nperm, rs_pv_mgm, 'Remission')
save(perm_lst_rm, file = '../sensitivity/intermed/perm_lst_rm.RData')
perm_df_rm = sens_perm_df(perm_lst_rm[['stat_arr']])
perm_fx_rm = sens_perm_plot(perm_df_rm, 'fx')
perm_xfx_rm = sens_perm_plot(perm_df_rm, 'xfx')
perm_x2fx_rm = sens_perm_plot(perm_df_rm, 'x2fx')
 
perm_fx_rm
perm_xfx_rm
perm_x2fx_rm

ggsave('../sensitivity/intermed/perm_fx_rm.png', perm_fx_rm, 
       height = 10, width = 12)
ggsave('../sensitivity/intermed/perm_xfx_rm.png', perm_xfx_rm, 
       height = 10, width = 12)
ggsave('../sensitivity/intermed/perm_x2fx_rm.png', perm_x2fx_rm, 
       height = 10, width = 12)

# load('../sensitivity/intermed/perm_lst_rm.RData')
# Get the observed p-values

load('../permut_data/intermed/perms_ge_rem.RData')
obs_pvals_ge_rm = get_pvals(perms_ge_rem[['stat_mat']])

# Arrange the data for plotting

pval_arr_rm = perm_lst_rm[['pval_arr']]
colnames(pval_arr_rm) = names(obs_pvals_ge_rm)
head(pval_arr_rm)

reps = dimnames(perm_lst[[1]])$replicate
cutoffs = get_reps(reps)
pval_df_rm = (pval_arr_rm
           %>% data.frame()
           %>% mutate(Replicate = reps,
                      AbsCutoff = as.numeric(cutoffs[['abs']]),
                      PresCutoff = as.numeric(cutoffs[['pres']]))
           %>% pivot_longer(c(-Replicate, -AbsCutoff, -PresCutoff), 
                            names_to  = 'Statistic', 
                            values_to = 'P')
           %>% mutate(Observed = case_when(AbsCutoff == 25 & PresCutoff == 90 ~
                                               'Used',
                                           TRUE ~ 'Sens')))
head(pval_df_rm)

sens_pval_rm = ggplot(pval_df_rm, aes(P)) +
    geom_histogram(bins = 4) +
    facet_wrap(~Statistic) +
    geom_linerange(aes(ymin = 0, ymax = 50),
               data = filter(pval_df, Observed == 'Used'),
               colour = 'red') +
    xlab('Permutation p-value') +
    scale_x_continuous(labels = calc_lab)
sens_pval_rm

ggsave('../sensitivity/sens_pval_rm.png', sens_pval_rm)

# PICK UP HERE

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
