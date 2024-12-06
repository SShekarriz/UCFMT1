library(tidyverse)
source('./permutation_functions.R')
source('./paperfig_functions.R')

# This mapfile is used for every non-16S feature type
mapfile_mgm <- read.csv("../data/UCFMT1_METAmap.txt", sep = "\t")

# These patient vectors are used for every non-16S feature type
tx_pv_mgm = pat_to_vect(mapfile_mgm, 'Fig_lab', 'Treatment')
rs_pv_mgm = pat_to_vect(filter(mapfile_mgm, Treatment == 'FMT'), 'Fig_lab', 
                        'Remission')
## Genes ####

### Run Engraftment Functions

### Import the data

marker_lvl_ge <- read.csv("../data/genes_lvl.csv")

# marker_lvl_tst = marker_lvl_ge[sample(1:nrow(marker_lvl_ge),
                                    # 10000),]

# Set the presence and absence cutoffs for patients. These aren't used for donor
# B because all the genes came from the donor B assemblies.
cutoff_abs_ge = seq(0, 40, 5)
cutoff_pres_ge = seq(70, 100, 5)

all_cutoff = expand.grid(cutoff_abs_ge, cutoff_pres_ge)
colnames(all_cutoff) = c('abs', 'pres')

### Get the engraftment matrix

print('Make long mark ge')
long_mark_ge = mark_to_long(marker_lvl_ge)

print('Make engr_ge')
engr_ge = get_engraft(long_mark_ge, mapfile_mgm, 25,
                      90, Treatment %in% c('FMT', 'Placebo'))
# engr_arr = array(dim = c(nrow(all_cutoff), 755662, 24),
#                  dimnames = list('replicate' = paste('abs', all_cutoff$abs,
#                                                      'pres', all_cutoff$pres,
#                                                      sep = '_'),
#                                  'genes' = rownames(engr_ge),
#                                  'samples' = colnames(engr_ge)))
# for (i in 1:nrow(all_cutoff)){
#     cp = all_cutoff$pres[i]
#     ca = all_cutoff$abs[i]
#     print(i)
#     engr_arr[i,,] = get_engraft(long_mark_ge, mapfile_mgm, ca, cp,
#                                Treatment %in% c('FMT', 'Placebo'))
# }

# save(engr_arr, file = '../sensitivity/intermed/engr_arr.RData')
print('load engr_arr')
load('../sensitivity/intermed/engr_arr.RData')
# 
cts_ge = count_engraft(engr_ge, tx_pv_mgm, 'Treatment')
cts_ge_rem = count_engraft(engr_ge, rs_pv_mgm, 'Remission')
# arr_cts = array(dim = c(nrow(all_cutoff), 755662, 2),
#                 dimnames = list('replicate' = dimnames(engr_arr)$replicate,
#                                 'genes' = dimnames(engr_arr)$genes,
#                                 'groups' = colnames(cts_ge)))
# arr_cts_rem = array(dim = c(nrow(all_cutoff), 755662, 2),
#                 dimnames = list('replicate' = dimnames(engr_arr)$replicate,
#                                 'genes' = dimnames(engr_arr)$genes,
#                                 'groups' = colnames(cts_ge_rem)))
# for (i in 1:nrow(all_cutoff)){
#     arr_cts[i,,] = count_engraft(engr_arr[i,,], tx_pv_mgm, 'Treatment')
#     arr_cts_rem[i,,] = count_engraft(engr_arr[i,,], rs_pv_mgm, 'Remission')
# }
# save(arr_cts, arr_cts_rem, file = '../sensitivity/intermed/cts_arrs.RData')
print('load cts_arr')
load('../sensitivity/intermed/cts_arrs.RData')

### Process Data for Visualization

# profile the feature types:
print('make profile ge')
profile_ge = get_profile_pat(long_mark_ge, mapfile_mgm, 25, 90,
                             Treatment %in% c('FMT', 'Placebo'))
# arr_prof = array(dim = c(nrow(all_cutoff), 755662, 24),
#                  dimnames = list('replicate' = dimnames(engr_arr)$replicate,
#                                  'genes' = rownames(profile_ge),
#                                  'samples' = colnames(profile_ge)))

# for (i in 1:nrow(all_cutoff)){
#     print(i)
#     tst = get_profile_pat(long_mark_ge, mapfile_mgm, 
#                                     all_cutoff$abs[i], all_cutoff$pres[i],
#                                     Treatment %in% c('FMT', 'Placebo'))
#     print('made the matrix')
#     arr_prof[i,,] = as.matrix(tst)
# }

# save(arr_prof, file = '../sensitivity/intermed/arr_prof.RData')
load('../sensitivity/intermed/arr_prof.RData')


profile_ge_don = get_profile_don(long_mark_ge, mapfile_mgm, 25, 90,
                                 Treatment %in% c('FMT','Placebo'))

# lst_prof_don = list()
# nms = paste('abs', all_cutoff$abs, 'pres', all_cutoff$pres, sep = '_')
# print('start the for loop')
# for (i in 1:length(nms)){
#     ca = all_cutoff$abs[i]
#     cp = all_cutoff$pres[i]
#     
#     print(i)
#     print(nms[i])
#     lst_prof_don[[nms[i]]] = get_profile_don(long_mark_ge, mapfile_mgm,
#                                              ca, cp, Treatment %in% c('FMT', 
#                                                                       'Placebo'))
# }
# 
# save(lst_prof_don, profile_ge_don, 
#      file = '../sensitivity/intermed/lst_prof_don.RData')
load('../sensitivity/intermed/lst_prof_don.RData')

# Wait, what is gene coverage actually doing anyway?

cov_hist = ggplot(long_mark_ge, aes(x = abundance)) +
    geom_histogram(bins = 100) +
    theme_bw()

cov_hist_log = cov_hist + scale_y_log10()
cov_hist_sqrt = cov_hist + scale_y_sqrt()
cov_hist
cov_hist_log
cov_hist_sqrt

ggsave('../sensitivity/gene_cov_hist.png', cov_hist)
ggsave('../sensitivity/gene_cov_hist_log.png', cov_hist_log)
ggsave('../sensitivity/gene_cov_hist_sqrt.png', cov_hist_sqrt)
