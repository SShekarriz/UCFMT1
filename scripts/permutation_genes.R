#!/usr/bin/env Rscript
#Shekarriz Sep11,2023

library(tidyverse)
library(foreach)
library(doParallel)

mapfile <- read.csv("../data/UCFMT1_METAmap.txt", sep = "\t")
marker_lvl_RA <- read.csv("../data/genes_lvl.csv")

cutoff = 90

## Functions
# Make the marker data frame long
mark_to_long = function(marker_lvl) {
  # This can be done once. It doesn't need to be repeated each time
  # convert marker-lvl to long format and make StudyID variable to match mapfile
  # and select only donor's strains
  marker_lvl %>%
    gather(sample, abundance, -Marker) %>%
    mutate(Study_ID = case_when(
      str_detect(sample,"PMCL") ~ paste(gsub("_.*", "",sample)),
      TRUE ~ paste(sample))) -> long_markerlvl
  return(long_markerlvl)
}

################################################################################
# a function to detect common engraftment based on strains (metaphlan markers)
################################################################################
get_engraft <- function(long_markerlvl, mapfile, cutoff=cutoff, ...){
  # complete list of engrafted
  engraft <- (long_markerlvl 
              %>% left_join(mapfile) 
              %>% filter(...)
              %>% spread(Timepoint, abundance, fill = 0) 
	      %>% filter(WK0 < 5 & WK6 >= cutoff) 
              %>% select(Marker, Study_ID, WK6) 
              %>% spread(Study_ID, WK6, fill = 0) 
              %>% gather(Study_ID, WK6, -Marker)
              %>% left_join(mapfile))
  return(engraft)
}

count_engraft <- function(engraft, group_cols){
  # Count how many markers engrafted in how many patients
  g1 = group_cols
  g2 = c('Marker', group_cols)
  ct_df = (engraft
           %>% group_by(.[,g2])
           %>% summarize(n = sum(WK6 > 0))
           %>% pivot_wider(names_from = {{ g1 }},
                           values_from = n)
           %>% ungroup()
           %>% mutate(Total = rowSums(across(2:ncol(.)))))
  # Uncomment the next bit if you want the cumulative sums, but it adds 0.8s
  # to the run time!
  # ct_df = (ct_df
  #          %>% mutate(n = rev(cumsum(rev(n))))) 
  return(ct_df)
  
}
################################################################################


## permutation of treatment FMT vs placebo (n=12)
# Number of permutations
num_permutations <- 10
it_res = list() # when not parallel
# Remove unnecessary columns from the mapfile to speed up 
# the loop (saves 10s per loop)
map_mod = (mapfile
           %>% select(Study_ID, Timepoint, Treatment, Fig_lab))
# Converting the markers to long format only needs to be done once, saves 5s per
# loop
markers_long = mark_to_long(marker_lvl_RA)

# Set up the cores for the parallelization and create the cluster
# cores = 6 # when testing parallel
# cores = 60 # when running parallel
# cl = makeCluster(cores) # when parallel

# Load required packages into the clusters
# clusterEvalQ(cl, library(tidyverse)) # when parallel

# Start parallel work
# registerDoParallel(cl) # when parallel
# it_res = foreach(i = 1:num_permutations) %dopar% { # when parallel
for (i in 1:num_permutations) {  # when not parallel
  # Permute the Treatment column
  permuted_data <- (map_mod
                    %>% select(Fig_lab, Treatment)
                    %>% unique()
                    %>% mutate(Treatment = sample(Treatment, replace = FALSE))
                    %>% right_join(select(map_mod, -Treatment))
                    %>% select(-Fig_lab)) # don't need this anymore
  tmplst = list()
  tmplst[['raw']] = permuted_data
  tmplst[['counts']] = (markers_long
                             %>% get_engraft(permuted_data, cutoff, 
                                             Treatment %in% c('FMT', 'Placebo'))
                             %>% count_engraft(group_cols = 'Treatment')
                             %>% mutate(Uniqueness = 
                                          case_when(FMT == 0 & Placebo == 0 ~ 'zero',
                                                    FMT == 0 ~ 'Placebo',
                                                    Placebo == 0 ~ 'FMT',
                                                    TRUE ~ 'Both'),
                                        it = paste('it', i, sep = '_')))
  tmplst[['Npts_tx']] = rbind(
    (tmplst[['counts']]
     %>% count(FMT)
     %>% rename(Npts = FMT)
     %>% mutate(Treatment = 'FMT')),
    (tmplst[['counts']]
     %>% count(Placebo)
     %>% rename(Npts = Placebo)
     %>% mutate(Treatment = 'Placebo')))
  tmplst[['Npts_total']] = (tmplst[['counts']]
                                 %>% count(Total, Uniqueness)
                                 %>% rename(Npts = Total)
                                 %>% mutate(Treatment = 'Total'))
  # tmplst  # when parallel
  it_res[[i]] = tmplst # when not parallel
}
# Stop parallel work
# stopCluster(cl)

save(it_res, file = '../permut_data/intermed/genes_it_res.RData')

### Permuted AUC
#approxfun creates an approximate function from the data using imputation
#integrate performs numerical integration to get the area under the curve
auc_null = matrix(NA, ncol = 2, nrow = num_permutations)
colnames(auc_null) = c('deltaAUC','iteration')

for (i in 1:num_permutations){
  fmt_fn = with(filter(it_res[[i]][['Npts_tx']], Treatment == 'FMT'),
                approxfun(Npts, n))
  plc_fn = with(filter(it_res[[i]][['Npts_tx']], Treatment == 'Placebo'),
                approxfun(Npts, n))
  fmt_auc = with(filter(it_res[[i]][['Npts_tx']], Treatment == 'FMT'),
                 integrate(fmt_fn, 1, max(Npts))$value)
  plc_auc = with(filter(it_res[[i]][['Npts_tx']], Treatment == 'Placebo'),
                 integrate(plc_fn, 1, max(Npts))$value)
  d_auc = fmt_auc - plc_auc 
  auc_null[i,] = c(d_auc, i)
  
}
write.table(auc_null, "../permut_data/genes_auc_null.txt",
            sep = "\t", row.names = F, quote = F)


## permutations of responders vs non-responders within FMT
# Number of permutations
num_permutations <- 10

# fmt_res = list()

# Remove unnecessary columns from the mapfile to speed up the loop (saves 10s per 
# loop)
map_mod = (mapfile
           %>% filter(Treatment == 'FMT')
           %>% select(Study_ID, Timepoint, Remission, Fig_lab))
# Loop through permutations
cores = 6
cl = makeCluster(cores)
clusterEvalQ(cl, library(tidyverse))
registerDoParallel(cl)
fmt_res = foreach(i = 1:num_permutations) %dopar% {
# for (i in 1:num_permutations) {
  # Permute the Treatment column
  permuted_data <- (map_mod
                    %>% select(Fig_lab, Remission)
                    %>% unique()
                    %>% mutate(Remission = sample(Remission, replace = FALSE))
                    %>% right_join(select(map_mod, -Remission))
                    %>% select(-Fig_lab)) # don't need this anymore
  tmplst = list()
  tmplst[['raw']] = permuted_data
  tmplst[['counts']] = (markers_long
                              %>% get_engraft(permuted_data, cutoff, 
                                              Remission %in% c('Res', 'NoRes'))
                              %>% count_engraft(group_cols = 'Remission')
                              %>% mutate(Uniqueness = 
                                           case_when(Res == 0 & NoRes == 0 ~ 'zero',
                                                     Res == 0 ~ 'NoRes',
                                                     NoRes == 0 ~ 'Res',
                                                     TRUE ~ 'Both'),
                                         it = paste('it', i, sep = '_')))
  tmplst[['Npts_tx']] = rbind(
    (tmplst[['counts']]
     %>% count(Res)
     %>% rename(Npts = Res)
     %>% mutate(Remission = 'Res')),
    (tmplst[['counts']]
     %>% count(NoRes)
     %>% rename(Npts = NoRes)
     %>% mutate(Remission = 'NoRes')))
  tmplst[['Npts_total']] = (tmplst[['counts']]
                                 %>% count(Total, Uniqueness)
                                 %>% rename(Npts = Total)
                                 %>% mutate(Remission = 'Total'))
  tmplst
  
}
stopCluster(cl)

save(tmplst, file = '../permute_data/intermed/gene_it_res.RData')

### Permuted AUC
#approxfun creates an approximate function from the data using imputation
#integrate performs numerical integration to get the area under the curve

rm_auc_null = matrix(NA, ncol = 2, nrow = num_permutations)
colnames(rm_auc_null) = c('deltaAUC','iteration')

for (i in 1:num_permutations){
  res_fn = with(filter(fmt_res[[i]][['Npts_tx']], Remission == 'Res'),
                approxfun(Npts, n))
  nrs_fn = with(filter(fmt_res[[i]][['Npts_tx']], Remission == 'NoRes'),
                approxfun(Npts, n))
  res_auc = with(filter(fmt_res[[i]][['Npts_tx']], Remission == 'Res'),
                 sum(integrate(res_fn, 0, 1)$value,
                     integrate(res_fn, 1, max(Npts))$value))
  nrs_auc = with(filter(fmt_res[[i]][['Npts_tx']], Remission == 'NoRes'),
                 sum(integrate(nrs_fn, 0, 1)$value,
                     integrate(nrs_fn, 1, max(Npts))$value))
  d_auc = res_auc - nrs_auc
  rm_auc_null[i,] = c(d_auc, i)
  
}

write.table(rm_auc_null, "../permut_data/genes_remission_auc_null.txt",
            sep = "\t", row.names = F, quote = F)
