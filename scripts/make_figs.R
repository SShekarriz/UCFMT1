#!/usr/bin/env Rscript

##### Setup #####

# Import packages, set working directory, source functions
library(tidyverse)
theme_set(theme_bw())
# setwd('~/Projects/Favours/Sharok/UCFMT1/scripts/')
source('./permutation_functions.R')

# Create the p-value storage matrix

##### 16S Setup #####

### Import the data ###

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


### Get the engraftment matrix ###

# convert the 16s data to long format
long_mark_16s = mark_to_long(marker_lvl_16s, strains16s)
# Create a matrix of engraftment counts per patient
engr_16s = get_engraft(long_mark_16s, mapfile, cutoff_abs_16S, cutoff_pres_16S,
                       Treatment %in% c('FMT', 'Placebo'))


#### 16S FMT vs Placebo ####

# I think the only sensible way to plot this is to plot the full data and just
# talk about what we did for the statistics in the methods section.

#### Fig 1 
cts_16s_plt = count_engraft(engr_16s, tx_pv, 'Treatment')

f1_16s = plot_fig_13(cts_16s_plt, 'Treatment')
ggsave('../plots/f1_16s.pdf',f1_16s, height = 5, width = 5.5)
f1_16s

#### Fig 2 

f2_16s = plot_fig_2(cts_16s_plt)
ggsave('../plots/f2_16s.pdf',f2_16s, height = 5, width = 5.5)
f2_16s

#### 16S Plot Res vs. NoRes ####

# I think the only sensible way to plot this is to plot the full data and just
# talk about what we did for the statistics in the methods section.

#### Fig 1 
cts_16s_plt = count_engraft(engr_16s, rs_pv, 'Remission')

f3_16s = plot_fig_13(cts_16s_plt, 'Remission')
f3_16s
ggsave('../plots/f3_16s.pdf',f3_16s, height = 5, width = 5.5)

#### Setup Metagenomics ####

# This mapfile is used for every non-16S feature type
mapfile <- read.csv("../data/UCFMT1_METAmap.txt", sep = "\t")

# These patient vectors are used for every non-16S feature type
tx_pv = pat_to_vect(mapfile, 'Fig_lab', 'Treatment')
rs_pv = pat_to_vect(filter(mapfile, Treatment == 'FMT'), 'Fig_lab', 
                          'Remission')

#### Species Setup ####

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

#### Fig 1
f1_sp = plot_fig_13(cts_sp, 'Treatment')
f1_sp
ggsave('../plots/f1_sp.pdf',f1_sp, height = 5, width = 5.5)

#### Fig 2 

f2_sp = plot_fig_2(cts_sp)
f2_sp
ggsave('../plots/f2_sp.pdf',f2_sp, height = 5, width = 5.5)

#### Species Res vs NoRes ####

f3_sp = plot_fig_13(cts_sp_rem, 'Remission')
f3_sp
ggsave('../plots/f3_sp.pdf',f3_sp, height = 5, width = 5.5)

#### Strains Setup ####

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

f1_st = plot_fig_13(cts_st, 'Treatment')
f1_st
ggsave('../plots/f1_st.pdf',f1_st, height = 5, width = 5.5)

f2_st = plot_fig_2(cts_st)
f2_st
ggsave('../plots/f2_st.pdf',f2_st, height = 5, width = 5.5)

#### Strains Res vs NoRes ####

f3_st = plot_fig_13(cts_st_rem, 'Remission')
f3_st
ggsave('../plots/f3_st.pdf',f3_st, height = 5, width = 5.5)

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

f1_mg = plot_fig_13(cts_mg, 'Treatment')
f1_mg
ggsave('../plots/f1_mg.pdf',f1_mg, height = 5, width = 5.5)

f2_mg = plot_fig_2(cts_mg)
f2_mg
ggsave('../plots/f2_mg.pdf',f2_mg, height = 5, width = 5.5)

#### MAGs Res vs NoRes ####

f3_mg = plot_fig_13(cts_mg_rem,'Remission')
f3_mg
ggsave('../plots/f3_mg.pdf',f3_mg, height = 5, width = 5.5)

#### Genes Setup ####

### Import the data

marker_lvl_ge <- read.csv("../data/genes_lvl.csv")

# Set the presence and absence cutoffs for patients. These aren't used for donor
# B because all the genes came from the donor B assemblies.
cutoff_pres_ge = 90
cutoff_abs_ge = 25

### Get the engraftment matrix

long_mark_ge = mark_to_long(marker_lvl_ge)
engr_ge = get_engraft(long_mark_ge, mapfile, cutoff_abs_ge, cutoff_pres_ge,
                       Treatment %in% c('FMT', 'Placebo'))


#### Genes FMT vs Placebo ####

f1_ge = plot_fig_13(cts_ge, 'Treatment')
f1_ge
ggsave('../plots/f1_ge.pdf',f1_ge, height = 5, width = 5.5)

f2_ge = plot_fig_2(cts_ge)
f2_ge
ggsave('../plots/f2_ge.pdf',f2_ge, height = 5, width = 5.5)

#### Genes Res vs NoRes ####

f3_ge = plot_fig_13(cts_ge_rem, 'Remission')
f3_ge
ggsave('../plots/f3_ge.pdf',f3_ge, height = 5, width = 5.5)
