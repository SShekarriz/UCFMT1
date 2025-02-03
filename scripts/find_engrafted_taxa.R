#!/usr/bin/env Rscript

# Setup #####

## Identify features that are multiply engrafted in FMT. These are the
## interesting ones. 

# Import packages, set working directory, source functions
library(tidyverse)

# 16S #####

## Import the data ####
load('./processed_data/eng_16s.rds')

## Identify engrafted taxa ####
fmt_eng_16s = (cts_16s_plt
               %>% data.frame()
               %>% filter(tx_tot >= 3)
               %>% select(tx_tot)
               %>% rename(nPat = tx_tot)
               %>% rownames_to_column('TaxInfo'))

## Write Files ####
write.csv(fmt_eng_16s, file = './results/engrafted_16s_taxa.csv')

# Species ####
load('./processed_data/eng_sp.rds')

## Identify engrafted taxa ####
fmt_eng_sp = (cts_sp
               %>% data.frame()
               %>% filter(tx_tot >= 3)
               %>% select(tx_tot)
               %>% rename(nPat = tx_tot)
               %>% rownames_to_column('TaxInfo'))

## Write Files ####
write.csv(fmt_eng_sp, file = './results/engrafted_sp_taxa.csv')

# Strains ####
load('./processed_data/eng_st.rds')

## Identify engrafted taxa ####
fmt_eng_st = (cts_st
               %>% data.frame()
               %>% filter(tx_tot >= 3)
               %>% select(tx_tot)
               %>% rename(nPat = tx_tot)
               %>% rownames_to_column('TaxInfo'))

## Write Files ####
write.csv(fmt_eng_st, file = './results/engrafted_st_taxa.csv')

# MAGs ####
load('./processed_data/eng_mg.rds')

## Identify engrafted taxa ####
fmt_eng_mg = (cts_mg
               %>% data.frame()
               %>% filter(tx_tot >= 3)
               %>% select(tx_tot)
               %>% rename(nPat = tx_tot)
               %>% rownames_to_column('TaxInfo'))

## Write Files ####
write.csv(fmt_eng_mg, file = './results/engrafted_mg_taxa.csv')

# Genes ####
load('./processed_data/eng_ge.rds')

## Identify engrafted taxa ####
fmt_eng_ge = (cts_ge
               %>% data.frame()
               %>% filter(tx_tot >= 3)
               %>% select(tx_tot)
               %>% rename(nPat = tx_tot)
               %>% rownames_to_column('TaxInfo'))

## Write Files ####
write.csv(fmt_eng_ge, file = './results/engrafted_ge_taxa.csv')
