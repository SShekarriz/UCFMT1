# This script is to do all data wrangling, organization, etc. to create the data
# structures that will be used by run_permutation.R, make_paperFigs.R, and any
# other analysis or visualizing in the UCFMT1 metagenomics project.

library(tidyverse)

### Run Engraftment Analysis ####

source('./permutation_functions.R')

#### 16S ####

# Import 16s mapfile and create named vectors of patient IDs/groups for permuting
mapfile_16s <- read.csv("../data/UCFMT1_16Smap.txt")
tx_pv_16s = pat_to_vect(mapfile_16s, 'Fig_lab', 'Treatment')
rs_pv_16s = pat_to_vect(filter(mapfile_16s, Treatment == 'FMT'), 'Fig_lab', 
                          'Remission')

# Import the 16s ASV data
marker_lvl_16s <- read.csv("../data/16s_lvl_RA.csv")

# Import donor B 16s ASVs and create a vector
strains16s <- read.csv("../data/jcsz_16s_core_DonB.csv")
strains16s %>%
  pull() -> strains16s

# Set the cutoffs for 16s presence and absence
cutoff_pres_16S = 0.00013
cutoff_abs_16S = 0


### Get the engraftment matrix 

# convert the 16s data to long format
long_mark_16s = mark_to_long(marker_lvl_16s, strains16s)
# Create a matrix of engraftment counts per patient
engr_16s = get_engraft(long_mark_16s, mapfile_16s, cutoff_abs_16S, cutoff_pres_16S,
                       Treatment %in% c('FMT', 'Placebo'))

# Count engraftments per group
cts_16s_plt = count_engraft(engr_16s, tx_pv_16s, 'Treatment')
cts_16s_plt_rem = count_engraft(engr_16s, rs_pv_16s, 'Remission')

##### Save 16S Engraftment Objects ####
eng_16s = c('mapfile_16s', 'tx_pv_16s', 'rs_pv_16s', 'marker_lvl_16s',
            'strains16s', 'cutoff_pres_16S', 'cutoff_abs_16S', 'long_mark_16s',
            'engr_16s', 'cts_16s_plt', 'cts_16s_plt_rem')
save(list = eng_16s, file = '../processed_data/eng_16s.rds')

#### Setup Metagenomics ####

# This mapfile is used for every non-16S feature type
mapfile_mgm <- read.csv("../data/UCFMT1_METAmap.txt", sep = "\t")

# These patient vectors are used for every non-16S feature type
tx_pv_mgm = pat_to_vect(mapfile_mgm, 'Fig_lab', 'Treatment')
rs_pv_mgm = pat_to_vect(filter(mapfile_mgm, Treatment == 'FMT'), 'Fig_lab', 
                          'Remission')

##### Save Metagenomics Objects ####

mgm_setup = c('mapfile_mgm', 'tx_pv_mgm','rs_pv_mgm')
save(list = mgm_setup, file = '../processed_data/mgm_setup.rds')

#### Species ####

# Read in the species-level markers
marker_lvl_sp <- read.csv("../data/species_lvl_RA.csv")

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

##### Save Species Engraftment ####

eng_sp = c('marker_lvl_sp','donorB_species_all', 'donorB_species', 
           'cutoff_pres_sp', 'cutoff_abs_sp', 'strains_sp', 'long_mark_sp',
           'engr_sp', 'cts_sp', 'cts_sp_rem')
save(list = eng_sp, file = '../processed_data/eng_sp.rds')

#### Strains ####

# Import the strain-level marker data
marker_lvl_st <- read.csv("../data/marker_lvl_RA.csv")

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

##### Save Strain Engraftment ####
eng_st = c('marker_lvl_st', 'donorB_strains_all', 'donorB_strains', 
           'strains_st', 'cutoff_pres_st', 'cutoff_abs_st', 'long_mark_st', 
           'engr_st', 'cts_st', 'cts_st_rem')
save(list = eng_st, file = '../processed_data/eng_st.rds')

#### MAG ####

##### Read in the contig IDs and their bins ####
c_inBin <- read.csv("../data/Assembly/contigs_inBins.txt", sep = "\t", header = F)
c_inMAG <- read.csv("../data/Assembly/contigs_inMAGs.txt", sep = "\t", header = F)
colnames(c_inBin) <- c("contig", "bin_id")
colnames(c_inMAG) <- c("contig", "bin_id")

# NB that c_inMAG is a proper subset of c_inBin, so I'm not sure we need both.

# Get a vector of the unique MAG ids.
MAGs = (c_inMAG 
        %>% select(bin_id) 
        %>% distinct() 
        %>% pull())


##### Read in the MAG coverage ####
# define colnames for marker_lvl_mg
heads_cover <- c("sample.id",
                 "contig", "startpos", "endpos", 
                 "numreads", "covbases", "coverage", "meandepth",
                 "meanbaseq", "meanmapq",
                 "ReadSampling", "Quality")

path="../data/perfect_mapping/DBD"
patt="_perfect_.*.cover" # these are sample coverage files. Rows are contigs.
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
marker_lvl_mg_all = (marker_lvl_mg 
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
marker_lvl_mg = (marker_lvl_mg_all 
                 %>% filter(str_detect(sample, "PMCL")) 
                 %>% select(sample, bin_id, bin_coverage) 
                 %>% spread(sample, bin_coverage) 
                 %>% rename(Marker= bin_id))

##### Save marker_lvl_mg ####

save(marker_lvl_mg, file = '../processed_data/marker_lvl_mg.rds')

# Set the presence and absence cutoffs for patients. DonorB doesn't need these
# because all the MAGs come from donor B
cutoff_pres_mg = 75
cutoff_abs_mg = 20


### Get the engraftment matrix

long_mark_mg = mark_to_long(marker_lvl_mg)
engr_mg = get_engraft(long_mark_mg, mapfile_mgm, cutoff_abs_mg, cutoff_pres_mg,
                       Treatment %in% c('FMT', 'Placebo'))
cts_mg = count_engraft(engr_mg, tx_pv_mgm, 'Treatment')
cts_mg_rem = count_engraft(engr_mg, rs_pv_mgm, 'Remission')

##### Save MAG Engraftment ####

eng_mg = c('marker_lvl_mg_all', 'marker_lvl_mg', 'cutoff_pres_mg', 
           'cutoff_abs_mg', 'long_mark_mg', 'engr_mg', 'cts_mg', 'cts_mg_rem',
           'MAGs')
save(list = eng_mg, file = '../processed_data/eng_mg.rds')

#### Genes ####

### Import the data

marker_lvl_ge <- read.csv("../data/genes_lvl.csv")

# Set the presence and absence cutoffs for patients. These aren't used for donor
# B because all the genes came from the donor B assemblies.
cutoff_pres_ge = 90
cutoff_abs_ge = 25

### Get the engraftment matrix

long_mark_ge = mark_to_long(marker_lvl_ge)
engr_ge = get_engraft(long_mark_ge, mapfile_mgm, cutoff_abs_ge, cutoff_pres_ge,
                       Treatment %in% c('FMT', 'Placebo'))
cts_ge = count_engraft(engr_ge, tx_pv_mgm, 'Treatment')
cts_ge_rem = count_engraft(engr_ge, rs_pv_mgm, 'Remission')

##### Save Gene Engraftment ####
eng_ge = c('marker_lvl_ge', 'cutoff_pres_ge', 'cutoff_abs_ge', 'long_mark_ge',
           'engr_ge', 'cts_ge', 'cts_ge_rem')
save(list = eng_ge, file = '../processed_data/eng_ge.rds')

### Prep Benchmarking Data ####

#### 16S ####

##### Donor B Samples ####

# Using the same core DonorB ASVs as the engraftment uses: strains16s

donorB_16s_samples <- read.csv("../data/16s_lvl_RA_donorB_samples.csv")
donorB_16s_samples %>%
  gather(sample, abundance, -Marker) %>%
 # any marker with abundance > 0
# markers that were present in at least one donor B sample
  filter(!is.na(abundance)) %>%
  filter(as.numeric(abundance) > 0) -> donorB_16s_samples

##### Save 16S Benchmarking Data ####
bm_16s = c('strains16s', 'donorB_16s_samples')
save(list = bm_16s, file = '../processed_data/bm_16s.rds')

### Prep Profile Data ####

### Prep Visualization Data ####