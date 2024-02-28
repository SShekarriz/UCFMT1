library(tidyverse)

#### Start by making sure the processed data are identical

#### Process the data Sharok's way ####

mapfile <- read.csv("../data/UCFMT1_METAmap.txt", sep = "\t")

org <- read.csv("../data/Metaphlan4/merged_metaphlan4.txt", 
                sep = "\t", header = F, 
                comment.char = "#")
org %>%
  filter(V1 != "clade_name") %>%
  # selecting only at strain level- contains all lineage- avoid duplicates
  filter(grepl("t__S", V1)) -> metaphlan
colnames(metaphlan) <- gsub(".mp.profile_sp", "", org[1,])

metaphlan %>%
  mutate(temp=clade_name) %>%
  separate(temp, c("Kingdom", "Phylum", "Class", 
                       "Order", "Family", "Genus", 
                       "Species", "marker"), sep = "[|]") %>%
  select(clade_name,Kingdom:marker) %>%
  mutate(sp_marker= gsub("t__", "", marker)) %>%
  mutate(sp_marker=gsub("_group", "", sp_marker))-> lineage

metaphlan %>%
  gather(sample, abundance, -clade_name) %>%
  #filtering donor B markers
  filter(str_detect(sample, "DonorB_D_")) %>%
  #any marker with abundance > 0
  #markers that were present in at least one donor B sample
  filter(!is.na(abundance)) %>%
  filter(as.numeric(abundance) > 0) %>%
  mutate(sample= gsub("DonorB_D_", "", sample)) -> donorB_species

mapfile %>% select(Study_ID, Treatment, Timepoint) -> map

### seting engraftment cut-off
# based on donor B benchmarking data: 0.005 (as shown above)
# minimum relative abundance detected maximum (top quartile) of species
cutoff= 0.005
####

# working with only donor B species that meeet cutoff
# > average of minimum detections
donorB_species %>%
  filter(as.numeric(abundance) > cutoff) -> donorB_species


donorB_species %>%
  select(clade_name) %>% distinct() %>%
  pull() -> B_species

metaphlan %>%
  rename(Marker = "clade_name") -> marker_lvl_RA

# This is added by Jake vvv
marker_lvl_RA = (marker_lvl_RA
                 %>% mutate(across(-Marker, as.numeric)))
# ^^^ this is added by Jake

# Make the marker data frame long
mark_to_long_ss = function(marker_lvl, strains) {
    # This can be done once. It doesn't need to be repeated each time
    # convert marker-lvl to long format and make StudyID variable to match mapfile
    # and select only donor's strains
    marker_lvl %>%
        filter(Marker %in% strains) %>%
        gather(sample, abundance, -Marker) %>%
        mutate(Study_ID = case_when(
            str_detect(sample,"PMCL") ~ paste(gsub("_.*", "",sample)),
            TRUE ~ paste(sample))) -> long_markerlvl
    return(long_markerlvl)
}


markers_long = mark_to_long_ss(marker_lvl_RA, B_species)

#### Rename Sharok's data objects ####

marker_lvl_ss = marker_lvl_RA
dB_ss = B_species
mapfile_ss = mapfile
map_ss = map
long_mark_ss = markers_long

#### Process the data Jake's way ####

mapfile <- read.csv("../data/UCFMT1_METAmap.txt", sep = "\t")
# These patient vectors are used for every non-16S feature type
pat_to_vect = function(mapfile, patcol, catcol){
    # Create a two-column data frame of unique patient IDs and the grouping column
    # (probably Treatment or Remission)
    to_perm = mapfile[,c(patcol,catcol)]
    to_perm = unique(to_perm)
    
    # Create a vector of the grouping variable and name it with the patient IDs
    # to make it easy to permute
    perm = to_perm[,catcol]
    names(perm) = to_perm$Fig_lab
    
    # Remove any NAs
    perm = na.omit(perm)
    return(perm)
}

tx_pv = pat_to_vect(mapfile, 'Fig_lab', 'Treatment')
rs_pv = pat_to_vect(filter(mapfile, Treatment == 'FMT'), 'Fig_lab', 
                          'Remission')


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

mark_to_long_jc = function(marker_lvl, donB = NULL) {
    
    # For some feature types there is not a resticted list of donor B features
    # because all features are from donor B
    if(!is.null(donB)){
        # But if there is such a list, remove any features that are absent from
        # it
        marker_lvl = filter(marker_lvl, Marker %in% donB)
    } 
    
    # Gather the marker data frame to a long format and fix the naming problems
    # in the sample names
    long_markerlvl = (marker_lvl
                      %>% gather(sample_orig, abundance, -Marker)
                 %>% mutate(sample = case_when(str_detect(sample_orig,"PMCL") ~ 
                                             paste(gsub("_.*", "",sample_orig)),
                                                    TRUE ~ paste(sample_orig))))
    return(long_markerlvl)
}


long_mark_sp = mark_to_long_jc(marker_lvl_sp, strains_sp)

#### Rename Jake's variables ####
marker_lvl_jc = marker_lvl_sp
dB_jc = strains_sp
mapfile_jc = mapfile
long_mark_jc = long_mark_sp

#### Are these identical? ####

# Start with marker_lvl

all.equal(marker_lvl_jc, marker_lvl_ss)
# For some reason, Sharok's data are character instead of numeric. Fix that and
# carry on

marker_lvl_ss = (marker_lvl_ss
                 %>% mutate(across(-Marker, as.numeric)))
all.equal(marker_lvl_jc, marker_lvl_ss)

# Okay, these become the same once that problem is solved. I'll go back and fix
# it in the data processing.

all.equal(long_mark_ss, long_mark_jc)

# The names are different. That's okay. 
colnames(long_mark_ss)
colnames(long_mark_jc)

long_jc_tst = long_mark_jc
colnames(long_jc_tst) = colnames(long_mark_ss)
all.equal(long_jc_tst, long_mark_ss)

# Okay. Except for the column names, the long_mark data frames are now identical.

# Let's take it to the plotting data frames.

#### Sharok's plotting data frame ####

get_engraft_ss <- function(long_markerlvl, mapfile, cutoff=cutoff, ...){
  # complete list of engrafted
  engraft <- (long_markerlvl 
              %>% left_join(mapfile) 
              %>% filter(...)
              %>% spread(Timepoint, abundance, fill = 0) 
              %>% filter(WK0 <= 0 & WK6 > cutoff) 
              %>% select(Marker, Study_ID, WK6) 
              %>% spread(Study_ID, WK6, fill = 0) 
              %>% gather(Study_ID, WK6, -Marker)
              %>% left_join(mapfile))
    return(engraft)
}

count_engraft_ss <- function(engraft, group_cols){
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

# Remove unnecessary columns from the mapfile to speed up the loop (saves 10s per 
# loop)
map_mod_ss = (mapfile_ss
           %>% select(Study_ID, Timepoint, Treatment, Fig_lab))

cts_ss = (long_mark_ss
       %>% get_engraft_ss(map_mod_ss, cutoff, 
                       Treatment %in% c('FMT', 'Placebo'))
       %>% count_engraft_ss(group_cols = 'Treatment')
       %>% mutate(Uniqueness = 
                    case_when(FMT == 0 & Placebo == 0 ~ 'zero',
                              FMT == 0 ~ 'Placebo',
                              Placebo == 0 ~ 'FMT',
                              TRUE ~ 'Both'),
                  it = 'observed'))
npts_tx_ss = rbind(
  (cts_ss
   %>% count(FMT)
   %>% rename(Npts = FMT)
   %>% mutate(Treatment = 'FMT')),
  (cts_ss
   %>% count(Placebo)
   %>% rename(Npts = Placebo)
   %>% mutate(Treatment = 'Placebo')))

plt_df_ss = filter(npts_tx_ss, Npts > 0)

f1_ss <- ggplot(plt_df_ss, aes(Npts, n, color=Treatment)) +
  geom_point() + geom_line() +
  # geom_ribbon(data = ribdat, aes(ymin=Placebo,ymax=FMT, x = Npts),
  #             linetype=2, alpha=0.1, inherit.aes = FALSE)+
  # annotate("text", x=9, y=90, size=2,
  #         label= paste("Permutation test,\n p=", auc_pval, sep = " ")) +
  scale_x_continuous(breaks = seq(0, max(plt_df_ss$Npts), 1)) +
  scale_color_manual(values = cols) +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_blank())
f1_ss

#### Jake's plotting data frame ####

get_engraft_jc = function(long_markerlvl, mapfile, cutoff_abs, cutoff_pres,  ...){
    engraft = (long_markerlvl
               # Join the feature count/coverage data with the metadata and
               # filter to the desired groups
               %>% left_join(mapfile, by = c('sample' = 'Study_ID'))
               %>% filter(...)
               # Reshape the table so we can determine what has been engrafted
               %>% select(Marker, Timepoint, abundance, Fig_lab)
               %>% spread(Timepoint, abundance, fill = 0)
               # Set engrafted column to 1 if engraftment criteria are met, 0
               # otherwise
               %>% mutate(engrafted = case_when(WK0 <= cutoff_abs & WK6 > cutoff_pres ~ 1,
                                                TRUE ~ 0))
               # Select the columns and reshape the data frame to create a wide
               # matrix of 1/0 engraftment values
               %>% select(Marker, Fig_lab,  engrafted)
               %>% pivot_wider(names_from = Fig_lab, values_from = engrafted,
                               values_fill = 0)
               # Make the Marker column the rownames
               %>% column_to_rownames('Marker'))
    
    # Turn it into a matrix and remove all rows that are 0 everywhere. This
    # removal has no effect on the test statistic, it's just tidier
    engraft = as.matrix(engraft)
    engraft = engraft[rowSums(engraft) > 0,]
    return(engraft)
}

count_engraft_jc = function(engraft, pat_vect, txrm){
    # Is this treatment or remission?
    if (txrm == 'Treatment'){
        bl = 'Placebo'
        tx = 'FMT'
    } else if (txrm == 'Remission'){
        bl = 'NoRes'
        tx = 'Res'
    }
    
    # Subset the engraftment matrix to just the columns identified as FMT or Res
    # (see above) in the patient vector and take the row sums to get the number
    # of times each feature was engrafted in the FMT (Res) group
    engr_tx = engraft[,names(pat_vect[pat_vect == tx])]
    tx_tot = rowSums(engr_tx)
    
    # As above, but with the baseline (Placebo/NoRes) group
    engr_bl = engraft[,names(pat_vect[pat_vect == bl])]
    bl_tot = rowSums(engr_bl)
    
    # create a data frame with the two sets of counts as columns
    engr_ct = cbind(tx_tot, bl_tot)
    return(engr_ct)
}

cutoff_pres_sp = 0.005
cutoff_abs_sp = 0
engr_jc = get_engraft_jc(long_mark_jc, mapfile_jc, cutoff_abs_sp, cutoff_pres_sp,
                       Treatment %in% c('FMT', 'Placebo'))


# Get the engraftment counts and calculate the test statistics
cts_jc = count_engraft_jc(engr_jc, tx_pv, 'Treatment')

#### Compare the counts ####

all.equal(cts_jc, cts_ss)
jc_cts_tst = data.frame(cts_jc)
ss_cts_tst = data.frame(cts_ss)
summary(jc_cts_tst)
summary(ss_cts_tst)

# reformat Sharok's to match Jake's
ss_cts_tst = data.frame(tx_tot = cts_ss$FMT,
                        bl_tot = cts_ss$Placebo)
rownames(ss_cts_tst) = cts_ss$Marker

all.equal(jc_cts_tst, ss_cts_tst)

# Okay, the counting is different. In particular, there are _more_ matches in
# Sharok's data than in mine. Let's check how that happened.

# Sharok's data had > for Wk6, mind had >=. I fixed that but it didn't help.

sum(rownames(ss_cts_tst) %in% rownames(jc_cts_tst))

# all 310 of the ones I catch are in Sharok's. But he catches 67 more.

missed = rownames(ss_cts_tst)[!(rownames(ss_cts_tst) %in% rownames(jc_cts_tst))]
length(missed)

missed[1:5]

head(long_mark_jc)

tst = long_mark_jc %>% filter(Marker == missed[1], !startsWith(sample, 'Don'))
tst

# This is from the old get_engraft function used in Sharok's scripts. I'm going
# through it line by line to see what it does differently from the new version.


  engraft2 <- (long_mark_ss
              %>% left_join(mapfile_ss)
              %>% filter(Treatment %in% c('FMT', 'Placebo'))
              %>% spread(Timepoint, abundance, fill = 0))
              # %>% filter(WK0 <= 0 & WK6 > cutoff)
              # %>% select(Marker, Study_ID, WK6) 
              # %>% spread(Study_ID, WK6, fill = 0) 
              # %>% gather(Study_ID, WK6, -Marker)
              # %>% left_join(mapfile))


# I don't understand what this broken spread is doing. Let's find out

tst_df = (long_mark_ss
          %>% left_join(mapfile_ss)
          %>% filter(sample %in% c('PMCL293_S293', 'PMCL297_S294'))
          %>% filter(Marker %in% long_mark_ss$Marker[1:4])
          %>% mutate(abundance = c(1,2,3,4,5,6,7,8)))
tst_df

tst_df_wide = (tst_df
               %>% spread(Timepoint, abundance, fill = 0))
filter(tst_df, Marker == tst_df$Marker[1])
filter(tst_df_wide, Marker == tst_df$Marker[1])

# Okay. I see what is happening. The spread in the old version of the
# get_engraft function doesn't put WK6 and WK0 on the same line. Instead, it
# keeps two lines for that marker/patient combo. In the first line, it keeps the
# WK0 value and puts 0 in WK6. In the second line it keeps the WK6 value and
# puts 0 in WK 0. So then when we do the filter, what it does is end up
# returning every marker where the WK6 value exceeds the WK6 cutoff, regardless
# of what that marker's value was at WK0.  


#### Testing it out ####

# I'll test this. If I get the same number of engrafted markers by simply
# filtering on WK6 > cutoff, I'm right.

get_engraft_tst <- function(long_markerlvl, mapfile, cutoff=cutoff, ...){
  # complete list of engrafted
  engraft <- (long_markerlvl 
              %>% left_join(mapfile) 
              %>% filter(...)
              %>% filter(Timepoint == 'WK6')
              %>% rename(WK6 = abundance)
              %>% filter(WK6 > cutoff) 
              %>% select(Marker, Study_ID, WK6) 
              %>% spread(Study_ID, WK6, fill = 0) 
              %>% gather(Study_ID, WK6, -Marker)
              %>% left_join(mapfile))
    return(engraft)
}

cts_tst = (long_mark_ss
       %>% get_engraft_tst(map_mod_ss, cutoff, 
                       Treatment %in% c('FMT', 'Placebo'))
       %>% count_engraft_ss(group_cols = 'Treatment')
       %>% mutate(Uniqueness = 
                    case_when(FMT == 0 & Placebo == 0 ~ 'zero',
                              FMT == 0 ~ 'Placebo',
                              Placebo == 0 ~ 'FMT',
                              TRUE ~ 'Both'),
                  it = 'observed'))
dim(cts_tst)
all.equal(cts_tst, cts_ss)

# Yep. We were not filtering on the WK0 value, even though we meant to.

# So, the upshot is that the old graphs (and p-values) are wrong. Instead of
# engrafted values, we had counts of just things that were high at week 6. So we
# should re-generate our plots and re-calculate our p-values using the new
# functions, which don't have this error.

