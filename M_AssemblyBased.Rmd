---
title: "Metabat_bin_analysis"
author: "Sharok"
date: "7/26/2019"
output: html_document
---

These are my analysis to detect MAG and MABs in CEMG and DMG data. Assign
taxonomic assignment to them and produce map files for further visulaization
in anvio platform.


```{r echo=FALSE, message=FALSE, warning=FALSE}
library(ggplot2)
library(tidyverse)
library(tidytext)
#library(taxize)#get taxonomy lineage
library(edgeR)
#library(Biostrings) #reading fasta file
library(seqinr) #reading fasta file


colours1 <- c("#02b6d3","#668400","#f1bc8c","#e58eff","#76b000","#ffa4e4","#bf4200","#c0019a","#c452d7","#feb873","#922d4f",
             "#763cc2","#9a3cc2","#ec0078","#c8cb82","#bc6c00","#ff7d99","#ff7bb4","#cb88a5","#63c9ff","#af003a","#b77668",
             "#9fb500","#006aaf","#ff5245","#c00067","#5fdabb","#f7bc5b","#913412","#ffa477","#009a2a","#01b97b","#de0042",
             "#109eff","#f07312","#816190","#f11067","#a11b2d","#835d00","#8857e0","#fd59d2","#a80225","#91d86a","#019d6d",
             "#eb2ea5","#62409a","#285d35","#007b19","#31e27c","#723c89","#006826","#ff8845","#ffab9a","#a2aeff","#f45a2b",
             "#ff3d54","#7b4429","#0294e5","#4343b6","#1b4eaa","#96d68b","#ff7369","#cb2d10","#74da9d","#dd7300","#843c41",
             "#81411c","#ff68b4","#4fbf36","#f3b3de","#c3cd65","#736a00","#ff73ef","#ff9e37","#8f1396","#ff9fab","#8486ff",
             "#9d175e","#957a00","#385e00","#978b56","#025eb5","#9c3100","#00714c","#5173fe","#008a77","#cc9000","#d0a400",
             "#b70072","#d10090","#86336c","#02a4ae","#63dd8f","#e1c28e","#8d3729","#1b5588","#524b85","#00b8ff","#cea7ff","#ff6e74")
colours = c('grey69',"#4f8579","#783fcc","#69d24d","#cb4bbd","#c6dc46",
					"#542871","#78d792","#cc4472","#83d7d0","#d44d33","#676fcd",
					"#ceb854","#403d57","#b97839","#84a4cb","#588038","#c68ac4",
					"#48472a","#c9c39c","#6e2b34","#c78889")

BinStats <- read.csv("Metabat_ChechM_edits.txt", sep = "\t", 
                     stringsAsFactors = FALSE)
Bin_species <- read.csv("Metabat_bracken_species.txt", sep = "\t",
                        stringsAsFactors = FALSE)
Bin_genus <- read.csv("Metabat_bracken_genus.txt", sep = "\t",
                      stringsAsFactors = FALSE)

##########################################################
## A function to make nested strip in ggplot
## https://github.com/teunbrand/ggnomics/blob/standalone-facet_nested/R/
source('~/Drive2/software/ggplot_faceting/facet_nested.R')
##########################################################

```

# Quality of bins
Defining MAG, and MABs based on completeness and contamination values reported
by CheckM and vosialize the qulaity and quantity of these bins

```{r echo=FALSE}
BinStats %>% filter(as.numeric(Completeness) >= 70 &
                      as.numeric(Contamination) <= 10) %>%
  mutate(Quality = rep("MAG")) -> MAGs

BinStats %>% filter(between(as.numeric(Completeness), 50, 70) &
                      as.numeric(Contamination) <= 10) %>%
  mutate(Quality = rep("MAB")) -> MABs

BinStats %>% filter(as.numeric(Completeness) < 50 |
                      as.numeric(Contamination) > 10) %>%
  mutate(Quality = rep("Bin")) -> Bins

GoodColumns <- c('Bin.Id', 'Marker.lineage', 'Completeness', 'Contamination', 
                 'Strain.heterogeneity', 'Quality', 'BinName') 


bind_rows(MAGs, MABs, Bins) %>%
  mutate(BinName = str_replace(Bin.Id, "bin.", "bin_")) %>%
  select(GoodColumns)-> BinsFine


BinsFine %>% select(BinName, Contamination, Strain.heterogeneity, 
                    Marker.lineage, Quality) -> CS_df
BinsFine %>% select(BinName, Completeness) -> C_df

CS_df %>%
  gather(AxisX, ValueX, -BinName, -Marker.lineage, -Quality) %>%
  left_join(C_df) %>%
  gather(AxisY, ValueY, -AxisX, -ValueX, -BinName, 
         -Marker.lineage, -Quality) -> BinFineLong
  

#########################################################################
#
#########################################################################

ContigInBins <- read.csv("contigs_inBins.txt", sep = "\t", 
                     stringsAsFactors = FALSE, header = FALSE)
colnames(ContigInBins) <- c("Contig", "BinName")

# selecting contigs in MAGs and MAB
ContigInBins %>%
  left_join(BinsFine, by = "BinName") %>%
  filter(Quality %in% c("MAG", "MAB")) %>%
  select(Contig, BinName) -> ContigsIn_MAG_MAB
write.table(ContigsIn_MAG_MAB, "contigs_in_MAG_MAB.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# source file to descibe source of MAGs and MABs

Source_color <- data.frame(Source = c("Common", "CEMG+", "DMG+", 
                                      "Shared", "OtherSamples"),
                           Color = c("#FF0000", "#009900", "#004C99", 
                                     "#4C0099", "#CC6600"))


```

# Taxonomic annotation of bins
Here are my work on Bracken and krakenHLL outputs. Bracken did the assignment
only at genus level which means we've got "unknown" information for cases
where we couldn't assign proper genus information. So I used krakenHll tool
and parsed the output of all krakenHLL files together using my own script to 
identify best hit by comparing multiple taxonomic level. It turns out we couldn't
assign phylum level information to some of the bins even using krakenHLL. So I 
used ChecKM lineage marker information for those bins that couldn't get any
phylum information. In summary, I tried hard to assign at least phylum information
to all of these bins by any mean but I'll generated another list of bins with
roboust taxonomomic information as well that can be used for the purpose of 
phylogeny. 

```{r echo=FALSE, message=FALSE, warning=FALSE}
###############################################################################
##################TAXONOMIC INFORMATION FROM KRAKENHLL#########################
# this is the NCBI full lineage information
full_lineage <- read.csv("~/Drive2/software/lineages-2019-02-20.csv",
                         stringsAsFactors = FALSE)
colnames(full_lineage)[1] <- "TaxID"
full_lineage %>% 
  select(TaxID, superkingdom, phylum, class, 
         order, family, genus, species) -> full_lineage

# Using BRACKEN Species and genus information:> Bracken output is fairly
# accurate for speceis/genus level data but not lower so for each Bin the
# Bracken output must be checken for species/genus level
Bin_species_brac <- read.csv("Metabat_bracken_species.txt", sep = "\t",
                        stringsAsFactors = FALSE)
Bin_genus_brac <- read.csv("Metabat_bracken_genus.txt", sep = "\t",
                      stringsAsFactors = FALSE)
#select only those taxonomy that have more than 55 % of bin assigned to st:
Bin_species_brac %>% filter(fraction_total_reads >= 0.55) %>%
  mutate(Brac_species = name) %>%
  select(BinName, Brac_species) -> Bin_species_brac
Bin_genus_brac %>% filter(fraction_total_reads >= 0.55) %>%
  mutate(Brac_genus = name) %>%
  left_join(Bin_species_brac, by = "BinName") %>%
  select(BinName, Brac_genus, Brac_species)-> Bin_bracken
  

# Importing krakenHll output using my own script to select best hit
Bin_krakenHll <- read.csv("All_bins_krakenHLL", sep = "\t",
                      stringsAsFactors = FALSE, header = FALSE,
                      na.strings = c("","NA"))
col <- c("BinName", "Percen", "Reads", "TaxRead", "kmers", 
         "Dup", "Cov", "TaxID", "Rank", "TaxName")
colnames(Bin_krakenHll) <- col


Bin_krakenHll %>%
  mutate(BinName = 
gsub("/home/sharok/Drive2/DonorB_metagenomics/DonorB_DMG_CEMG_bin_assignment/krakenHLL//", 
     "", BinName)) %>%
  mutate(BinName = gsub(".fa.krakenHLL", "", BinName)) %>%
  mutate(BinName = gsub("bin.", "bin_", BinName)) %>%
  mutate(BinName = trimws(BinName)) %>%
  replace_na(list(TaxName = "unknown",
                  Rank = "unknown")) %>%
  #removing the Uknown bins
  filter(!TaxName == "unknown") %>%
  select(BinName, Rank, TaxID, TaxName, Percen) %>%
  #removing [] from taxa name
  mutate(TaxName = gsub("[[]", "", TaxName)) %>%
  mutate(TaxName = gsub("[]]", "", TaxName)) %>%
  left_join(full_lineage, by = "TaxID") -> df_krakenHll


BinsFine %>%
  select(BinName, Marker.lineage) %>%
  separate(Marker.lineage, c("Rank", "Name"), sep = "__") %>%
  filter(!is.na(Name)) %>%
  separate(Name, c("TaxName", "UID"), sep = " ") -> CheckM_taxa
  CheckM_taxa %>%
  distinct(TaxName) -> CheckM_UID 
# Manually find the checkM taxIDs so that it can be merged with krakenHLL
# and can be used in case of not finding proper taxa in krakenHll output
c("2", "1224", "28221", "186802", "28890", "186826", "976", "909929",
"186803", "171549", "1300", "31953", "201174", "80864", "186801", "543", 
"1301", "712", "838", "29547") -> CheckM_UID$TaxID
CheckM_UID %>% 
  mutate(TaxID = as.numeric(TaxID)) %>%
  right_join(CheckM_taxa, by = "TaxName") %>%
  left_join(full_lineage, by = "TaxID") %>%
  select(BinName, TaxName, superkingdom, phylum, class, order, 
         family, genus) -> CheckM_complete
colnames(CheckM_complete) <- c("BinName","checkm_TaxName",
                               "checkm_superkingdom",
                               "checkm_phylum","checkm_class",
                               "checkm_order", "checkm_family", 
                               "checkm_genus")

#########################################
### FINAL TAXONOMY FOR BIN: this script combine krakenHLL, checkM, and Bracken
### taxa assignment and generate best consensous information from these sources
#########################################
 #Merging the lineage info with Bin quality file
 BinsFine %>%
  #ADD the checkm taxonomic information 372bins
  left_join(CheckM_complete, by = "BinName") %>%
  # ADD krakenHLL taxa info 370bins
  left_join(df_krakenHll, by = "BinName") %>%
  # ADD bracken taxa assignment 254bins 
  left_join(Bin_bracken, by = "BinName") %>%
  # for those with N
  #A assignment add checkm marker lineage
  mutate(phylum = case_when(is.na(phylum) ~ as.character(checkm_phylum),
                            phylum == "" ~ as.character(checkm_phylum),
                                TRUE ~ as.character(phylum))) %>%
   # get the checkM superkingdom for those that don't even have phylum
  mutate(phylum = case_when(is.na(phylum) ~ paste("k__", checkm_superkingdom,
                                                  sep = ""),
                            phylum == "" ~ paste("k__", checkm_superkingdom,
                                                  sep = ""),
                                TRUE ~ as.character(phylum))) %>%
   # if order is empty get it from checkM order
   mutate(order = case_when(is.na(order) ~ as.character(checkm_order),
                            order == "" ~ as.character(checkm_order),
                                TRUE ~ as.character(order))) %>%
   #get the phylum information if it's NA or empty
   mutate(order = case_when(is.na(order) ~ paste("p__", phylum,
                                                  sep = ""),
                            order == "" ~ paste("p__", phylum,
                                                  sep = ""),
                                TRUE ~ as.character(order))) %>%
   mutate(order = gsub("p__k__", "k__", order)) %>%
   # if family is empty get it from checkM family first
   mutate(family = case_when(is.na(family) ~ as.character(checkm_family),
                            family == "" ~ as.character(checkm_family),
                                TRUE ~ as.character(family))) %>%
   #get the order information if it's NA or empty
   mutate(family = case_when(is.na(family) ~ paste("o__", order,
                                                  sep = ""),
                            family == "" ~ paste("o__", order,
                                                  sep = ""),
                                TRUE ~ as.character(family))) %>%
   mutate(family = gsub("o__p__", "p__", family)) %>%
   mutate(family = gsub("o__k__", "k__", family)) %>%
   # if genus is empty get it from bracken genus level
   mutate(genus = case_when(is.na(genus) ~ as.character(Brac_genus),
                            genus == "" ~ as.character(Brac_genus),
                                TRUE ~ as.character(genus))) %>%
   # if genus is empty get it from family
    mutate(genus = case_when(is.na(genus) ~ paste("f__", family, sep = ""),
                           genus == "" ~  paste("f__", family, sep = ""),
                           TRUE ~ as.character(genus))) %>%
    mutate(genus = gsub("f__p__", "p__", genus)) %>%
    mutate(genus = gsub("f__k__", "k__", genus)) %>%
    mutate(genus = gsub("f__o__", "o__", genus)) %>%
   
  # if species is empty get it from bracken species level
   mutate(species = case_when(is.na(species) ~ as.character(Brac_species),
                            species == "" ~ as.character(Brac_species),
                                TRUE ~ as.character(species))) %>% 
  # if species is empty get it from genus column 
    mutate(species = case_when(is.na(species) ~ paste("g__", genus, sep = ""),
                             species == "" ~  paste("g__", genus, sep = ""),
                           TRUE ~ as.character(species))) %>%
    mutate(species = gsub("g__p__", "p__", species)) %>%
    mutate(species = gsub("g__k__", "k__", species)) %>%
    mutate(species = gsub("g__o__", "o__", species)) %>%
    mutate(species = gsub("g__f__", "f__", species)) %>%
   select(-c("checkm_TaxName", "checkm_superkingdom",
              "checkm_phylum","checkm_class","checkm_order", 
              "checkm_family", "checkm_genus", "class",
              "Brac_genus", "Brac_species")) -> BinsFine_taxa
 

 ####################################################
 # MAG_MAB_Complete
 ####################################################
 # This is a table contain complete information about MAG_MAB
 # in contains complete checkM and taxonomy, and length of MAGs etc
 
 # importing anivo summary of bins here
path1 <- paste("~/Drive2/UC_FMT/anvio_setup_DMG-CEMG-all/SAMPLES-MERGED",
"_UCpatients_2500_SUMMARY_METABAT_BINS/bins_summary.txt",
sep = "")
BINS_length <- read.csv(path1, sep = "\t")
colnames(BINS_length)[1] <- "BinName"

BinsFine_taxa %>%
   left_join(BINS_length, by = "BinName") -> BinsFine_complete
 write.table(BinsFine_complete, 
            "METABAT_447BINS_complete_information.txt", 
             sep = "\t", 
             row.names = FALSE, col.names = FALSE, quote = FALSE)

BinsFine_taxa %>%
   filter(Quality %in% c("MAG", "MAB")) %>%
   left_join(BINS_length, by = "BinName") -> MAG_MAB_complete


#######################################################
#### MAKING METABAT_MAG_MAB_phylum colours
Phylum_color <- data.frame(phylum = BinsFine_taxa %>%
                             filter(Quality %in% c("MAG", "MAB")) %>%
                             distinct(phylum),
                           Color = c("#FF0000", "#009900", "#004C99", 
                                     "#4C0099", "#CC6600",
                                     "#FF33FF", "#FFFF33"))
# selecting only MAG_MABs
 BinsFine_taxa %>%
   filter(Quality %in% c("MAG", "MAB")) %>%
   left_join(Phylum_color, by = "phylum") %>%
   select(BinName, phylum, Color) -> MAG_MAB_phylum_info
 write.table(MAG_MAB_phylum_info, 
            "contigs_in_MAG_MAB_PhylumInfo.txt", 
             sep = "\t", 
             row.names = FALSE, col.names = FALSE, quote = FALSE)

 
#### MAKING METABAT_MAG_MAB_order colours
Order_color <- data.frame(phylum = BinsFine_taxa %>%
                             filter(Quality %in% c("MAG", "MAB")) %>%
                             distinct(order),
                           Color = colours1[1:21])

 # selecting only MAG_MABs 
 BinsFine_taxa %>%
   filter(Quality %in% c("MAG", "MAB")) %>%
   left_join(Order_color, by = "order") %>%
   select(BinName, order, Color) -> MAG_MAB_order_info
  write.table(MAG_MAB_order_info, 
            "contigs_in_MAG_MAB_OrderInfo.txt", 
             sep = "\t", 
             row.names = FALSE, col.names = FALSE, quote = FALSE)


#removing those that couldn't assign any taxa via krakenHLL OR checkM
# I'm removing these because wanna make a clean phylogeny tree and it's gonna 
# be weired if we wont be able to assign even a phylum to a bacteria.
 
############################################################################ 
# selecting only clean MAG_MABs with proper taxa at phylum level
 BinsFine_taxa %>%
   filter(Quality %in% c("MAG", "MAB")) %>%
   left_join(Phylum_color, by = "phylum") %>%
   select(BinName, phylum, Color) %>%
   filter(!phylum == "k__Bacteria")-> MAG_MAB_clean_phylum_info
  write.table(MAG_MAB_clean_phylum_info, 
            "contigs_in_MAG_MAB_clean_PhylumInfo.txt",
             sep = "\t", 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
#new bin collection for clean bins (those assigned at least to a phylum)
  ContigInBins %>% 
  right_join(MAG_MAB_clean_phylum_info, by = "BinName") %>% 
  select(Contig, BinName) -> ContigsIn_MAG_MAB_clean
write.table(ContigsIn_MAG_MAB_clean, 
            "contigs_in_MAG_MAB_clean.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


 # selecting only clean MAG_MABs with proper taxa at order level
 BinsFine_taxa %>%
   filter(Quality %in% c("MAG", "MAB")) %>%
   left_join(Order_color, by = "order") %>%
   select(BinName, order, Color) %>%
   filter(!order %in% c("k__Bacteria"))-> MAG_MAB_clean_order_info
 write.table(MAG_MAB_clean_order_info, 
            "contigs_in_MAG_MAB_clean_OrderInfo.txt",
             sep = "\t", 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
 
##############################################################################
 # selecting only clean MAG_MABs with proper taxa at order level
 BinsFine_taxa %>%
   filter(Quality %in% c("MAG", "MAB")) %>%
   left_join(Order_color, by = "order") %>%
   select(BinName, order, Color) %>%
   filter(!order %in% c("k__Bacteria")) %>%
   filter(!str_detect(order, "p__*"))-> MAG_MAB_Superclean_order_info
 write.table(MAG_MAB_Superclean_order_info, 
            "contigs_in_MAG_MAB_Superclean_OrderInfo.txt",
             sep = "\t", 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
 
 #new bin collection for super clean bins (those assigned at least to an order)
  ContigInBins %>% 
  right_join(MAG_MAB_Superclean_order_info, by = "BinName") %>% 
  select(Contig, BinName) -> ContigsIn_MAG_MAB_Superclean
write.table(ContigsIn_MAG_MAB_Superclean, 
            "contigs_in_MAG_MAB_Superclean.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

```


OK, so to visualize tree of DonorB refined genomes I think the below anvio 
profile using METABAT_MAG_MAB_clean at phylum level can be used:
anvi-interactive -c contigs-fixed.db -p SAMPLES-MERGED_DonorB/PROFILE.db 
-C METABAT_MAG_MAB_clean 
-t METABAT_INFORMATION/Campbel_et_al_Ribosomalgenes_METABAT_MAG_MAB_clean_
protein_tree.txt --server-only


#Engraftment
Here are my codes to combine the bin information with the Patient's mapfile
and identify those bins that were not present before FMT and become abundant
after FMT. Those bins will be labelled as engrafted category of bins. Also I'll 
select a list of bins that were present before and after FMT and be labelled
as no change category

```{r}
#Import the mapfile from the UCFMT1
mapfile <- read.csv("MapFile_metagenomics_July17.txt", sep = "\t",
                    stringsAsFactors = FALSE) 
mapfile %>% filter(!Status == "Healthy") %>%
  mutate(Fig_lab = paste("pt", Trial_no, sep = "")) -> mapfile

path2 <- paste("~/Drive2/UC_FMT/anvio_setup_DMG-CEMG-all/SAMPLES-MERGED_UCpatients",
"_2500_newSamples_SUMMARY_MAG_MAB_genecoverage/bins_across_samples/detection.txt",
sep = "")

#Import the detection file across bins
detection <- read.csv(path2, sep = "\t", stringsAsFactors = FALSE)
colnames(detection) <- gsub("_MAP_SORTD", "", colnames(detection))

path3 <- paste("~/Drive2/UC_FMT/anvio_setup_DMG-CEMG-all/SAMPLES-MERGED_UCpatients",
"_2500_newSamples_SUMMARY_MAG_MAB_genecoverage/bins_across_samples/mean_coverage.txt",
sep = "")

################################################################################
#Import gtdbtk assignment on Bins
################################################################################
gtdb <- read.csv("gtdbtk.bac120.summary.tsv", sep = "\t")
gtdb %>% separate(classification, c("Kingdom", "Phylum", "Class",
                                   "Order", "Family", "Genus",
                                   "Species"), ";.__") %>%
  select(user_genome, Kingdom:Species) -> gtdb

detection %>%
  mutate(user_genome= gsub("bin_", "bin.", bins))%>%
  select(user_genome, bins)%>%
  left_join(gtdb, by="user_genome") %>%
  left_join(BinsFine_taxa %>%
            mutate(user_genome=Bin.Id) %>%
            select(user_genome, family, genus)%>%
            rename(Krak_fam=family, Krak_gen=genus),
            by = "user_genome") -> gtdb

################################################################################
# Done
################################################################################

#Import the coverage file across bins, make it long and add mapfile
coverage <- read.csv(path3, sep = "\t", stringsAsFactors = FALSE)
colnames(coverage) <- gsub("_MAP_SORTD", "", colnames(coverage))
coverage %>%
  gather(SampleName, coverage_value, -bins) %>%
  left_join(mapfile, by = "SampleName") %>%
  #removing duplicate samples (those that their depth increased)
  filter(!SampleName %in% c("s721_CCGTCC", "s797_GTGGCC")) %>%
  select(c("bins", "coverage_value", "Timepoint", "Initial")) %>%
  #adding "Time" term so it wont become inteager
  mutate(Timepoint = paste("Time", Timepoint, sep = "")) -> long_coverage
colnames(long_coverage)[1] <- "BinName"


path3_2 <- paste("~/Drive2/UC_FMT/anvio_setup_DMG-CEMG-all/SAMPLES-MERGED_UCpatients",
"_2500_newSamples_SUMMARY_MAG_MAB_genecoverage/bins_across_samples/variability.txt",
sep = "")

#Import the coverage file across bins, make it long and add mapfile
variability <- read.csv(path3_2, sep = "\t", stringsAsFactors = FALSE)
colnames(variability) <- gsub("_MAP_SORTD", "", colnames(variability))
variability %>%
  gather(SampleName, variability_value, -bins) %>%
  left_join(mapfile, by = "SampleName") %>%
  #removing duplicate samples (those that their depth increased)
  filter(!SampleName %in% c("s721_CCGTCC", "s797_GTGGCC")) %>%
  select(c("bins", "variability_value", "Timepoint", "Initial")) %>% 
  #adding "Time" term so it wont become inteager
  mutate(Timepoint = paste("Time", Timepoint, sep = "")) -> long_variability
colnames(long_variability)[1] <- "BinName"



Patient_fate <- data.frame(Initial = c("PM","AF","DAB","KG","WAF","BN",
                                       "TJD","PT","MP","EB","ES"),
                           Fate = c("R", "R", "P", "R", "N", "R", "N",
                                    "R", "N", "N", "R"))

# converting to long format so that can be merged with mapfile
# clean up and convert back to wide for engraftment detection
 detection %>% 
  gather(SampleName, detect_value, -bins) %>%
  left_join(mapfile, by = "SampleName") %>%
  #removing duplicate samples (those that their depth increased)
  filter(!SampleName %in% c("s721_CCGTCC", "s797_GTGGCC")) %>%
  select(c("bins", "detect_value", "Timepoint", "Initial", "Fig_lab")) %>%
  #adding "Time" term so it wont become inteager
  mutate(Timepoint = paste("Time", Timepoint, sep = "")) %>%
  # make it wide to easily select egrafted bins
  spread(key=Timepoint, value=detect_value) %>%
  group_by(Initial) %>%
  # define the type of MAG in terms of their change over FMT
  mutate(model= case_when(Time0 <= 0.1 & Time6 >= 0.6 ~ "Engraft",
                          Time0 >= 0.6 & Time6 <= 0.1 ~ "Lost",
                          Time0 >= 0.6 & Time6 >= 0.6 ~ "Present",
                          Time0 <= 0.6 & Time6 <= 0.6 ~ "NotPresent",
                            TRUE ~ "MinorChange")) %>%
  # make it long again
  gather(Timepoint, detect_value, -bins, -Initial, -model, -Fig_lab) %>% 
  # adding patient fate  
  left_join(y=Patient_fate, by = "Initial") %>%
  mutate(BinName = bins) %>% select(-bins) %>%
  #add coverage information to the final table
  left_join(long_coverage, by = c("BinName", "Timepoint", "Initial")) %>%
  #add variability information to the final table
  left_join(long_variability, by = c("BinName", "Timepoint", "Initial")) %>% 
  # add krakenHLL taxonomic information 
  left_join(BinsFine_taxa, by = "BinName") -> MAG_MAB_type_FMT
 
 #### COMBINE SNV with previous model, generate DS_model and add to the data>
 MAG_MAB_type_FMT %>%
 left_join(
 MAG_MAB_type_FMT %>%
   select(BinName, Timepoint, Initial, model, variability_value) %>%
   spread(key=Timepoint, value=variability_value) %>%
   group_by(Initial) %>%
   mutate(SNVmodel= case_when(Time6 <= 0.5 & Time0 >=1 ~ "StrainReplace",
                            TRUE ~ "PatientStrain")) %>%
   mutate(DS_model= case_when(model == "Present" | model == "MinorChange" ~ SNVmodel,
                              TRUE ~ model)) %>%
   gather(Timepoint, variability_value, -BinName, -Initial, -model,
          -SNVmodel, -DS_model) %>% 
   select(BinName, Timepoint, Initial, DS_model), 
   by = c("BinName", "Timepoint", "Initial")) -> MAG_MAB_type_FMT
   
MAG_MAB_type_FMT$Timepoint <- factor(MAG_MAB_type_FMT$Timepoint,
                                     levels = c("Time6", "Time0"))

###############################################################################
# New dataset the contsins all levels of data together.
###############################################################################

d_v_cols <-  c(        "Detect_0-0.2" = "#ffffcc",
                       "Detect_0.2-0.4" = "#c2e699",
                       "Detect_0.4-0.6" = "#78c679",
                       "Detect_0.6-0.8" = "#31a354", 
                       "Detect_0.8-1" = "#006837",
                      ########## SNV #################
                     "variability_LowDetect" = "NA",
                     "variability_0-0.25" = "#d73027",
                     "variability_0.25-0.5" = "#f46d43",
                     "variability_0.5-1" = "#fdae61",
                     "variability_1-10"  =  "#abd9e9", 
                     "variability_SNV>10" = "#2c7bb6")


# creating a bin ordering file.
MAG_MAB_type_FMT %>%
  mutate(BinName= gsub("bin_", "M", BinName)) %>%
  select(BinName, model) %>%
  group_by(BinName, model) %>% tally() %>% arrange(model) %>%
  spread(model, n, fill=0) %>%
  mutate(model = case_when(Engraft >= 2  ~ "Engraft",
                           Engraft == 0 & Lost >= 2  ~ "Lost",
                           TRUE ~ "Nochange")) %>%
  arrange(model) %>%
  select(BinName, model) -> MAG_MAB_items_order


MAG_MAB_type_FMT$Timepoint <- factor(MAG_MAB_type_FMT$Timepoint,
                                     levels = c("Time6", "Time0"))

MAG_MAB_type_FMT %>%
  # selecting data columns
  select(Initial, Fig_lab, model, DS_model, Timepoint, Fate, BinName,
         detect_value, variability_value) %>%
  mutate(BinName= gsub("bin_", "M", BinName)) %>%
  # Detection characters
mutate(D=case_when(between(detect_value, 0, 0.2) ~ paste("0-0.2"),
                              between(detect_value, 0.2, 0.4) ~ paste("0.2-0.4"),
                              between(detect_value, 0.4, 0.6) ~ paste("0.4-0.6"),
                              between(detect_value, 0.6, 0.8) ~ paste("0.6-0.8"),
                              between(detect_value, 0.8, 1) ~ paste("0.8-1"),
                              TRUE ~ as.character(detect_value))) %>% 
  # add label to characters
  mutate(D = paste("Detect", D, sep = "_")) %>%
  # Variability characters
mutate(S= case_when(  detect_value <= 0.2 ~ paste("LowDetect"),
                      between(variability_value, 0, 0.25) ~ paste("0-0.25"),
                      between(variability_value, 0.25, 0.5) ~ paste("0.25-0.5"),
                      between(variability_value, 0.5, 1) ~ paste("0.5-1"),
                      between(variability_value, 1, 10) ~ paste("1-10"),
                              variability_value >= 10 ~ paste("SNV>10"),
                              TRUE ~ as.character(variability_value))) %>%
  # add labels to characters
  mutate(S=paste("variability", S, sep = "_")) %>%
  # NOW removing the actual detection and variability values
  # in order to combine the these two character data together
  select(-detect_value, -variability_value) %>%
  # now combine
  gather(value_type, value,
         -Initial, -Fig_lab, -model, -DS_model, 
         -Timepoint, -Fate, -BinName) -> MAG_MAB_type_FMT_viz


MAG_MAB_type_FMT_viz$model <- factor(MAG_MAB_type_FMT_viz$model, 
                                     levels = c("Present", "MinorChange",
                                            "Lost", "Engraft",
                                            "NotPresent"))
MAG_MAB_type_FMT_viz$DS_model <- factor(MAG_MAB_type_FMT_viz$DS_model, 
                                    levels = c("PatientStrain", "StrainReplace",
                                            "Lost", "Engraft",
                                            "NotPresent"))

MAG_MAB_type_FMT_viz$BinName <- factor(MAG_MAB_type_FMT_viz$BinName,
                       levels = MAG_MAB_items_order %>% pull(BinName))

```



```{r}
ggplot(MAG_MAB_type_FMT_viz, 
       aes(BinName, Timepoint)) +
   geom_tile(aes(fill = value), alpha = 8/10) +
   scale_fill_manual(values= d_v_cols) +
   facet_grid(Fate+Fig_lab+value_type~., scale = "free", space = "free", switch = c("y")) +
   theme_bw()+ 
   #scale_x_discrete(position = "top")+
   theme(text = element_text(size=8),
        legend.position="bottom",
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 0),
        strip.text.y = element_text(angle = 0))
ggsave("METABAT_bin_analysis_figures/MAGMAB_D_S_all.png", 
       heigh = 18, width = 19, units = "cm")
```


```{r}
ggplot(MAG_MAB_type_FMT_viz, 
       aes(BinName, Timepoint)) +
   geom_tile(aes(fill = value), alpha = 8/10) +
   scale_fill_manual(values= d_v_cols) +
   facet_grid(Fate+Fig_lab+value_type~model, scale = "free", space = "free", switch = c("y")) +
   theme_bw() + 
   #scale_x_discrete(position = "top") +
   theme(text = element_text(size=10),
        legend.position="bottom",
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(angle = 0),
        strip.text.y = element_text(angle = 0))


```

```{r}
ggplot(MAG_MAB_type_FMT_viz %>% 
      filter(model!= "NotPresent"), 
       aes(BinName, Timepoint)) +
   geom_tile(aes(fill = value), alpha = 8/10) +
   scale_fill_manual(values= d_v_cols) +
   facet_grid(Fate+Fig_lab+value_type~., scale = "free", space = "free", switch = c("y")) +
   theme_bw()+ 
   #scale_x_discrete(position = "top")+
   theme(text = element_text(size=10),
        legend.position="bottom",
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(angle = 0),
        strip.text.y = element_text(angle = 0))
ggsave("METABAT_bin_analysis_figures/MAGMAB_D_S_all_3.png", 
    heigh = 18, width = 25, units = "cm")

```

```{r}

ggplot(MAG_MAB_type_FMT_viz %>% 
       filter(model!= "NotPresent"), 
       aes(BinName, Timepoint)) +
   geom_tile(aes(fill = value), alpha = 8/10) +
   scale_fill_manual(values= d_v_cols) +
   facet_grid(Fate+Fig_lab+value_type~model, scale = "free", space = "free", switch = c("y")) +
   theme_bw() + 
   #scale_x_discrete(position = "top") +
   theme(text = element_text(size=10),
        legend.position="bottom",
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text( angle = 90),
        strip.text.x = element_text(angle = 0),
        strip.text.y = element_text( angle = 0))



```


```{r}
ggplot(MAG_MAB_type_FMT_viz %>% 
      filter(model!= "NotPresent"), 
       aes(BinName, Timepoint)) +
   geom_tile(aes(fill = value), alpha = 8/10) +
   scale_fill_manual(values= d_v_cols) +
   facet_grid(Fate+Fig_lab+value_type~DS_model, scale = "free", space = "free", switch = c("y")) +
   theme_bw() + 
   #scale_x_discrete(position = "top") +
   theme(text = element_text(size=15),
        legend.position="bottom",
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90),
        strip.text.x = element_text(size = 25, angle = 0),
        strip.text.y = element_text(size = 15, angle = 0))

```


```{r}

################################################################################
# adding gtdb taxonomoy to the engrafted MAGs taxonomy
################################################################################

  

top_fam <- c("Lachnospiraceae", "Bacteroidaceae", "Oscillospiraceae", 
             "Ruminococcaceae","Acutalibacteraceae")

MAG_MAB_type_FMT_viz %>% 
            filter(DS_model %in% c("Engraft", "StrainReplace")) %>%
            left_join(gtdb %>%
                      mutate(BinName= gsub("bin.", "M", user_genome)) %>%
                      select(BinName, Order, Family, Genus, Krak_fam, Krak_gen),
                      by = "BinName") %>%
            # the 3 archea genomes were not classified with gtdb database
            # just replace those NAs with their the kraken assignment
            mutate(Family = case_when(is.na(Family) ~ paste(Krak_fam),
                                      TRUE ~ as.character(Family))) %>%
            # do the same for genus as well
            mutate(Genus = case_when(is.na(Genus) ~ paste(Krak_gen),
                                      TRUE ~ as.character(Genus))) %>%
            # making a new column for family labs: only showing the 
            # abundant familes and name the rest others:
            mutate(Family_lab = case_when(Family %in% top_fam ~ paste(Family),
                                          TRUE ~ paste("Other"))) -> tbl
  
  
 tbl$Fig_lab <- factor(tbl$Fig_lab,
                       levels=c("pt74", "pt25", "pt84", "pt85",
                                "pt79","pt75", "pt56", "pt80", "pt60",
                                "pt10", "pt4"))
  
DS_model_colors2 <- c("Engraft" = "#006d2c",
                     "StrainReplace" = "#e41a1c",
                     "PatientStrain" = "#377eb8",
                      "Lost" = "#ff7f00",
                      "NotPresent" = "#bdbdbd")

tbl$Family_lab <- factor(tbl$Family_lab, 
                         levels = rev(c("Lachnospiraceae", "Bacteroidaceae", 
                                    "Oscillospiraceae", 
             "Ruminococcaceae","Acutalibacteraceae", "Other")))

tbl %>% filter(Timepoint=="Time6") %>%
group_by(BinName) %>% tally() %>% arrange(n) %>% pull(BinName) -> BinName_ord

tbl$BinName <- factor(tbl$BinName, levels = BinName_ord)
  
ggplot(tbl %>% filter(Timepoint=="Time6"), aes(BinName,Fig_lab)) +
       geom_tile(aes(fill=DS_model)) +
       facet_grid(Fig_lab+Fate~Family_lab, scales = "free", 
                  space = "free", switch = c("y")) +
       scale_fill_manual(values = DS_model_colors2) +
       scale_x_reordered() +
       scale_y_reordered() +
       theme_bw() +
       theme(text = element_text(size=7),
        legend.position="none",
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, size=3),
        strip.text.x = element_blank(),
        #strip.text.x = element_text(angle = 90),
        strip.text.y.left = element_text(angle = 0))
ggsave("METABAT_bin_analysis_figures/MAGMAB_Main_2A.png", 
heigh = 8, width = 14, units = "cm")


MAG_MAB_type_FMT_viz %>% 
            ungroup() %>%
            filter(DS_model %in% c("Engraft", "StrainReplace")) %>%
            select(Fig_lab, BinName, DS_model) %>%
            distinct() %>%
            group_by(Fig_lab, DS_model) %>% tally() -> tbl_bar

tbl_bar$Fig_lab <- factor(tbl_bar$Fig_lab,
                      levels=rev(c("pt74", "pt25", "pt84", "pt85",
                                   "pt79","pt75", "pt56", "pt80", "pt60",
                                   "pt10", "pt4")))

ggplot(tbl_bar, aes(n, Fig_lab, fill = DS_model)) +
       geom_bar(stat = "identity") +
       scale_fill_manual(values = DS_model_colors2) +
       theme_classic() +
       scale_x_reverse() +
       scale_y_discrete(position = "right")  +
       theme(text = element_text(size=10),
             legend.title = element_blank(),
             legend.position=c(0.4,0.8),
        axis.title = element_blank(), 
        axis.text.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 0))
ggsave("METABAT_bin_analysis_figures/MAGMAB_Main_2B.png",
       width = 4, height = 8, units = "cm")



```




```{r}


MAG_MAB_type_FMT_viz %>%
left_join(gtdb %>%
mutate(BinName= gsub("bin.", "M", user_genome)) %>%
select(BinName, Order, Family, Genus, Krak_fam, Krak_gen),
by = "BinName") %>%
mutate(Family = case_when(is.na(Family) ~ paste(Krak_fam),
                                      TRUE ~ as.character(Family))) %>%
            # do the same for genus as well
            mutate(Genus = case_when(is.na(Genus) ~ paste(Krak_gen),
                                      TRUE ~ as.character(Genus))) %>%
            # making a new column for family labs: only showing the 
            # abundant familes and name the rest others:
            mutate(Family_lab = case_when(Family %in% top_fam ~ paste(Family),
                                          TRUE ~ paste("Other"))) -> tbl

  gtdb_col<- c("Lachnospiraceae" = "#f0027f", 
             "Bacteroidaceae" = "#386cb0", 
             "Oscillospiraceae" = "#ffff99",
             "Ruminococcaceae" = "#7fc97f",
             "Acutalibacteraceae" = "#bf5b17",
             "Other" = "#a6cee3")

tbl$Family_lab <- factor(tbl$Family_lab, 
                         levels = rev(c("Lachnospiraceae", "Bacteroidaceae", 
                                    "Oscillospiraceae", 
             "Ruminococcaceae","Acutalibacteraceae", "Other")))


tbl %>%
  ungroup() %>%
  select(BinName, Family_lab) %>%
  distinct() %>% group_by(Family_lab) %>% tally() %>% 
  mutate(Donor= "DonorB") %>%
  ggplot(aes(Donor,n, fill=Family_lab)) +
  geom_bar(stat = "identity") +
  coord_polar("y") +
  scale_fill_manual(values = gtdb_col) +
  theme_classic() +
  geom_text(aes(label=n), position = position_stack(vjust = 0.5),
            size = 3) +
  theme(text = element_text(size = 5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())

tbl %>%
  ungroup() %>%
  select(BinName, Family_lab) %>%
  distinct() %>% group_by(Family_lab) %>% tally() %>% 
  mutate(Donor= "DonorB") %>%
  ggplot(aes(Donor,n, fill=Family_lab)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = gtdb_col) +
  scale_y_continuous(breaks  = c(0, 100, 200, 255)) +
  theme_classic() +
  geom_text(aes(label=n), position = position_stack(vjust = 0.5),
            size = 3) +
  theme(text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank())
ggsave("METABAT_bin_analysis_figures/DonorB_pie.png",
       width = 4, height = 8, units = "cm")


tbl$Fig_lab <- factor(tbl$Fig_lab,
                      levels=rev(c("pt74", "pt25", "pt84", "pt85",
                                   "pt79","pt75", "pt56", "pt80", "pt60",
                                   "pt10", "pt4")))

R=c("pt74", "pt85", "pt78", "pt56", "pt79", "pt10", "pt4")
PN=c("pt84", "pt75", "pt80", "pt60", "pt25")

ggplot(tbl %>% 
            ungroup() %>%
            select(Fig_lab, BinName, DS_model, Fate) %>%
            distinct() %>%
            group_by(Fig_lab, DS_model, Fate) %>% tally() %>% ungroup() %>%
            mutate(group= case_when(Fig_lab %in% R ~ paste("R"),
                                    Fig_lab %in% PN ~ paste("PN"))),
aes(n, Fig_lab, fill = factor(DS_model, 
       levels=rev(c("NotPresent","PatientStrain",
                    "Engraft", "StrainReplace", "Lost"))))) +
       geom_bar(stat = "identity") +
       facet_grid(group~., scales = "free", space = "free") +
       scale_x_continuous(breaks  = c(0, 100, 200, 255),
                          limits = c(0,255)) +
       scale_fill_manual(values = DS_model_colors2) +
       theme_classic() +
  theme(text = element_text(size=10),
        legend.position="none",
        legend.title = element_blank(),
        axis.title = element_blank(), 
        strip.text= element_blank())
ggsave("METABAT_bin_analysis_figures/DonorB_Main1.png",
       width = 6, height = 8, units = "cm")

tbl %>% 
ungroup() %>%
select(Fig_lab, BinName, DS_model, Fate) %>%
distinct() %>%
group_by(Fig_lab, DS_model, Fate) %>% tally() %>% ungroup() -> report

report %>%
group_by(DS_model) %>% summarise(n=sum(n)) %>%
mutate(Avg= n/11,
       Perent= (n/11) / 255 * 100)

```




```{r}

###############################################################################
# Comparing Mayo score with Engraftment/Replacement in each patient
###############################################################################

mapfile %>%
  filter(!X.SampleID %in% c("FMT721-b", "FMT797-b")) %>%
  select(Fig_lab, Timepoint, Mayoscore) %>%
  spread(Timepoint, Mayoscore) %>%
  mutate(Imp = `0` - `6`) %>%
  left_join(mapfile %>%
           filter(!X.SampleID %in% c("FMT721-b", "FMT797-b")) %>%
           mutate(Remission = case_when(Remission == 0 ~ paste("N"),
                               Remission == 1 ~ paste("R"))) %>%
          filter(Timepoint == 6) %>%
          select(Fig_lab, Remission), by = "Fig_lab") -> MayoScoreImprovement

MAG_MAB_type_FMT_viz %>% 
  filter(DS_model %in% c("Engraft", "StrainReplace")) %>%
  filter(Timepoint == "Time6") %>%
  ungroup() %>%
  select(Fig_lab, DS_model, BinName) %>%
  distinct() %>%
  group_by(Fig_lab, DS_model) %>% summarise(counts= n()) %>%
  left_join(MayoScoreImprovement, by = "Fig_lab") -> Mayo_vs_strains

ggplot(Mayo_vs_strains, aes(x = Imp, y = counts)) +
  geom_point() +
  stat_smooth() +
  #facet_grid(~DS_model) +
  theme_bw()

ggplot(Mayo_vs_strains, aes(x = Imp, y = counts)) +
  geom_point() +
  stat_smooth() +
  facet_grid(~DS_model) +
  theme_bw()


```


```{r}
##############################################################################
#################################
# Finding common engrafted Bins:
#################################
# FINDING BINS THAT ENGRAFTED AT LEAST FOR 3 PATIENTS
MAG_MAB_type_FMT_viz %>%
#filter(DS_model %in% c("Engraft")) %>%
filter(DS_model %in% c("Engraft", "StrainReplace")) %>%
filter(Timepoint == "Time6") %>%
select(Fig_lab, DS_model, Fate, BinName) %>% distinct() %>%
group_by(BinName) %>% tally() %>% arrange(n) -> common_eng

MAGengraftment_curv <- data.frame(Npts = c("1", "2", "3", "4", "5", "6"),
    eMAGs = c(length(common_eng %>% filter(n >= 1) %>% pull(BinName)),
              length(common_eng %>% filter(n >= 2) %>% pull(BinName)),
              length(common_eng %>% filter(n >= 3) %>% pull(BinName)),
              length(common_eng %>% filter(n >= 4) %>% pull(BinName)),
              length(common_eng %>% filter(n >= 5) %>% pull(BinName)),
              length(common_eng %>% filter(n >= 6) %>% pull(BinName)) ))

ggplot(MAGengraftment_curv, aes(Npts, eMAGs, group ="1", colour = "1")) +
   geom_point()+
   geom_line()+
   scale_color_manual(
   values = c("#006837"),
   labels = c("MAGs (pts=11)")) +
   #scale_y_log10() +
   scale_y_continuous(breaks = c(0, 9, 25, 50, 75, 100)) +
   #geom_hline(yintercept=15, linetype='dashed', size=0.5)+
   geom_hline(yintercept=9, linetype='dashed')+
   theme_bw()+
   theme(text = element_text(size = 10, face = "bold"),
         legend.position=c(.6,.85),
        legend.title=element_blank(),
        strip.text = element_blank()) +
   xlab("Number of patients") +
   ylab("Engrafted MAGs")
ggsave('METABAT_bin_analysis_figures/MAGengraftment_curv_fig2.png', 
  heigh = 8, width = 8, units = "cm")
 
```

```{r}
common_eng %>%
  filter(n>=3) %>% pull(BinName)-> top_MAG
  as.character(top_MAG) -> top_MAG
  
gtdb %>%
mutate(BinName= gsub("bin.", "M", user_genome)) %>%
filter(BinName %in% top_MAG) %>%
select(BinName, Order, Family, Genus, Krak_fam, Krak_gen) %>%
mutate(Fig_gen = case_when(is.na(Genus) ~ paste(Krak_gen),
                           Genus == "UMGS1375" ~ paste("Lachnospiraceae"),
                           Genus == "UBA11524" ~ paste("Faecalibacterium"),
                           Genus == "Clostridium_A" ~ paste("Clostridium"),
                           TRUE ~ as.character(Genus))) -> top_gtdb

MAG_MAB_type_FMT_viz %>%
filter(DS_model %in% c("Engraft", "StrainReplace")) %>%
filter(Timepoint == "Time6") %>%
select(Fig_lab, DS_model, Fate, BinName) %>% distinct() %>%
filter(BinName %in% top_MAG) %>%
left_join(top_gtdb, by = "BinName") %>%
ggplot(aes(BinName, Fig_lab)) +
  geom_tile(aes(fill=DS_model)) +
  facet_grid(Fate~Fig_gen, scales = "free", 
                  space = "free", switch = c("y")) +
  scale_fill_manual(values = DS_model_colors2)+
       theme_classic() +
       theme(text = element_text(size=8),
        legend.position="none",
        legend.title = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.text.y.left = element_text(angle = 0))
ggsave("METABAT_bin_analysis_figures/CommonMAGs_eng.png", 
heigh = 6, width = 12, units = "cm")


```



```{r echo=FALSE, message=FALSE, warning=FALSE}
########################################################
# creating a bin ordering file.
MAG_MAB_type_FMT %>%
  select(BinName, model) %>%
  group_by(BinName, model) %>% tally() %>% arrange(model) %>%
  spread(model, n, fill=0) %>%
  mutate(model = case_when(Engraft >= 2  ~ "Engraft",
                           Engraft == 0 & Lost >= 2  ~ "Lost",
                           TRUE ~ "Nochange")) %>%
  arrange(model) %>%
  select(BinName, model) -> MAG_MAB_items_order
write.table(MAG_MAB_items_order, 
            "MAG_MAB_Engraftment_items-order.txt", 
             sep = "\t", 
             row.names = FALSE, col.names = FALSE, quote = FALSE)

# a bin ordering file for "clean" bins, they defined here earlier as those
# bins that have proper phylum-level at least
MAG_MAB_items_order %>%
  right_join(MAG_MAB_clean_phylum_info) %>%
  select(BinName, model) %>%
  arrange(model) -> MAG_MAB_items_order_clean

MAG_MAB_items_order %>%
  right_join(MAG_MAB_clean_phylum_info) %>%
  select(BinName, model) %>%
  arrange(model) %>%
  select(BinName) -> MAG_MAB_items_order_clean_forANVIO
write.table(MAG_MAB_items_order_clean_forANVIO, 
            "MAG_MAB_Engraftment_items-order_clean.txt", 
             sep = "\t", 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
```

# Changes in bacterial Genes/functions
These are my work towards finding engrafted genes and functions. Here I've
made a one single bin for all contigs and summarized data through anvio for 
EVERYTHING profile. I have to figure out what type of information do I have for
each file. I am going to find the list of engrafted/lost genes across responders
and non-responder group and make sure to subset this by placebo samples.


```{r}
############################################################################
path7 <- paste("~/Drive2/UC_FMT/anvio_setup_DMG-CEMG-all/",
      "SAMPLES-MERGED_UCpatients_2500_newSamples_SUMMARY_DEFAULT_genecoverage/",
      "bin_by_bin/EVERYTHING/EVERYTHING-gene_detection.txt", sep = "")
gene_detect <- read.csv(path7, sep = "\t", stringsAsFactors = FALSE)
colnames(gene_detect) <- gsub("_MAP_SORTD", "", colnames(gene_detect))

#############################################################################
path8 <- paste("/media/Disk1/UC_FMT/anvio_setup_DMG-CEMG-all/Kaiju_taxonomy/",
               "gene_calls_nr.names", sep = "")
gene_taxa <- read.csv(path8, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
kaiju_cols <- c("Assignment", "gene_callers_id", "TaxId_CloseHit",
                "MatchScore", "TaxonIDS", "AccessionNumbers", "MatchingAA",
                "Lineage")
colnames(gene_taxa) <- kaiju_cols
gene_taxa %>%
  select(gene_callers_id, Lineage) %>%
  separate(Lineage, c("superkingdom","phylum","order","class",
                      "family","genus","species")) -> gene_taxa_clean

###############################################################################
path9 <- paste("/media/Disk1/UC_FMT/anvio_setup_DMG-CEMG-all/",
"SAMPLES-MERGED_UCpatients_2500_newSamples_SUMMARY_DEFAULT_genecoverage/",
"bin_by_bin/EVERYTHING/EVERYTHING-gene_calls.txt", sep = "")
gene_func <- read.csv(path9, sep = "\t", stringsAsFactors = FALSE)
gene_func %>% select(gene_callers_id, dna_sequence) -> gene_dna_seq
gene_func$dna_sequence <- NULL
colnames(gene_func)[2] <- "Contig"

gene_func %>%
  mutate(length= stop - start) -> gene_func_length

################################################################################
path10 <- paste("/media/Disk1/UC_FMT/anvio_setup_DMG-CEMG-all/",
"SAMPLES-MERGED_UCpatients_2500_newSamples_SUMMARY_DEFAULT_genecoverage/",
"bin_by_bin/EVERYTHING/EVERYTHING-gene_non_outlier_coverages.txt", sep = "")
gene_cove <- read.csv(path10, sep = "\t", stringsAsFactors = FALSE)
colnames(gene_cove) <- gsub("_MAP_SORTD", "", colnames(gene_cove))

gene_cove %>%
  gather(SampleName, coverage_value, -gene_callers_id) %>% 
  left_join(mapfile, by = "SampleName") %>%
  #removing duplicate samples (those that their depth increased)
  filter(!SampleName %in% c("s721_CCGTCC", "s797_GTGGCC")) %>%
  select(c("gene_callers_id", "coverage_value", "Timepoint", "Initial")) %>% 
  #adding "Time" term so it wont become inteager
  mutate(Timepoint = paste("Time", Timepoint, sep = "")) -> gene_cove_long

################################################################################

Patient_fate <- data.frame(Initial = c("PM","AF","DAB","KG","WAF","BN",
                                       "TJD","PT","MP","EB","ES"),
                           Fate = c("R", "R", "P", "R", "N", "R", "N",
                                    "R", "N", "N", "R"))

# converting to long format so that can be merged with mapfile
# clean up and convert back to wide for engraftment detection
 gene_detect %>% 
  gather(SampleName, detect_value, -gene_callers_id) %>%
  left_join(mapfile, by = "SampleName") %>%
  #removing duplicate samples (those that their depth increased)
  filter(!SampleName %in% c("s721_CCGTCC", "s797_GTGGCC")) %>%
  select(c("gene_callers_id","detect_value", "Timepoint", "Initial")) %>%
  #adding "Time" term so it wont become inteager
  mutate(Timepoint = paste("Time", Timepoint, sep = "")) %>%
  # make it wide to easily select egrafted bins
  spread(key=Timepoint, value=detect_value) %>%
  group_by(Initial) %>%
  # define the type of bins in terms of their change over FMT
  mutate(model= case_when(Time0 <= 0 & Time6 >= 0.6 ~ "Engraft",
                            Time0 >= 0.6 & Time6 <= 0 ~ "Lost",
                            TRUE ~ "Nochange")) %>%
  # make it long again
  gather(Timepoint, detect_value, -gene_callers_id, -Initial, -model) %>%
  # adding patient fate  
  left_join(y=Patient_fate, by = "Initial") %>%
  # adding the non-outlier coverage information
  left_join(gene_cove_long, by = c("gene_callers_id", 
                                   "Initial", "Timepoint")) -> gene_type_FMT
 

```



```{r}

  ########################## 
  # FOR RESPONDERS
  ##########################
  # Subseting to engrafted AND
  # only to responders patient and placebo sample
 gene_type_FMT %>% filter(model == "Engraft" & !Fate == "N") %>%
   filter(Timepoint == "Time6") %>%
   # filter only those genes that engrafted with minimum 5X coverage
   filter(coverage_value >= 5) %>%
   group_by(gene_callers_id) %>%
   select(-coverage_value) %>%
   spread(Initial, detect_value) %>%
   #convert NA To 0 for all
   #filter(DAB >= 0.6) -> test 
   replace(., is.na(.), 0) %>%
   select(-Fate) -> gene_engraft_res 


  # the engraftment gene must be common across at least 3 Responders
  # with no overlap with placebo
  # which means half of responder patients
 gene_engraft_res[rowSums(gene_engraft_res[,4:10] > 0.6) >= 3, ] -> gene_engraft_res_com
 gene_engraft_res_com %>%
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id" ) %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") -> gene_engraft_res_com_complete
 
################################################################################ 
  ######################### 
  # FOR NON-RESPONDERS ####
  #########################
  # Subseting to engrafted AND
  # only to non-responders patient and placebo sample
 
 gene_type_FMT %>% filter(model == "Engraft" & !Fate == "R") %>%
   filter(Timepoint == "Time6") %>%
   filter(coverage_value >= 5) %>%
   group_by(gene_callers_id) %>%
   select(-coverage_value) %>%
   spread(Initial, detect_value) %>%
   #filter(DAB >= 0.6) -> test 
   replace(., is.na(.), 0) %>%
   select(-Fate) -> gene_engraft_Nores
 
 # the engraftment gene must be common across at least 3 Non-Responders
  # with no overlap with placebo
 gene_engraft_Nores[rowSums(gene_engraft_Nores[,4:8] > 0.6) >= 3, ] -> gene_engraft_Nores_com
 
 gene_engraft_Nores_com %>%
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id" ) %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") -> gene_engraft_Nores_com_complete

  #############################################################################
  ####################################
  # FOR ALL Patients ENGRAFTED #######
  ####################################
  # Subseting to engrafted for all
 gene_type_FMT %>% filter(model == "Engraft") %>%
   filter(Timepoint == "Time6") %>%
   filter(coverage_value >= 5) %>%
   group_by(gene_callers_id) %>%
   select(-coverage_value) %>%
   spread(Initial, detect_value) %>%
   #convert NA To 0 for all
   replace(., is.na(.), 0) %>%
   select(-Fate) -> gene_engraft_all

 ##############################################################################
 # TEST DATASET, TO Make engraftment curve**
 # 2patients:
 gene_engraft_all[rowSums(gene_engraft_all[,4:14] > 0.6) >= 2, ] -> gene_engraft_all_com_2
 # 3patients:
 gene_engraft_all[rowSums(gene_engraft_all[,4:14] > 0.6) >= 3, ] -> gene_engraft_all_com_3
 # 4patients:
 gene_engraft_all[rowSums(gene_engraft_all[,4:14] > 0.6) >= 4, ] -> gene_engraft_all_com_4
 # 5patients:
 gene_engraft_all[rowSums(gene_engraft_all[,4:14] > 0.6) >= 5, ] -> gene_engraft_all_com_5
 # 6patients:
 gene_engraft_all[rowSums(gene_engraft_all[,4:14] > 0.6) >= 6, ] -> gene_engraft_all_com_6
 
 engraftment_curv <- data.frame(Npts = c("1", "2", "3", "4", "5", "6"),
    Egenes = c(length(gene_engraft_all$gene_callers_id), 
               length(gene_engraft_all_com_2$gene_callers_id), 
               length(gene_engraft_all_com_3$gene_callers_id), 
               length(gene_engraft_all_com_4$gene_callers_id), 
               length(gene_engraft_all_com_5$gene_callers_id), 
               length(gene_engraft_all_com_6$gene_callers_id)))
 
 
 #engraftment_curv_fig <- 
ggplot(engraftment_curv, aes(Npts, Egenes, group ="1", colour = "1")) +
   geom_point()+
   geom_line()+
   scale_color_manual(
   values = c("#006837"),
   labels = c("Genes (pts=11)"))+
   #scale_y_sqrt() +
   scale_y_log10() +
   geom_hline(yintercept=267, linetype='dashed')+
   #geom_hline(yintercept=log(7), linetype='dashed', size=1)+
   #geom_hline(yintercept=log(267), linetype='dashed', size=0.5)+
   #geom_hline(yintercept=log(13092), linetype='dashed', size=0.5)+
   #scale_y_continuous(breaks = c(0, 7, 267, 13092, 139535)) +
   theme_bw()+
   theme(text = element_text(size = 10, face = "bold"),
         legend.position=c(.6,.85),
        legend.title=element_blank(),
        strip.text = element_blank()) +
   xlab("Number of patients") +
   ylab("Engrafted genes")
 ggsave('METABAT_bin_analysis_figures/Common_engraftmentCurve.png', 
  heigh = 8, width = 8, units= "cm")



```



```{r}

###############################################################################
  # this is a big table of all data together:
  #It's not related to this analysis but I wanted a table of gene, contig, and
  # bin information:
  gene_type_FMT %>%
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id" ) %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") %>%
   mutate(length = as.numeric(stop) - as.numeric(start))
  

  
 # the engraftment gene must be common across at least 3 patient
 # with no overlap with placebo
 gene_engraft_all[rowSums(gene_engraft_all[,4:14] > 0.6) >= 3, ] -> gene_engraft_all_com
 
 gene_engraft_all_com %>%
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id" ) %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") %>%
   mutate(length = as.numeric(stop) - as.numeric(start)) -> gene_engraft_all_com_complete
 
 # summary of Kaiju assignment on the genes
 gene_engraft_all_com_complete %>%
   group_by(genus.x) %>% tally() %>% top_n(10, n)
 
 # summary of contigs for the engrafted genes:
 gene_engraft_all_com_complete %>%
   group_by(Contig) %>% tally() %>% top_n(10, n)
 
 #Summary of MAG/MAB for these genes:
 gene_engraft_all_com_complete %>%
   select(length)
   

```


```{r}

 gene_engraft_all %>%
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id" ) %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") %>%
   mutate(length = as.numeric(stop) - as.numeric(start)) -> gene_engraft_all_NonCom_complete
   
  # summary of Kaiju assignment on the genes
gene_engraft_all_NonCom_complete %>%
   group_by(genus.x) %>% tally() %>% top_n(10, n)
 
 
 ###############################################################################
# write the 267 genes ids into a file to get aa sequences from DonB database:
 gene_engraft_all_com_complete %>%
   select(gene_callers_id) -> gene_callers_ids 
 write.csv(gene_callers_ids, "gene_engraft_all_com_only-ids.txt",
           sep = ",", col.names = FALSE, row.names = FALSE,
           quote = FALSE)
 ###############################################################################
# write all the 139,535 genes ids into a file to get aa sequences from DonB database:
 gene_engraft_all_NonCom_complete %>%
   select(gene_callers_id) -> gene_callers_ids_NonCom
 write.csv(gene_callers_ids_NonCom, "gene_engraft_all_NonCom-ids.txt",
           sep = ",", col.names = FALSE, row.names = FALSE,
           quote = FALSE)
 
 #############################################################################
 ######## WRITING THE COMMON ENGRAFTED GENES INFO TABLE AND SEQUENCE FILE
 write.table(gene_engraft_all_com_complete,
           "gene_engraft_all_com_complete.txt", row.names = FALSE, 
           quote = FALSE, sep = "\t")
 gene_engraft_all_com_complete %>%
   select(gene_callers_id) %>%
   left_join(gene_dna_seq, by = "gene_callers_id") -> gene_engraft_dna_seq 
 
 #############################################################################
 ####### Writing a clean table of COMMON ENGRAFTED GENES for manual curation
 
 coloumns_of_int <- c("gene_callers_id", "Contig", "start","stop","direction",
                      "COG_FUNCTION","COG_FUNCTION..ACCESSION.","KeggGhostKoala",
                      "KeggGhostKoala..ACCESSION.","COG_CATEGORY",
                      "COG_CATEGORY..ACCESSION.",
                      "phylum.x","order.x","family.x","genus.x","species.x",
                      "BinName","Bin.Id","Quality","Rank", "TaxName",
                      "phylum.y", "order.y", "family.y","genus.y","species.y",
                      "length")
 
 gene_engraft_all_com_complete %>%
   select(coloumns_of_int) %>%
   mutate(Header= paste(genus.x, gene_callers_id, sep = "_")) %>%
   mutate(Header= paste(Header, COG_FUNCTION..ACCESSION., 
                        sep = "|cog:")) %>%
   mutate(Header= paste(Header, KeggGhostKoala..ACCESSION., 
                        sep = "|kegg:")) %>%
   mutate(gene_callers = paste(genus.x, gene_callers_id, sep = "_")) %>%
   arrange(gene_callers) -> gene_engraft_all_com_complete_edited
 
 write.table(gene_engraft_all_com_complete_edited,
           "gene_engraft_all_com_complete_edited.txt", row.names = FALSE, 
           quote = FALSE, sep = "\t")
 


 #############################################################################
 ######## WRITING ALL THE ENGRAFTED GENES INFO TABLE AND SEQUENCE FILE
 write.table(gene_engraft_all_NonCom_complete,
           "gene_engraft_all_NonCom_complete.csv", row.names = FALSE, 
            quote = FALSE, sep = "\t")
gene_engraft_all_NonCom_complete %>%
   select(gene_callers_id) %>%
   left_join(gene_dna_seq, by = "gene_callers_id") -> gene_engraft_dna_seq_NonCom 


 # A FUNCTION TO WRITE FASTA FILE FOR GENE
 writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", 
                                                  data[rowNum,"gene_callers_id"], 
                                                  sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"dna_sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
 }
 
 #writeFasta(gene_engraft_dna_seq_NonCom, "gene_engraft_all_NonCom_complete.fasta")
 
 ###############################################################################
 # READING AA seqs of the commonly engrafted genes to make "faa"
 ###############################################################################
 gene_engraft_aa_seq <- ampir::read_faa("gene_engraft_all_com_aa.fasta")
 gene_engraft_aa_seq %>%
   dplyr::rename(gene_callers_id=seq_name,
          dna_sequence=seq_aa) -> gene_engraft_aa_seq
gene_engraft_aa_seq$gene_callers_id <- as.integer(gene_engraft_aa_seq$gene_callers_id)
####################################
 writeFasta(
 gene_engraft_all_com_complete %>%
   mutate(Header= paste(genus.x, gene_callers_id, sep = "_")) %>%
   mutate(Header= paste(Header, COG_FUNCTION..ACCESSION., 
                        sep = "|cog:")) %>%
   mutate(Header= paste(Header, KeggGhostKoala..ACCESSION., 
                        sep = "|kegg:")) %>%
   select(gene_callers_id, Header) %>%
   left_join(gene_dna_seq, by = "gene_callers_id") %>%
   ungroup(gene_callers_id) %>%
   mutate(gene_callers_id = Header),
   "gene_engraft_all_com_complete.fasta")
#################################### 
writeFasta(
 gene_engraft_all_com_complete %>%
   mutate(Header= paste(genus.x, gene_callers_id, sep = "_")) %>%
   mutate(Header= paste(Header, COG_FUNCTION..ACCESSION., 
                        sep = "|cog:")) %>%
   mutate(Header= paste(Header, KeggGhostKoala..ACCESSION., 
                        sep = "|kegg:")) %>%
   select(gene_callers_id, Header) %>%
   left_join(gene_engraft_aa_seq, by = "gene_callers_id") %>%
   ungroup(gene_callers_id) %>%
   mutate(gene_callers_id = Header),
  "gene_engraft_all_com_complete.faa")

```



```{r}

##############################################################################
 mapfile %>%
   select(Initial, Fig_lab) -> mapfile_B
 ##############################################################################
 ### RESULTS: 7 COMMON GENES (AT LEAST IN 4 PATIENTS) ARE ENGRAFTED IN TOTAL
 
 gene_engraft_all_com_4 %>%
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id" ) %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") %>%
   mutate(length = as.numeric(stop) - as.numeric(start)) -> gene_engraft_all_com_4_complete
 

 gene_engraft_all_com_4 %>%
   gather(Initial, detect_value, -gene_callers_id, -model, -Timepoint) %>% 
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id" ) %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") %>%
   ungroup() %>%
   mutate(gene_callers_id = paste("g_", gene_callers_id, "")) %>%
   left_join(Patient_fate, by = "Initial") %>%
   left_join(mapfile_B, by = "Initial") -> gene_engraft_all_4pt_visual


ggplot(gene_engraft_all_4pt_visual, 
       aes(gene_callers_id, Initial)) +
   geom_tile(aes(fill = detect_value)) +
   facet_grid(Fate+Fig_lab~BinName+Quality+genus.x,  space = "free", scales = "free") +
   scale_fill_gradient(low = "white", high = "#00468b") +
   theme_bw() +
   theme(#axis.title.x = element_blank(),
        text = element_text(size=10),
        #legend.position="bottom",
        legend.position="none",
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        strip.text.y = element_text(size = 10, angle = 0)) +
  ggtitle("Common engrafted genes (minimum 4 individuals) across all PT")+
  xlab("7 genes")

```




```{r}

###############################################################################
 ###RESULT: 267 COMMON GENES (AT LEAST IN 3 PATIENTS) ARE ENGRAFTED IN TOTAL
 ########## 265 OF THOSE ARE PRESENT IN RESPONDERS PATIENTS AND ONLY 2 IN NON-
 ########## RESPONDER GROUP. ALSO, NO COMMON ENGRAFTED GENE FOUND FOR PLACEBO
 ########## SAMPLE.

 #Visualizing all the engrafted genes (those common at least across 3 PT)
 #for all patients at timepoint6
 gene_engraft_all_com %>%
   gather(Initial, detect_value, -gene_callers_id, -model, -Timepoint) %>%
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id" ) %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") %>%
   ungroup() %>%
   mutate(gene_callers_id = paste("g_", gene_callers_id, "")) %>%
   left_join(Patient_fate, by = "Initial") %>%
   left_join(mapfile_B, by = "Initial") %>%
   mutate(Detection= case_when(detect_value == 0 ~ paste("0.0"),
                              between(detect_value, 0.2, 0.4) ~ paste("0.2 - 0.4"),
                              between(detect_value, 0.4, 0.6) ~ paste("0.4 - 0.6"),
                              between(detect_value, 0.6, 0.8) ~ paste("0.6 - 0.8"),
                              between(detect_value, 0.8, 1) ~ paste("0.8 - 1"),
                              TRUE ~ as.character(detect_value))) -> gene_engraft_all_visual

 
 
 ##Creating a dataset contain both Time0 and Time6 for engrafted genes
 gene_engraft_all_com %>% select(gene_callers_id) %>%
   ungroup() -> gene_of_interest
 dplyr::pull(gene_of_interest, gene_callers_id) -> gene_of_interest
 
 
 ##############################################################################
 # FIGURE 8 - engrafted genes - only Timeppoint 6 for each pt
 ##############################################################################
 
Detection_cols3 <-  c(#"0.0" = "#ffffcc",
                     "0.2 - 0.4" = "#a1dab4",
                     "0.4 - 0.6" = "#41b6c4",
                     "0.6 - 0.8" = "#2c7fb8", 
                     "0.8 - 1" = "#253494")

ggplot(gene_engraft_all_visual %>%
          mutate(family.x = case_when(is.na(family.x) ~ "NA",
                                      TRUE ~ as.character(family.x))) %>%
          mutate(family.x = case_when(family.x == "NA" ~ paste("Unknown_family"),
                                      TRUE ~ as.character(family.x))) %>%
          mutate(gene_callers_id = paste(genus.x, gene_callers_id, sep = "_")) %>%
          mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)), 
       aes(gene_callers_id, Fig_lab)) +
  geom_tile(aes(fill= Detection), alpha = 8/10) +
     scale_fill_manual(values= Detection_cols3) +
   facet_grid(~Fate+Fig_lab~family.x, scale = "free", space = "free", switch = "y") +
    theme(#axis.title.x = element_blank(),
        text = element_text(size=15),
        legend.position="bottom",
        #legend.position="none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90),
        strip.text.x = element_text(size = 8, angle = 90),
        strip.text.y = element_text(size = 15, angle = 0),
         strip.placement = "none")
 ggtitle("Common engrafted genes (minimum 3 individuals) across all PT")

  


```



```{r}

  
#############################################################################
# engrafted genes - only Timeppoint 6 for each pt- ** ALSO THOSE THE ALREADY
# PRESENTS
#############################################################################
  
  
  gene_engraft_all_visual %>%
    select(-detect_value) %>%
    #select(gene_callers_id, Initial, Fig_lab, Fate, Detection) %>%
    distinct() %>%
    mutate(Detection = paste("Engraft", Detection, sep = "_")) %>%
    mutate(gene_callers_id = gsub(" ", "", gene_callers_id)) %>%
    # Whcih of these genes were already present?**********
    left_join(gene_type_FMT %>%
                filter(Timepoint == "Time6") %>%
                filter(model == "Nochange") %>%
                filter(detect_value >= 1) %>%
                filter(coverage_value >= 10) %>%
                mutate(gene_callers_id = paste("g_", gene_callers_id, sep = "")) %>%
                select(gene_callers_id, Initial, detect_value)) %>%
    mutate(detect_value= case_when(detect_value == 0 ~ paste("0.0"),
                              between(detect_value, 0.2, 0.4) ~ paste("0.2 - 0.4"),
                              between(detect_value, 0.4, 0.6) ~ paste("0.4 - 0.6"),
                              between(detect_value, 0.6, 0.8) ~ paste("0.6 - 0.8"),
                              between(detect_value, 0.8, 1) ~ paste("0.8 - 1"),
                              TRUE ~ as.character(detect_value))) %>%
    mutate(detect_value = paste("Present", detect_value, sep = "_")) %>%
    mutate(Value = case_when(Detection == "Engraft_0.0" & 
                               detect_value != "Present_NA" ~ paste(detect_value),
                             TRUE ~ as.character(Detection))) %>%
    # Highlighting those genes that engrafted in 4 patients#*********
    left_join(
    gene_engraft_all_4pt_visual %>%
    select(gene_callers_id, Initial, Fig_lab) %>%
    distinct() %>%
    mutate(gene_callers_id = gsub(" ", "", gene_callers_id)) %>%
    mutate(GeneType= "SuperGene" )) %>%
    mutate(Value = case_when(GeneType == "SuperGene" &
                               Value == "Engraft_0.8 - 1" ~ paste("*", Value, sep = ""),
                             TRUE ~ as.character(Value))) ->  gene_engraft_present_visual
  
                
  
  Value_cols3 <-  c(#"Engraft_0.0" = "#ffffcc",
                     "Engraft_0.2 - 0.4" = "#a1dab4",
                     "Engraft_0.4 - 0.6" = "#41b6c4",
                     "Engraft_0.6 - 0.8" = "#2c7fb8", 
                     "Engraft_0.8 - 1" = "#253494",
                     "*Engraft_0.8 - 1" = "red",
                     "Present_0.8 - 1" = "#252525")

ggplot(gene_engraft_present_visual %>%
          mutate(family.x = case_when(is.na(family.x) ~ "NA",
                                      TRUE ~ as.character(family.x))) %>%
          mutate(family.x = case_when(family.x == "NA" ~ paste("Unknown_family"),
                                      TRUE ~ as.character(family.x))) %>%
          mutate(gene_callers_id = paste(genus.x, gene_callers_id, sep = "_")) %>%
          mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)), 
       aes(gene_callers_id, Fig_lab)) +
  geom_tile(aes(fill= Value), alpha = 8/10) +
     scale_fill_manual(values= Value_cols3) +
   facet_grid(~Fate+Fig_lab~family.x, scale = "free", space = "free", switch = "y") +
    theme(#axis.title.x = element_blank(),
        text = element_text(size=15),
        legend.position="bottom",
        #legend.position="none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90),
        strip.text.x = element_text(size = 8, angle = 90),
        strip.text.y = element_text(size = 15, angle = 0),
         strip.placement = "none")
 ggtitle("Common engrafted genes (minimum 3 individuals) across all PT")

  
```





```{r}
ggplot(gene_engraft_all_visual, 
       aes(Fig_lab, gene_callers_id)) +
   geom_tile(aes(fill = detect_value)) +
   scale_fill_gradient(low = "white", high = "#00468b") +
   #facet_nested(Quality+BinName+TaxName~Fate, scales = "fixed", space = "fixed") +
   facet_grid(Quality+BinName+TaxName~Fate, scales = "free", space = "free") +
   theme_bw() +
   #coord_flip() +
   theme(#axis.title.x = element_blank(),
        text = element_text(size=15),
        legend.position="right",
        #legend.position="none",
        #axis.text.y = element_blank(), 
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5),
        #axis.text.x = element_text(size = 10, angle = 35),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8, angle = 0)) +
  ggtitle("Common engrafted genes (minimum 3 individuals) across all PT by Bin taxa")



```


```{r}


gene_engraft_all_visual %>%
  mutate(FUNCTION = case_when(KeggGhostKoala== "" &
                                    COG_FUNCTION=="" ~ paste("Unknown KEGG or COG"),
                                    !KeggGhostKoala=="" &
                                    !COG_FUNCTION=="" ~ paste("KEGG + COG"),
                                    !KeggGhostKoala=="" ~ paste("KEGG"),
                                    !COG_FUNCTION=="" ~ paste("COG"),
            TRUE ~ as.character(KeggGhostKoala))) -> gene_engraft_all_visual2
  

 ggplot(gene_engraft_all_visual2, 
       aes(Initial, gene_callers_id)) +
   geom_tile(aes(fill = detect_value)) +
   scale_fill_gradient(low = "white", high = "#00468b") +
   facet_grid(~FUNCTION~Fate, scale = "free", space = "free") +
   theme_bw()+
   theme(text = element_text(size=15),
        legend.position="right",
        #legend.position="none",
        #axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5),
         strip.text.x = element_text(angle=0),
         strip.text.y = element_text(size = 8, angle = 0))+
  ggtitle("Common engrafted genes (minimum 3 individuals) across all PT by KEGG")

  

```


```{r}
 
 #############################################################################
  ####################################
  # FOR ALL Patients LOST #######
  ####################################
  # Subseting to engrafted for all
 gene_type_FMT %>% filter(model == "Lost") %>%
   filter(Timepoint == "Time0") %>%
   filter(coverage_value >= 5) %>%
   group_by(gene_callers_id) %>%
   select(-coverage_value) %>%
   spread(Initial, detect_value) %>%
   #convert NA To 0 for all
   replace(., is.na(.), 0) %>%
   select(-Fate) -> gene_lost_all

 
 # the engraftment gene must be common across at least 3 patient
 # with no overlap with placebo
 gene_lost_all[rowSums(gene_lost_all[,4:14] > 0.6) >= 3, ] -> gene_lost_all_com
 
 gene_lost_all_com %>%
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id") %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") %>%
   mutate(length = as.numeric(stop) - as.numeric(start)) -> gene_lost_all_com_complete
 
 summary(gene_lost_all_com_complete$length)
 
 gene_lost_all_com %>%
   gather(Initial, detect_value, -gene_callers_id, -model, -Timepoint) %>%
   left_join(gene_func, by = "gene_callers_id") %>%
   left_join(gene_taxa_clean, by = "gene_callers_id" ) %>%
   left_join(ContigInBins, by = "Contig") %>%
   left_join(BinsFine_complete, by = "BinName") %>%
   ungroup() %>%
   mutate(gene_callers_id = paste("g_", gene_callers_id, "")) %>%
   left_join(Patient_fate, by = "Initial") %>%
   mutate(FUNCTION = case_when(KeggGhostKoala== "" &
                                    COG_FUNCTION=="" ~ paste("Unknown KEGG or COG"),
                                    !KeggGhostKoala=="" &
                                    !COG_FUNCTION=="" ~ paste("KEGG + COG"),
                                    !KeggGhostKoala=="" ~ paste("KEGG"),
                                    !COG_FUNCTION=="" ~ paste("COG"),
            TRUE ~ as.character(KeggGhostKoala))) -> gene_lost_all_visual
 
 
ggplot(gene_lost_all_visual, 
       aes(Initial, gene_callers_id)) +
   geom_tile(aes(fill = detect_value)) +
   scale_fill_gradient(low = "white", high = "darkred") +
   facet_grid(Quality+BinName+TaxName~Fate, scale = "free", space = "free") +
   theme_bw()+
    theme(text = element_text(size=15),
        legend.position="right",
        #legend.position="none",
        #axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
         axis.text.y = element_text(size = 5),
         strip.text.x = element_text(angle=0),
         strip.text.y = element_text(size = 10, angle = 0))+
  ggtitle("Common lost genes (minimum 3 individuals) across all PT by bin")



```

```{r}

ggplot(gene_lost_all_visual, 
       aes(Initial, gene_callers_id)) +
   geom_tile(aes(fill = detect_value)) +
   scale_fill_gradient(low = "white", high = "darkred") +
   facet_grid(~FUNCTION~Fate, scale = "free", space = "free") +
   theme_bw()+
    theme(text = element_text(size=15),
        legend.position="right",
        #legend.position="none",
        #axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5),
         strip.text.x = element_text(angle=0),
         strip.text.y = element_text(size = 10, angle = 0))+
  ggtitle("Common lost genes (minimum 3 individuals) across all PT by Function")


```



```{r echo=FALSE, message=FALSE, warning=FALSE}


###############################################
# Extracting 5000 Up/down stream genes across *commonly engrafted gene*:

  gene_engraft_all_com_complete %>%
    select(gene_callers_id, Contig, start, stop, length) %>%
    mutate(ContigInfo=Contig) %>%
    separate(ContigInfo, c("S1", "S2", "S3", "S4", "S5", "S6"), "_") %>%
    dplyr::rename(ContigLength=S4) %>%
    select(gene_callers_id, Contig, ContigLength, start, stop, length) %>%
    # Calculate Up/Down stream positions
    mutate(DownStart=start - 5000) %>%
    mutate(Upstop=stop + 5000) %>%
    # For those contig that their length is less than 10,000 > grab the whole
    mutate(DownStart = case_when(as.numeric(ContigLength) <= 10000 ~ paste(1),
                                 TRUE ~ as.character(DownStart))) %>%
    mutate(Upstop = case_when(as.numeric(ContigLength) <= 10000 ~ 
                                   paste(as.numeric(ContigLength)),
                                 TRUE ~ as.character(Upstop))) %>%
    # For those contigs which have less than 5000 downstream genes  
    mutate(DownStart = case_when(as.numeric(DownStart) < 0 ~ paste(1),
                                 TRUE ~ as.character(DownStart))) %>%
    # For those contigs which have less than 5000 upstream genes
    mutate(Upstop= case_when(as.numeric(Upstop) > as.numeric(ContigLength) ~ 
                               paste(as.numeric(ContigLength)),
                             TRUE  ~ as.character(Upstop))) %>%
    mutate(OperonLength = as.numeric(Upstop) - as.numeric(DownStart)) -> UpDownGenes_complete

# Write the new start and stop position to a new file so that can be extracted  
ungroup(UpDownGenes_complete) %>% 
  mutate(Orient= paste(Contig, DownStart, sep = ":")) %>%
  mutate(Orient= paste(Orient, Upstop, sep = "-")) %>% 
  select(Orient) -> UpDownGene_f 
  
write.table(UpDownGene_f, "gene_engraft_all_com_UpDownStream.txt", sep = "\t",
          row.names = FALSE, col.names = FALSE, quote = FALSE)
  
# then run this below code to extract the oprean
# samtools faidx contigs-fixed.fa -r UpDownGene-f -o UpDown.fasta
  
  
```

# Functiona/PATHWAY
Here I retrieved the amino acid sequeces for the 267 genes that engrafted
among three patients and ran them against uniref90 via diamond blastp.
Then the diamond output imported here, low evalue hits removed and the uniref90
ids written to file to be uploaded in the uniprot website. The Functions and 
pathway as well as other inforamtion retrieved and saved into a file. The final
result updated adain and a result table generated. 

by default evalue cut-off for Diamond is 0.001 the reported result are already
more stringent than blast. So I don't trim of any result and upload all these
Unire90 ids to the Uniprot website.


```{r}

dmnd_gene_engraft_all_com <- read.csv("gene_engraft_all_com_complete_faa_maxtarget1.dimndoute",
                                      sep = "\t", header = FALSE, stringsAsFactors = F)
dmnd_cols <- c("qseqid", "sseqid", "pident", "length", 
"mismatch", "gapopen", "qstart", "qend", "sstart", "send",
"evalue", "bitscore")
colnames(dmnd_gene_engraft_all_com) <- dmnd_cols


dmnd_gene_engraft_all_com %>%
  filter(pident >= 40) %>% #set this min to remove 25 duplication of a single gene
  select(sseqid) -> gene_engraft_all_com_uniref90
  #mutate(sseqid = gsub("UniRef90_", "", sseqid)) -> gene_engraft_all_com_uniref90
write.table(gene_engraft_all_com_uniref90,
          "gene_engraft_all_com_aa_uniref90-ids.txt", sep = ",",
          row.names = FALSE, col.names = FALSE, quote = FALSE)

path11 <- "Uniprot_outputs/gene_engraft_all_com_aa_uniref90_uniprot.txt"
Uniprot_gene_engraft_all <- read.csv(path11, sep = "\t", 
                                     stringsAsFactors = FALSE)
Uniprot_gene_engraft_all %>%
  mutate(sseqid= paste("UniRef90", Entry, sep = "_")) -> Uniprot_gene_engraft_all

gene_engraft_all_com_complete_edited %>% #the complete common engrafted list
left_join( 
  dmnd_gene_engraft_all_com %>%
  filter(pident >= 40) %>%#set this min to remove 25 duplication of a single gene
  select(qseqid, sseqid, pident) %>%
  left_join(Uniprot_gene_engraft_all, by = "sseqid") %>%
  select(-Entry, -Entry.name, -Status) %>%
  dplyr::rename(Header=qseqid),  by = "Header") -> gene_engraft_all_com_complete_edited_uniprot

write.table(gene_engraft_all_com_complete_edited_uniprot,
          "gene_engraft_all_com_complete_edited_uniprot.txt", sep = "\t",
          row.names = FALSE, quote = FALSE)


#Manually editi the  taxonomy using all information  and  read the file
edited_uniprot2 <- read.csv("gene_engraft_all_com_complete_edited_uniprot_Mtaxa.txt",
                            sep = "\t",  stringsAsFactors = F)

edited_uniprot3 <- read.csv("gene_engraft_all_com_complete_edited_uniprot_Mtaxa_pfam.txt",
                            sep = "\t",  stringsAsFactors = F)



```



```{r}

gene_engraft_all_com_complete_edited_uniprot %>%
  select(gene_callers_id:direction, COG_FUNCTION, COG_CATEGORY) %>%
  left_join(edited_uniprot3 %>% select(gene_callers_id, Protein.names)) %>%
  mutate(protein.type= case_when(is.na(Protein.names) ~ paste("Hypothetical"),
                         Protein.names== "Uncharacterized protein" ~  paste("Hypothetical"),
                         TRUE ~ paste("Known protein"))) %>%
  mutate(temp= case_when(COG_CATEGORY == "" ~ paste(protein.type),
                         TRUE ~ paste(COG_CATEGORY))) %>%
  mutate(temp= gsub("R!!!", "", temp)) %>%
  mutate(temp= gsub("!!!R", "", temp)) %>%
  mutate(temp= gsub("S!!!", "", temp)) %>%
  mutate(temp= gsub("!!!S", "", temp)) %>%
  mutate(COG_code= gsub("!!!.*", "", temp)) %>%
  group_by(COG_code) %>% tally() %>%
  # case R or S from COG are general function or unknown function, similar to
  # known protein category. we know their protein name but not funciton;
  mutate(COG_code= case_when(COG_code %in% c("R", "S") ~ paste("Known protein"),
                             TRUE ~ paste(COG_code))) %>%
  mutate(COG_edit= case_when(n < 5 ~ paste("Other"),
                             TRUE ~ paste(COG_code))) %>%
  group_by(COG_edit) %>% summarise(n= sum(n)) -> COG_pie

COG_pie$COG_edit <- factor(COG_pie$COG_edit, 
                           levels = c("Hypothetical", "Known protein","K", "J",
                                      "E", "L", "G", "H","Other"))

COG_COL <- c("Hypothetical" = "#bdbdbd", 
             "Known protein"= "#80b1d3",
             "K" = "#fdb462", 
             "J" = "#8dd3c7",
             "E" = "#ffffb3", 
             "L" = "#bebada", 
             "G" = "#fb8072", 
             "H" = "#fccde5",
             "Other" = "#b3de69")


COG_pie %>%
  mutate(Type= paste("CEGs")) %>%
  ggplot(aes(Type,n, fill=COG_edit)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = COG_COL) +
  #scale_y_continuous(breaks  = c(0, 100, 200, 255)) +
  theme_classic() +
  geom_text(aes(label=n), position = position_stack(vjust = 0.5),
            size = 3) +
  theme(text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())
ggsave("METABAT_bin_analysis_figures/COG_pie.png", width = 9, height = 8,
       units = "cm")



```




```{r}

##############################################################################
 # FIGURE 8 - engrafted genes - only Timeppoint 6 for each pt
 ##############################################################################
 
Detection_cols3 <-  c(#"0.0" = "#ffffcc",
                     "0.2 - 0.4" = "#c2e699",
                     "0.4 - 0.6" = "#78c679",
                     "0.6 - 0.8" = "#31a354", 
                     "0.8 - 1" = "#006837")


 #gene_engraft_Figure72_families <- 
   ggplot(gene_engraft_all_visual %>%
         #remove whitespace
         mutate(gene_callers_id= gsub(" ", "", gene_callers_id)) %>%
         #join uniprot info
          left_join(
          edited_uniprot2 %>%
          select(gene_callers_id, Protein.names, Family, Family_label, Genus) %>%
          mutate_if(is.character, trimws) %>%
          mutate(Family_label = factor(Family_label, levels = c("Lachnospiraceae",
                                                                "Ruminococcaceae",
                                                                "Other"))) %>%
          mutate(gene_callers_id = paste("g", gene_callers_id, sep = "_")),
          by = "gene_callers_id") %>%
          #mutate(Family = case_when(family.x == "NA" ~ paste("Unknown_family"),
          #                            TRUE ~ as.character(family.x))) %>%
          mutate(gene_callers_id = paste(Genus, gene_callers_id, sep = "_")) %>%
          mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)), 
          #mutate(gene_callers_id = paste(Protein.names, gene_callers_id, sep = "_")) %>%
          #mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)), 
        aes(gene_callers_id, Fig_lab)) +
        geom_tile(aes(fill= Detection), alpha = 8/10) +
        scale_fill_manual(values= Detection_cols3) +
   facet_grid(~Fate+Fig_lab~Family_label, scale = "free", space = "free", switch = "y") +
   theme_classic()  +
    theme(#axis.title.x = element_blank(),
        text = element_text(size=10),
        legend.position="bottom",
        #legend.position="none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text= element_text(face ="bold"),
         strip.placement = "none")
  ggsave('METABAT_bin_analysis_figures/gene_engraft_Figure72_families.png', 
  heigh = 12, width = 18, units = "cm")

  

```


```{r}

#############################################################################
# engrafted genes - only Timeppoint 6 for each pt- ** ALSO THOSE THE ALREADY
# PRESENTS
#############################################################################

  Value_cols3 <-  c(#"Engraft_0.0" = "#ffffcc",
                     "Engraft_0.2 - 0.4" = "#c2e699",
                     "Engraft_0.4 - 0.6" = "#78c679",
                     "Engraft_0.6 - 0.8" = "#31a354", 
                     "Engraft_0.8 - 1" = "#006837",
                     "*Engraft_0.8 - 1" = "red",
                     "Present_0.8 - 1" = "#252525")

ggplot(gene_engraft_present_visual %>%
          
           #join uniprot info
          left_join(
          edited_uniprot2 %>%
          select(gene_callers_id, Protein.names, Family, Family_label, Genus) %>%
          mutate_if(is.character, trimws) %>%
          mutate(Family_label = factor(Family_label, levels = c("Lachnospiraceae",
                                                                "Ruminococcaceae",
                                                                "Other"))) %>%
          mutate(gene_callers_id = paste("g", gene_callers_id, sep = "_")),
          by = "gene_callers_id") %>%
          #mutate(Family = case_when(family.x == "NA" ~ paste("Unknown_family"),
          #                            TRUE ~ as.character(family.x))) %>%
          mutate(gene_callers_id = paste(Genus, gene_callers_id, sep = "_")) %>%
          mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)),                                          
       aes(gene_callers_id, Fig_lab)) +
  geom_tile(aes(fill= Value), alpha = 8/10) +
     scale_fill_manual(values= Value_cols3) +
   facet_grid(~Fate+Fig_lab~Family_label, scale = "free", space = "free", switch = "y") +
    theme_classic() +
     theme(#axis.title.x = element_blank(),
        text = element_text(size=10),
        legend.position="bottom",
        #legend.position="none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text= element_text(face ="bold"),
         strip.placement = "none")
 ggtitle("Common engrafted genes (minimum 3 individuals) across all PT")


```


```{r}

#############################################################################
# engrafted genes - identified proteins
#############################################################################
 
  Detection_cols3 <-  c(#"0.0" = "#ffffcc",
                     "0.2 - 0.4" = "#c2e699",
                     "0.4 - 0.6" = "#78c679",
                     "0.6 - 0.8" = "#31a354", 
                     "0.8 - 1" = "#006837")


ggplot(gene_engraft_all_visual %>%
          
         #remove whitespace
         mutate(gene_callers_id= gsub(" ", "", gene_callers_id)) %>%
         #join uniprot info
          left_join(
          edited_uniprot2 %>%
          select(gene_callers_id, Protein.names, Family, Family_label, Genus) %>%
          mutate_if(is.character, trimws) %>%
          mutate(Family_label = factor(Family_label, levels = c("Lachnospiraceae",
                                                                "Ruminococcaceae",
                                                                "Other"))) %>%
          mutate(gene_callers_id = paste("g", gene_callers_id, sep = "_")),
          by = "gene_callers_id") %>%
          #mutate(Family = case_when(family.x == "NA" ~ paste("Unknown_family"),
          #                            TRUE ~ as.character(family.x))) %>%
          mutate(gene_callers_id = paste(Genus, gene_callers_id, sep = "_")) %>%
          mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)), 
          #mutate(gene_callers_id = paste(Protein.names, gene_callers_id, sep = "_")) %>%
          #mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)), 
        aes(Fig_lab, Protein.names)) +
        geom_tile(aes(fill= Detection), alpha = 8/10) +
        scale_fill_manual(values= Detection_cols3) +
  #facet_grid(~Fate+Fig_lab~Genus, scale = "free", space = "free", switch = "y") +
        facet_grid(Family~Fate+Fig_lab, scales = "free", space = "free", switch = "x") +
        theme_bw()  +
        theme(#axis.title.x = element_blank(),
        text = element_text(size=15),
        legend.position="bottom",
        #legend.position="none",
        axis.text.y = element_text(size = 7, angle = 0, hjust = T),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 8, angle = 90),
        strip.text.y = element_text(size = 8, angle = 0),
         strip.placement = "none") +
         #coord_flip() +
 ggtitle("Common engrafted genes (minimum 3 individuals) across all PT")

  

```

```{r}
  #############################################################################
# engrafted genes - identified proteins
#############################################################################
 
  Detection_cols3 <-  c(#"0.0" = "#ffffcc",
                     "0.2 - 0.4" = "#c2e699",
                     "0.4 - 0.6" = "#78c679",
                     "0.6 - 0.8" = "#31a354", 
                     "0.8 - 1" = "#006837")


ggplot(gene_engraft_all_visual %>%
          
         #remove whitespace
         mutate(gene_callers_id= gsub(" ", "", gene_callers_id)) %>%
         #join uniprot info
          left_join(
          edited_uniprot2 %>%
          select(gene_callers_id, Protein.names, Family, Family_label, Genus) %>%
          mutate_if(is.character, trimws) %>%
          mutate(Family_label = factor(Family_label, levels = c("Lachnospiraceae",
                                                                "Ruminococcaceae",
                                                                "Other"))) %>%
          mutate(gene_callers_id = paste("g", gene_callers_id, sep = "_")),
          by = "gene_callers_id") %>%
          filter(Protein.names !="Uncharacterized protein") %>%
          filter(!is.na(Protein.names)) %>%
          #mutate(Family = case_when(family.x == "NA" ~ paste("Unknown_family"),
          #                            TRUE ~ as.character(family.x))) %>%
          mutate(gene_callers_id = paste(Genus, gene_callers_id, sep = "_")) %>%
          mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)), 
          #mutate(gene_callers_id = paste(Protein.names, gene_callers_id, sep = "_")) %>%
          #mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)), 
        aes(Fig_lab, Protein.names)) +
        geom_tile(aes(fill= Detection), alpha = 8/10) +
        scale_fill_manual(values= Detection_cols3) +
  #facet_grid(~Fate+Fig_lab~Genus, scale = "free", space = "free", switch = "y") +
        facet_grid(Family+Genus~Fate+Fig_lab, scales = "free", space = "free", switch = "x") +
        theme_bw()  +
        theme(#axis.title.x = element_blank(),
        text = element_text(size=15),
        legend.position="bottom",
        #legend.position="none",
        axis.text.y = element_text(size = 7, angle = 0, hjust = T),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 8, angle = 90),
        strip.text.y = element_text(size = 8, angle = 0),
         strip.placement = "none") +
         #coord_flip() +
 ggtitle("Common engrafted genes (minimum 3 individuals) across all PT")

  

```


```{r}
  
  #####################################################################
    Detection_cols3 <-  c(#"0.0" = "#ffffcc",
                     "0.2 - 0.4" = "#c2e699",
                     "0.4 - 0.6" = "#78c679",
                     "0.6 - 0.8" = "#31a354", 
                     "0.8 - 1" = "#006837")


gene_engraft_all_visual %>%
          
         #remove whitespace
         mutate(gene_callers_id= gsub(" ", "", gene_callers_id)) %>%
         #join uniprot info
          left_join(
          edited_uniprot3 %>%
          select(gene_callers_id, Protein.names, Family, Family_label, Genus, Pfam) %>%
          separate(Pfam, c("pfam1", "pfam2"), sep = ";") %>%
          mutate(ylab= paste(paste0("[", Protein.names, sep=""), 
                             paste0(pfam1, "]", sep=""), sep = "][")) %>%
          mutate_if(is.character, trimws) %>%
          mutate(Family_label = factor(Family_label, levels = c("Lachnospiraceae",
                                                                "Ruminococcaceae",
                                                                "Other"))) %>%
          mutate(gene_callers_id = paste("g", gene_callers_id, sep = "_")),
          by = "gene_callers_id") %>%
          filter(Protein.names !="Uncharacterized protein") %>%
          filter(!is.na(Protein.names)) %>%
          #mutate(Family = case_when(family.x == "NA" ~ paste("Unknown_family"),
          #                            TRUE ~ as.character(family.x))) %>%
          mutate(gene_callers_id = paste(Genus, gene_callers_id, sep = "_")) %>%
          mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)) %>%
          #mutate(gene_callers_id = paste(Protein.names, gene_callers_id, sep = "_")) %>%
          #mutate(gene_callers_id = gsub("_g_", "_", gene_callers_id)), 
        ggplot(aes(Fig_lab, ylab)) +
        geom_tile(aes(fill= Detection), alpha = 8/10) +
        scale_fill_manual(values= Detection_cols3) +
  #facet_grid(~Fate+Fig_lab~Genus, scale = "free", space = "free", switch = "y") +
        facet_grid(Family+Genus~Fate+Fig_lab, scales = "free", space = "free", switch = "x") +
        theme_bw()  +
        theme(#axis.title.x = element_blank(),
        text = element_text(size=15),
        legend.position="bottom",
        #legend.position="none",
        axis.text.y = element_text(size = 12, angle = 0, hjust = T),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 12, angle = 90),
        strip.text.y = element_text(size = 12, angle = 0),
         strip.placement = "none")

  

  


```











