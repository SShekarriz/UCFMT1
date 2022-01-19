---
title: "UCFM1_16S_paperAnalysis"
author: "Sharok"
date: "27/08/2019"
output: html_document
---


```{r echo=FALSE,message=FALSE}
library(dplyr)
library(tidyr)
library(ggplot2)
library(pipeR)
library(rlang)
library(phyloseq)
library(knitr)
library(kableExtra)
library(cowplot)
library(ape)
library(vegan)
library(stringi)
library(tidytext)
library(lme4)
library(lmerTest)
library(hash)


colours1 <- c("#02b6d3","#668400","#f1bc8c","#e58eff","#76b000","#ffa4e4","#bf4200","#c0019a","#c452d7","#feb873","#922d4f",
             "#763cc2","#9a3cc2","#ec0078","#c8cb82","#bc6c00","#ff7d99","#ff7bb4","#cb88a5","#63c9ff","#af003a","#b77668",
             "#9fb500","#006aaf","#ff5245","#c00067","#5fdabb","#f7bc5b","#913412","#ffa477","#009a2a","#01b97b","#de0042",
             "#109eff","#f07312","#816190","#f11067","#a11b2d","#835d00","#8857e0","#fd59d2","#a80225","#91d86a","#019d6d",
             "#eb2ea5","#62409a","#285d35","#007b19","#31e27c","#723c89","#006826","#ff8845","#ffab9a","#a2aeff","#f45a2b",
             "#ff3d54","#7b4429","#0294e5","#4343b6","#1b4eaa","#96d68b","#ff7369","#cb2d10","#74da9d","#dd7300","#843c41",
             "#81411c","#ff68b4","#4fbf36","#f3b3de","#c3cd65","#736a00","#ff73ef","#ff9e37","#8f1396","#ff9fab","#8486ff",
             "#9d175e","#957a00","#385e00","#978b56","#025eb5","#9c3100","#00714c","#5173fe","#008a77","#cc9000","#d0a400",
             "#b70072","#d10090","#86336c","#02a4ae","#63dd8f","#e1c28e","#8d3729","#1b5588","#524b85","#00b8ff","#cea7ff",
             "#ff6e74")
colours = c('grey69',"#4f8579","#783fcc","#69d24d","#cb4bbd","#c6dc46",
					"#542871","#78d792","#cc4472","#83d7d0","#d44d33","#676fcd",
					"#ceb854","#403d57","#b97839","#84a4cb","#588038","#c68ac4",
					"#48472a","#c9c39c","#6e2b34","#c78889")

# Loading JAke's package
devtools::load_all('~/Drive2/software/aftersl1p/')

```


##export data into phyloseq:
Here a table of the mapfile will be created: mapFile_dada2out.txt if the mapFile_phyloseq.txt
will not be present in the directory then mapFile_dada2out will be used to import into phyloseq

```{r}
#check the ASV table 
#asd_table <- read.csv("mergtab_nochim_UCFMT1_16SJosie_transposed.csv")

#importing MapFile
Mapfile <- read.csv("Map0vs6_Dec4_JL_editJuly17.txt", sep = "\t", stringsAsFactors = F)
#remove all the white space from the MpFile
Mapfile <- data.frame(lapply(Mapfile, trimws), stringsAsFactors = FALSE)
#edit a colname for MapFile:
colnames(Mapfile)[29] <- "MajorDonor"
# Editing Timepoint_Treatment_PatientNo variable. it was not correct
# also edit the donor variable; had to add initial + Patient No to make it 
# Uniqu
Mapfile %>%
mutate(Donor = case_when(!PatientNo == "DNA" ~ 
                             paste(Donor, PatientNo, sep = ""),
                         TRUE ~ as.character(Donor)))  %>%
mutate(Remission= case_when(Remission == "No" ~ paste("NoRes"),
                            Remission == "Yes" ~ paste("Res"),
                            TRUE ~ as.character(Remission))) %>%
mutate(Timepoint= case_when(Timepoint == "0" ~ paste("WK0"),
                            Timepoint == "6" ~ paste("WK6"),
                            TRUE ~ as.character(Timepoint))) %>%
mutate(Treatment= case_when(Treatment == "P" ~ paste("Placebo"),
                            Treatment == "A" ~ paste("FMT"),
                            TRUE ~ as.character(Treatment))) %>%
mutate(Timepoint_Treatment= case_when(Timepoint_Treatment == "A0" ~ paste("F0"),
                            Timepoint_Treatment == "A6"~ paste("F6"),
                    TRUE ~ as.character(Timepoint_Treatment))) %>%
mutate(Timepoint_Treatment_PatientNo = case_when(!Timepoint_Treatment == "DNA" ~ 
                                        paste(Timepoint_Treatment, "_pt",
                                        PatientNo, sep = ""),
                    TRUE ~ as.character(Timepoint_Treatment))) %>%
mutate(Fig_lab= case_when(!PatientNo == "DNA" ~ 
                            paste("pt", PatientNo, sep = ""),
                           PatientNo == "DNA" ~ 
                            paste(gsub("Donor", "Don", Donor)))) %>%
mutate(MajorDonor = case_when(Treatment == "Placebo" ~ paste("Placebo"),
                              TRUE ~ as.character(MajorDonor))) %>%
mutate(MajorDonor = case_when(MajorDonor == "" ~ paste("Uknown"),
                              TRUE ~ as.character(MajorDonor))) %>%
mutate_all(funs(stringr::str_replace(., "DNA", "Slurry"))) %>%
  # new variable for ordination plot
  mutate(Rem_TimeTreat= case_when(Remission == "Slurry" ~ paste("Donor"),
                                  TRUE ~ paste(Remission, Timepoint_Treatment,
                                               sep = "_"))) %>%
 mutate(TimeTreat= case_when(Fig_lab == "DonB" ~ paste(Fig_lab),
            Timepoint_Treatment != "Slurry" ~ paste(Timepoint_Treatment),
            TRUE ~ paste("Donor"))) %>%
 mutate(donorB= case_when(Fig_lab == "DonB" ~ paste(Fig_lab),
                          Fig_lab != "DonB" ~ paste("Other")))-> Mapfile


#Importing the Dada2 sequece variants file
seqtab.nochim <- readRDS("mergetab_nochim_UCFMT1_16SJosie.rds")
# Importing taxonomic assignment
taxa <- readRDS("taxa_UCFMT1_16SJosie_gg2013.rds")

# add rowname to mapfile so phyloseq can read
rownames(Mapfile) <- Mapfile$X.SampleID

#import the data into phyloseq:
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(Mapfile), 
               tax_table(taxa))
ps

############################################################
#### REMOVING SAMPLES THAT ARE NOT PAIRED ##################
############################################################
# PTs AA64 and SW68 do not have timpoint6 samples so taking them out!
# PT NS73 has only timepoint6 so taking it out as well
# 3 X SAMPLES REMOVED HERE!
sample_data(ps) <- subset_samples(ps, !Donor %in% c("AA64", "SW68", "NS73"))
############################################################
# REMOVING SAMPLES THAT HAVE BELOW 5K READS
############################################################
# HERE IS THE READ PLOT
depth_plt <- plot_read_depth(ps) +
	theme_minimal()
depth_plt
# find the samples with low read numbers anything below 10,000
#data.frame(sample_sums(ps)) %>%
#   tibble::rownames_to_column("X.SampleID") %>%
#    rename(read_length = sample_sums.ps.)) %>%
#    right_join(data.frame(sample_data(ps)), by = "X.SampleID") %>% 
#    filter(read_length <= 10000) -> low_read_numbers

########## THERE ARE THREE SAMPLES BELOW 10,000 READS:#####################
# PBH48: 212reads > A0, Non-responder, DonorB <<< REMOVED THIS FOR NOW >>>
# KPM41: 2959reads> A0, Responders,    DonorD <<< REMOVED THIS FOR NOW >>>
# CM44:  8883reads> A0, Non-responder, DonorB
# ND1:  14728reads> A6, Non-responder, DonorC
sample_data(ps) <- subset_samples(ps, !Donor %in% c("PBH48", "KPM41"))
depth_plt_trim <- plot_read_depth(ps) +
	theme_minimal() + ggtitle("Samples above 5k reads")
depth_plt_trim
# STILL SOUNDS LIKE WE HAVE TWO OUTLIER!? SHOULD I TAKE THEM OUT?
# DB16: 357838reads> A6, Non-responder, Not clear donor
###########################################################################

#asd_table <- data.frame(otu_table(ps))
kable(sample_data(ps), "html", 
      caption = "Table 1: Provided metadata within phyloseq object") %>%
  kable_styling(full_width = F, font_size = 15) %>%
  scroll_box(width = "25cm", height = "10cm") %>%
    kable_styling("striped")

# Show available ranks in the dataset
# calculate rel. abundance for ps_F0 
ps_rel <- ps %>%
transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

ps_rel

##############################################################################
# A NEW PHYLOSEQ OBJECT: ALL PTS WHO GET DONORB, Placebo AND DONORB SAMPLES
##############################################################################

ps_DonB <- subset_samples(ps, MajorDonor %in% c("B", "Placebo") | Donor == "DonorB")
ps_DonB_rel <- ps_DonB %>%
transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

ps_DonB_rel

```

##Description of UCFMT1 16S data:
Ok, I carefully looked at the mapfile and here is what I see:
We had 180 x 16s samples, but had to take out 7 samples because they where not
paired info for each patients (PTs AA64 and SW68 do not have timpoint6 
samples so taking them out. PT NS73 has only timepoint6 so taking it out as well)
or they had few reads for patients: PBH48, KPM41
From the remained 173 (16S) samples, 124 Patients and 49 Donor samples.

124 samples from Patients + 49 Donors
124 = 62(2xtimepoint) : baseline and after FMT
62 = 9(Responder) + 53(Non-responder)
62 = 31(FMT) + 31(Placebo)
31(FMT) = 7(Responder) + 24(Non-responder)
31(Placebo) = 2(Responder) + 29(Non-responder)

49 = 7(DonorA) + 34(DonorB) + DonorC(2) +
     2(DonorE) + 2(DonorF) + 2(DonorG)

Please note that this above data is not (may not) follow tha original paper
stat because we are missing some of the samples. Either because stool was not
collected or because we had to take out non- paired samples.

##Description of DonorB 16S data:
These are all the patients who receive FMT from donorB (Major donor) and all the
16S samples collected from DonorB as well as placebo samples.

40 samples from Patients + 34 DonorB + 62 placebo
40= 20(2xtimepoint) : baseline and after FMT
20= 14 Non-responders patients + 6 Responder patients
31(Placebo) = 2(Responder) + 29(Non-responder)

## Taxa stack bar charts
prep the data to Genus level and fix the taxonomic assignment for "Others" then
subset dataset to smaller chunk just for the sake of barplot visualization
```{r}
#Takes the taxon assignment from the lowest assigne
#level in an OTU's assignment and fills it in to all the lower, unresolved fiels
ps_gen_df = make_phy_df(ps_rel,'Genus', 0.001, TRUE, prop = FALSE)
ps_gen_df %>% 
  mutate(Family_lab= gsub("Other", "Taxa < 1%", Family)) -> temp

#select the top 20 families
temp %>%
  count(Family_lab, wt = Abundance) %>%
  arrange(desc(n)) %>%
  pull(Family_lab) -> tfam
  tfam20 <- tfam[1:20]
temp %>%
  mutate(Family_lab= case_when(Family_lab %in% tfam20 ~ paste(Family_lab),
                               TRUE ~ paste("Other"))) -> ps_gen_df

```


```{r}

family_col= c(
"f__Lachnospiraceae" = "#6a3d9a",
"Taxa < 1%" = "#33a02c",
"f__Ruminococcaceae" = "#1f78b4",
"f__Clostridiaceae" = "#ff7f00",
"f__Bacteroidaceae" = "#b2df8a", 
"f__Enterobacteriaceae" = "#e31a1c",
"f__Streptococcaceae" = "#fdbf6f",
"f__Bifidobacteriaceae" = "#a6cee3",
"f__Prevotellaceae" = "#ffff99",
"f__Veillonellaceae" = "#f781bf",
"f__Erysipelotrichaceae" = "#cab2d6",
"f__Coriobacteriaceae" = "#fb9a99",
"f__" = "#808080",
"Other" = "#808080",
"f__Lactobacillaceae" = "#b15928",
"f__Christensenellaceae" = "#8dd3c7",
"f__Peptostreptococcaceae" = "#fb8072",
"f__Enterococcaceae" = "#fccde5",
"f__Rikenellaceae" = "#80b1d3",
"f__[Tissierellaceae]" = "#b3de69",
"f__Verrucomicrobiaceae" = "#ccebc5"
  
)

```


#Stacked barplot of all sample (Family-level)

```{r}

###############################################################################
########################### FAMILY  ###########################################
# sum the abundance for each class, across all IDs, & sort the result
  ps_gen_df%>%
  count(Family_lab, wt = Abundance) %>%
  arrange(desc(n)) %>%
  pull(Family_lab) -> sort.family
# get ID order, sorted by each ID's abundance in the most abundant class
  ps_gen_df %>%
  filter(Family_lab == sort.family[1]) %>%
  arrange(desc(Abundance)) %>% distinct(X.SampleID) %>%
  pull(X.SampleID) -> ID.order
  
# factor ID / Class in the desired order
ps_gen_df %>%
  mutate(ID = factor(X.SampleID, levels = ID.order)) %>%
  mutate(Family_lab = factor(Family_lab, levels = rev(sort.family))) %>%
  ggplot(aes(x = ID, y = Abundance, fill = Family_lab)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = family_col) +
  scale_y_continuous(limits = c(0,1.01), expand = c(0, 0)) +
  facet_grid(~TimeTreat, space = "free", scales = "free") +
  theme_classic() + theme(legend.position = "bottom",
                          axis.text.x = element_blank(),
                          legend.title = element_blank(),
                          axis.title.x = element_blank(),
                          text = element_text(size = 10),
                          legend.text = element_text(size = 5),
                          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")) + # get rid of legend panel bg) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1))
ggsave("figs/StackBar_Family_All.png", width = 16, height = 12, units = c("cm"),
       bg = "transparent")



```


#Stacked barplot of Donor B dataset (Family-level)
```{r echo=FALSE}

# factor ID / Class in the desired order
ps_gen_df %>%
  mutate(ID = factor(X.SampleID, levels = ID.order)) %>%
  mutate(Family_lab = factor(Family_lab, levels = rev(sort.family))) %>%
  filter(MajorDonor %in% c("B", "Placebo") | Donor == "DonorB") %>%
  ggplot(aes(x = ID, y = Abundance, fill = Family_lab)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = family_col) +
  scale_y_continuous(limits = c(0,1.01), expand = c(0, 0)) +
  facet_grid(~TimeTreat, space = "free", scales = "free") +
  theme_classic() + theme(legend.position = "bottom",
                          legend.title = element_blank(),
                          axis.text.x = element_blank(),
                          axis.title.x = element_blank(),
                          text = element_text(size = 10),
                          legend.text = element_text(size = 5),
                          legend.background = element_rect(fill = "transparent"), 
    legend.box.background = element_rect(fill = "transparent")) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1))
ggsave("figs/StackBar_Family_DonB.png", width = 16, height = 8, units = c("cm"),
       bg = "transparent")

```


```{r}

## Functions to calculate alpha diversity
library(emmeans)
UCFMT1_adiv <- function(phyloOBJ) {
  
  adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(phyloOBJ, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phyloOBJ, measures = "Shannon"),
  "TimeTreat" = phyloseq::sample_data(phyloOBJ)$TimeTreat)
  adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon"))) %>%
  ggplot(aes(x = TimeTreat, y = value)) +
  geom_boxplot(outlier.color = NA) +
  facet_wrap(~metric, scales = "free", ncol = 1) +
  geom_jitter(aes(color = TimeTreat), height = 0, width = .2) +
   scale_color_manual(values = c("#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4",
                                  "#7570B3", "#252525"),
                       breaks = c("F0", "F6", "P0", "P6", "DonB", "Donor")) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none", axis.title = element_blank())
}

  UCFMT_adiv_stat <- function(phyloOBJ) {
     
#Observed    
table <- data.frame(estimate_richness(phyloOBJ, measures = "Observed"),
                            sample_data(phyloOBJ))
lmer(Observed ~ Timepoint + (1 | PatientNo),
      data = table %>% filter(TimeTreat %in% c("F0", "F6"))) -> lm_FMT_obs
FMT1 = paste("FMT, WK0 vs WK6: p=", 
      format(round(anova(lm_FMT_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lmer(Observed ~ Timepoint + (1 | PatientNo),
      data = table %>% filter(TimeTreat %in% c("P0", "P6"))) -> lm_Pla_obs
Placebo1= paste("Placebo, WK0 vs WK6: p=", 
      format(round(anova(lm_Pla_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lm(Observed ~ donorB ,
      data = table) -> lm_DonB_obs
DonB1= paste("DonB vs All other samples: p=", 
      format(round(anova(lm_DonB_obs)$"Pr(>F)"[1], 3), nsmall = 3), sep = "")

#Shannon
table <- data.frame(estimate_richness(phyloOBJ, measures = "Shannon"),
                            sample_data(phyloOBJ))
lmer(Shannon ~ Timepoint + (1 | PatientNo),
      data = table %>% filter(TimeTreat %in% c("F0", "F6"))) -> lm_FMT_obs
FMT2 = paste("FMT, WK0 vs WK6: p=", 
      format(round(anova(lm_FMT_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lmer(Shannon ~ Timepoint + (1 | PatientNo),
      data = table %>% filter(TimeTreat %in% c("P0", "P6"))) -> lm_Pla_obs
Placebo2= paste("Placebo, WK0 vs WK6: p=", 
      format(round(anova(lm_Pla_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lm(Shannon ~ donorB ,
      data = table) -> lm_DonB_obs
DonB2= paste("DonB vs All other samples: p=", 
      format(round(anova(lm_DonB_obs)$"Pr(>F)"[1], 3), nsmall = 3), sep = "")

data.frame(Observed= c(FMT1, Placebo1, DonB1),
           Shannon=  c(FMT2, Placebo2, DonB2))
  
}


```


# Shannon Measurment
Sample are rarefied to calculate these values, here I'm using two dataset;
One compose of compelete samples one ONLY donorB dataset
```{r echo=FALSE, message=FALSE}
##############################################################
####Do rarefication normalization
#Subsample reads
ps_rare <- phyloseq::rarefy_even_depth(ps, sample.size = min(sample_sums(ps)),
                                       rngseed = 123)
UCFMT1_adiv(ps_rare)
ggsave('figs/Alpha_all.png', heigh = 12, width = 6, units = c("cm"))
#print stat
UCFMT_adiv_stat(ps_rare) -> Alpha_all_stat
Alpha_all_stat

################################################################
# DonB dataset ps_DonB
ps_DonB_rare <- phyloseq::rarefy_even_depth(ps_DonB, 
                                       sample.size = min(sample_sums(ps)),
                                       rngseed = 123)
UCFMT1_adiv(ps_DonB_rare)
ggsave('figs/Alpha_DonB.png', heigh = 12, width = 6, units = c("cm"))
# print stat
UCFMT_adiv_stat(ps_DonB_rare) -> Alpha_DonB_stat
Alpha_DonB_stat

```


```{r}
# Functions for Adonis test on distances (bray-Curtis and Aitchinson)
adonis_ord <- function(phyloOBJ, ..., METHOD) {
  
sub <- subset_samples(phyloOBJ, ...)
sub_DIS <- phyloseq::distance(sub, method = METHOD)
# adonis on timepoint
set.seed(10)
adon_sub = adonis(sub_DIS ~  PatientNo + Timepoint,
                     data = data.frame(sample_data(sub)))
# extract R2 and pvalue from adonis result
R2=format(round(adon_sub$aov.tab$R2[2], 
                2), nsmall = 3)
pval=adon_sub$aov.tab$`Pr(>F)`[2]
sub_pval=paste("R2=", R2, "(p=", pval, ")", sep = "")
sub_pval

}

adonis_ord_donB <- function(phyloOBJ, METHOD) {
# adonis
set.seed(10)
donB_adonis <- adonis(phyloseq::distance(phyloOBJ, method = METHOD) ~  donorB,
                     data = data.frame(sample_data(phyloOBJ)))
# extract R2 and pvalue from adonis result
R2=format(round(donB_adonis$aov.tab$R2[1], 
                2), nsmall = 3)
pval=donB_adonis$aov.tab$`Pr(>F)`[1]
sub_pval=paste("R2=", R2, "(p=", pval, ")", sep = "")
sub_pval

}

```



```{r}
# ordination plot function for FMT

UCFMT1_ordinate <- function(phyloOBJ, ordinat, METHOD) {

## all ###################################################################
p0 <- plot_ordination(phyloOBJ, ordinat, color = 'TimeTreat') +
   scale_color_manual(values = c("#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4",
                                  "#7570B3", "#252525"),
                       breaks = c("F0", "F6", "P0", "P6", "DonB", "Donor")) +
    #geom_line(aes(group= PatientNo), color = "#bdbdbd", linetype = "dashed") +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")
# FMT ###################################################################
p <- plot_ordination(phyloOBJ, ordinat, color = 'Timepoint_Treatment') +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#bdbdbd", "#bdbdbd",
                                  "#bdbdbd"),
                       breaks = c("F0", "F6", "P0", "P6", "Slurry"), 
	                   labels = c("FMT0", "FMT6", 
	                              "Placebo0", "Placebo6", "Donor"))
p1 <- p + stat_ellipse(data= subset(p$data, 
                                    Timepoint_Treatment %in% c("F0", "F6"))) +
  annotate(geom = 'text', size=3, fontface= "bold",
           # adonis result
           label = paste("FMT:",
           adonis_ord(phyloOBJ, Treatment %in% c("FMT"), METHOD = METHOD),
           sep = " "),
           x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")

# Placebo ##############################################################
p <- plot_ordination(phyloOBJ, ordinat, color = 'Timepoint_Treatment') +
 scale_color_manual(values = c("#bdbdbd", "#bdbdbd", "#A6CEE3", "#1F78B4",
                                  "#bdbdbd"),
                       breaks = c("F0", "F6", "P0", "P6", "Slurry"), 
	                   labels = c("FMT0", "FMT6", 
	                              "Placebo0", "Placebo6", "Donor"))
p2 <- p + stat_ellipse(data= subset(p$data, 
                                    Timepoint_Treatment %in% c("P0", "P6")))+
  annotate(geom = 'text', size=3, fontface= "bold",
           # adonis result
           label = paste("Placebo:", 
           adonis_ord(phyloOBJ, Treatment %in% c("Placebo"), METHOD = METHOD),
           sep = " "),
           x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")
   

# Responders #############################################################
p <- plot_ordination(phyloOBJ, ordinat, color = 'Rem_TimeTreat') +
  scale_color_manual(
  values = c("#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd",
             "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4",
             "#bdbdbd"),
  breaks = c("NoRes_F0", "NoRes_F6", "NoRes_P0", "NoRes_P6",
             "Res_F0", "Res_F6", "Res_P0", "Res_P6",
             "Donor"))
p3 <- p + stat_ellipse(data= subset(p$data, Rem_TimeTreat %in% c("Res_F0",
                        "Res_F6"))) +
           annotate(geom = 'text', size=3, fontface= "bold",
           # adonis result
           label = paste("FMT Res:",
           adonis_ord(phyloOBJ, Remission %in% c("Res") & 
                      Treatment == "FMT", METHOD = METHOD),
           sep = " "),
           x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")
  
# NoResponders ###########################################################
p <- plot_ordination(phyloOBJ, ordinat, color = 'Rem_TimeTreat') +
  scale_color_manual(
  values = c("#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4",
             "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd",
             "#bdbdbd"),
  breaks = c("NoRes_F0", "NoRes_F6", "NoRes_P0", "NoRes_P6",
             "Res_F0", "Res_F6", "Res_P0", "Res_P6",
             "Donor"))
p4 <- p + stat_ellipse(data= subset(p$data, 
                                    Rem_TimeTreat %in% c("NoRes_F0",
                        "NoRes_F6"))) +
  annotate(geom = 'text', size=3, fontface= "bold",
           # adonis result
           label = paste("FMT NoR:",
           adonis_ord(phyloOBJ, Remission %in% c("NoRes") &
                                Treatment == "FMT", METHOD = METHOD),
           sep = " "),
           x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")
# Donor B samples #########################################################

p <- plot_ordination(phyloOBJ, ordinat, color = 'TimeTreat') +
  scale_color_manual(values = c("#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd",
                                  "#7570B3", "#bdbdbd"),
                       breaks = c("F0", "F6", "P0", "P6", "DonB", "Donor"))

p5 <- p + stat_ellipse(data= subset(p$data, TimeTreat %in% c("DonB"))) +
  annotate(geom = 'text', size=3, fontface= "bold",
           # adonis result
           label =  paste("DonB:", 
           adonis_ord_donB(phyloOBJ, METHOD = METHOD),
           sep = " "),
           x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")
  
# save_plot
plot_grid(p0, p1, p3, p5, p2, p4, 
      ncol = 3, align = "hv")
}

```


#Aitchison’s distanc
Using raw count data to transform the data to centeriod log ratio and visualize
it in via euclidian distance in PCA space.

```{r }

#CLR transform
ps_clr <- microbiome::transform(ps, "clr")
psclr_dist= phyloseq::distance(ps_clr, method = "euclidean")
psclr_ord= ordinate(ps_clr, method = "RDA", distance = psclr_dist)
# save_plot
UCFMT1_ordinate(ps_clr, psclr_ord, "euclidean")
ggsave('figs/PCOA_Aithchinson_ps_all.png', heigh = 12, width = 18, units = c("cm"))


```

# Bray-Curtis
Using relative abundance information to plot PCOA and UPCGMA Tree for All samples
and only DonorB datasets.

```{r}

bray_dis_ps = phyloseq::distance(ps_rel, method = 'bray')
bray_ord_ps = ordinate(ps_rel, method = 'PCoA', distance = bray_dis_ps)
# save_plot
UCFMT1_ordinate(ps_rel, bray_ord_ps, "bray")
ggsave('figs/PCOA_bray_ps_all.png', heigh = 12, width = 18, units = c("cm"))


```


#Aitchison’s distanc for donor B dataset
Using raw count data to transform the data to centeriod log ratio and visualize
it in via euclidian distance in PCA space.

```{r}

#CLR transform ### ONLY DONOR B SAMPLES ################
ps_clr <- microbiome::transform(ps_DonB, "clr")

psclr_dist= phyloseq::distance(ps_clr, method = "euclidean")
psclr_ord= ordinate(ps_clr, method = "PCoA", distance = psclr_dist)
# save_plot
UCFMT1_ordinate(ps_clr, psclr_ord, "euclidean")
ggsave('figs/PCOA_Aithchinson_ps_DonorB.png', heigh = 12, width = 18, units = c("cm"))

```


# Bray-Curtis for Donor B samples
Using relative abundance information to plot PCOA and UPCGMA Tree for All samples
and only DonorB datasets.

```{r}

bray_dis_ps_DonB = phyloseq::distance(ps_DonB_rel, method = 'bray')
bray_ord_ps_DonB = ordinate(ps_DonB_rel, method = 'PCoA', 
                            distance = bray_dis_ps_DonB)
# save_plot
UCFMT1_ordinate(ps_DonB_rel, bray_ord_ps_DonB, "bray")
ggsave('figs/PCOA_bray_ps_DonorB.png', heigh = 12, width = 18, units = c("cm"))

```


#Engraftment Curve for donor B (based on relative abund)
Extracting donorB ASVs and look at DonorB treated patients and placebo

```{r}

#extracting donorB's ASVs:
subset_samples(ps_DonB_rel, donorB %in% c("DonB")) -> onlyDonB
psmelt(onlyDonB) -> df_DonB

df_DonB %>%
  add_rownames("ASVn") %>%
  select(ASVn, OTU, Phylum, Family, Genus, Sample, Abundance) %>%
  spread(Sample, Abundance, fill=0) -> df_DonB_wide

# selecting ASVs that are present at least in one DonorB samples:
#df_DonB_wide[rowSums(df_DonB_wide[,6:39] > 0) >= 1, ] %>%
df_DonB_wide[rowSums(df_DonB_wide[,6:39] > 0) >= 1, ] %>%
  pull(OTU)-> core_DonB

# Finding DonB ASVs in the whole DonorB dataset
subset_taxa(ps_DonB_rel, 
            rownames(tax_table(ps_DonB_rel)) %in% core_DonB) -> core_DonB_ps
core_DonB_ps <- prop_tax_down(core_DonB_ps, TRUE, dbig = FALSE)

###############################################################################
# A function to get the commonly engrafted ASVs
###############################################################################
Find_eASV <- function(core_DonB_ps, ...) {
#############Find the engrafted ASvs in FMT
#removing DonB ASvs
subset_samples(core_DonB_ps, Fig_lab != "DonB") %>%
#Selecting patients received FMT
subset_samples(...) %>%
  psmelt() %>%
  select(Fig_lab, Timepoint, Abundance, OTU, Phylum, Family, Genus) %>%
   spread(key=Timepoint, value=Abundance) %>%
     #mutate(model= case_when(WK0 <= 0 & WK6 > 0 ~ paste("Engraft"),
     mutate(model= case_when(WK0 <= 0 & WK6 > 0.001 ~ paste("Engraft"),
                          WK0 > 0.001 & WK6 <= 0 ~ paste("Lost"),
                          WK0 > 0.001 & WK6 > 0.001 ~ paste("Present"),
                          WK0 <= 0 & WK6 <= 0 ~ paste("NotPresent"),
                            TRUE ~ paste("MinorChange"))) -> ASV_models
# Classified donor B ASVs in patients:
table(ASV_models$model)
# Engrafte ASV in WK6
ASV_models %>%
  filter(model == "Engraft") %>%
  select(-WK0) %>%
  spread(Fig_lab, WK6, fill = 0) -> ASV_engraft
# Forloop to count the number of common engraftment across patients
EngraftOUT <- list()
EngraftCount <- list()
for (i in 1:length(unique(ASV_models$Fig_lab))) {
# selecting any ASV > 0 Across all the pts
EngraftOUT[[i]] <- ASV_engraft[rowSums(ASV_engraft[,6:ncol(ASV_engraft)] > 0) >= i, ]
EngraftCount[[i]] <- data_frame(Npts = paste(i),
                                eASV = length(EngraftOUT[[i]]$OTU))
}
do.call(rbind.data.frame, EngraftCount) }
###############################################################################
########## DonorB treated
Find_eASV(core_DonB_ps, Treatment %in% c("FMT")) -> DonB_curve

####### FMT vs Placebo
# Finding DonB ASVs in the whole ps dataset
subset_taxa(ps_rel, 
            rownames(tax_table(ps_rel)) %in% core_DonB) -> core_ps
core_ps <- prop_tax_down(core_ps, TRUE, dbig = FALSE)
Find_eASV(core_ps, Treatment %in% c("FMT")) -> FMT_curve
##############################################################################
#### randomly select 20 placebo pts and find eASVs
##############################################################################
data.frame(sample_data(core_DonB_ps)) %>%
  filter(Treatment == "Placebo") %>%
  select(Fig_lab) %>% distinct() %>% pull(Fig_lab) -> placebo_pts

pt_vector <- list()
engraftRes <- list()
for (i in 1:100){
  pt_vector[[i]] <- sample(placebo_pts, 20, replace = F)
  engraftRes[[i]] <- Find_eASV(core_ps, Fig_lab %in% pt_vector[[i]])
  engraftRes[[i]]$dataset <- paste("data", i, sep = "_")
}
do.call(rbind, engraftRes) -> engraftRes_unlist

engraftRes_unlist %>%
spread(dataset, eASV) %>%
gather(dataset, Value, -Npts) %>%
group_by(Npts) %>%
summarise(
    MU = mean(Value),
    SD = sd(Value, na.rm = TRUE),
    N = sum(!is.na(Value)),
    upper_limit = MU + SD/sqrt(N),
    lower_limit = MU - SD/sqrt(N)
  ) -> Placebo_curve
##############################################################################

Placebo_curve %>%
  select(Npts, MU, SD) %>%
  # rename the avg ASV engraftment for each patient
  rename(eASV=MU) %>%
  mutate(Treatment= paste("Placebo")) %>%
  bind_rows(DonB_curve %>%
              mutate(Treatment= paste("DonorB"))) %>%
  #filter(Treatment == "Placebo") %>%
  filter(as.numeric(Npts) <= 6) %>%
  ggplot(aes(as.numeric(Npts), eASV , colour = Treatment)) +
  #facet_wrap(~Dataset, scales = "free", ncol = 1) +
  #geom_point(aes(group=Treatment))+#, size=3.5) +
  geom_line(aes(group=Treatment))+
  geom_ribbon(aes(ymin=eASV-SD,ymax=eASV+SD), linetype=2, alpha=0.1)+
  #, size=3) +
  #scale_y_sqrt() +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6)) +
  scale_color_manual(
  values = c("#bdbdbd","#006837"), 
  breaks = c("Placebo", "DonorB"), 
  labels = c("Placebo(n=20)","DonorB(n=20)")) +
  theme_bw() +
  theme(text = element_text(size = 10),
         legend.position=c(.65,.85),
        legend.title=element_blank(),
        strip.text = element_blank()) +
   xlab("Number of patients") +
   ylab("Commonly engrafted ASVs")
ggsave('figs/ASVengraft_curve.png', heigh = 8, width = 6, units = c("cm"))


```


```{r}
# finding eASV for 3 patiens:
#removing DonB ASvs
subset_samples(core_DonB_ps, Fig_lab != "DonB") %>%
  #Selecting patients received FMT
 subset_samples(Treatment %in% c("FMT")) %>%
  psmelt() %>%
  select(Fig_lab, Timepoint, Abundance, OTU, Phylum, Family, Genus) %>%
   spread(key=Timepoint, value=Abundance) %>%
     #mutate(model= case_when(WK0 <= 0 & WK6 > 0 ~ paste("Engraft"),
     mutate(model= case_when(WK0 <= 0 & WK6 > 0.001 ~ paste("Engraft"),
                          WK0 > 0.001 & WK6 <= 0 ~ paste("Lost"),
                          WK0 > 0.001 & WK6 > 0.001 ~ paste("Present"),
                          WK0 <= 0 & WK6 <= 0 ~ paste("NotPresent"),
                            TRUE ~ paste("MinorChange"))) -> ASV_models
# Classified donor B ASVs in patients:
table(ASV_models$model)
# Engrafte ASV in WK6
ASV_models %>%
  filter(model == "Engraft") %>%
  select(-WK0) %>%
  spread(Fig_lab, WK6, fill = 0) -> ASV_engraft
# Forloop to count the number of common engraftment across patients
EngraftOUT <- list()
EngraftCount <- list()
for (i in 1:length(unique(ASV_models$Fig_lab))) {
# selecting any ASV > 0 Across all the pts
EngraftOUT[[i]] <- ASV_engraft[rowSums(ASV_engraft[,6:ncol(ASV_engraft)] > 0) >= i, ]
EngraftCount[[i]] <- data_frame(Npts = paste(i),
                                eASV = length(EngraftOUT[[i]]$OTU))
}

data.frame(EngraftOUT[3]) -> eASV_3pt
colnames(eASV_3pt)
eASV_3pt %>%
  gather(Fig_lab, Abundance, -OTU, -Phylum, -Family, -Genus, -model, na.rm = TRUE) %>%
  mutate(model= case_when(Abundance >= 0.001 ~ paste("Engraft"),
                          Abundance <= 0.001 ~ paste("Other"))) %>%
  left_join(
  data.frame(sample_data(ps_DonB_rel)) %>%
  filter(Treatment == "FMT") %>%
  select(Fig_lab, Remission) %>% distinct(), by = "Fig_lab") %>%
  left_join(
    eASV_3pt %>% select(OTU) %>% distinct() %>%
    add_rownames("Num"), by = "OTU") %>%
  mutate(ASV=paste(Num, Family, sep = "_")) %>%
  ggplot(aes(ASV, Fig_lab, fill=model)) +
  facet_grid(Remission~Family, scales = "free", space = "free") +
  scale_fill_manual(values = c("#006d2c", "NA"),
                    breaks = c("Engraft", "Other")) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 5),
        strip.text.x = element_blank(),
        legend.position = "none") + xlab("ASV") + ylab("Donor B treated patients") 

ggsave('figs/eASV_3pt.png', heigh = 10, width = 10, units = c("cm"))

```


#UPGMA Tree
Making UPGMA tree for all dataset and donorB dataset. Also creating iTOL mapfile
information to plot Timepoint, Treatment, Remission and Patients ids. Also,
visualize the barplot on the clustering trees.

```{r echo=FALSE, message=FALSE}
# UPGMA Tree for all PTS
bray_upgma_ps = as.phylo(hclust(bray_dis_ps, method = 'average'))
write.tree(phy = bray_upgma_ps, file = 'ITOL/bray_upgma_allPT.nwk')
# UPGMA Tree for only DonorB
bray_upgma_ps_DonB = as.phylo(hclust(bray_dis_ps_DonB, method = 'average'))
write.tree(phy = bray_upgma_ps_DonB, file = 'ITOL/bray_upgma_DonB.nwk')

```


Make the iTOL data files
```{r, echo=FALSE, message=FALSE}
# To make the taxa barcharts I need to glom and melt to phylum
ps_phyl_df = make_phy_df(ps_rel, rank = 'Phylum', cutoff = 0.01, prop = FALSE)
save(ps_phyl_df, file = 'ITOL/AllPTs_TaxaObjs_phyl_df.RData')
unique(ps_phyl_df$Phylum)

ps_DonB_phyl_df = make_phy_df(ps_DonB_rel, rank = 'Phylum', cutoff = 0.01, prop = FALSE)
save(ps_DonB_phyl_df, file = 'ITOL/DonB_PTs_TaxaObjs_phyl_df.RData')
unique(ps_DonB_phyl_df$Phylum)

```

Making iTOL map file for all the samples:

```{r echo=FALSE, message=FALSE}

Time_Treat_color <- read.csv("colours/Timepoint_Treatment.txt",
                             sep = "\t")
Remission_color <- read.csv("colours/Remission.txt",
                             sep = "\t")
patient_color <- read.csv("colours/Patients.txt", sep = "\t")
###############################################################################
# Generate the sample labels file: 'Donor_Timepoint'
###############################################################################
samps = sample_data(ps_rel)
# Generate the tip label file
labs = data.frame(SampleID = samps$X.SampleID, 
				  #Label = paste(samps$Donor, 
				  #              samps$Timepoint, sep = '_'))
          Label = samps$Fig_lab)
ID_lab_head = c('LABELS',
				'SEPARATOR TAB',
				'DATASET_LABEL,IDc_labels',
				'DATA')
id_f = file('ITOL/AllPTs_ID_labs.txt')
writeLines(ID_lab_head, id_f)
close(id_f)
id_f = file('ITOL/AllPTs_ID_labs.txt', open = 'at')
write.table(labs, col.names = FALSE, file = id_f,
			sep = '\t', row.names = FALSE, qmethod = NULL, quote = FALSE)
close(id_f)
###############################################################################
# Generate the Timepoint_Treatment colouring file
###############################################################################
Time_Treat = data.frame(SampleID = samps$X.SampleID,
				  TT = samps$Timepoint_Treatment)
Time_Treat %>% 
  mutate(TT = case_when(TT == "Slurry" ~ paste("Donor"),
                        TRUE ~ as.character(TT))) -> Time_Treat

Time_Treat %>% left_join(Time_Treat_color, by = "TT") %>%
	select(-TT) -> Time_Treat_cols
TT_head = c('DATASET_COLORSTRIP',	
			'SEPARATOR TAB',
			'DATASET_LABEL\tTimepoint/Treatment',
			'COLOR_BRANCHES\t0',	
			'LEGEND_TITLE\tTimepoint/Treatment',
			'LEGEND_SHAPES\t1\t1\t1\t1\t1',
		    'LEGEND_COLORS\t#B2DF8A\t#33A02C\t#A6CEE3\t#1F78B4\t#7570B3',
		    'LEGEND_LABELS\tFMT0\tFMT6\tPlacebo0\tPlacebo6\tDonor',
		    'DATA')
TT_f = file('ITOL/AllPTs_Timepoint_Treatment.txt')
writeLines(TT_head, TT_f)
close(TT_f)
TT_f = file('ITOL/AllPTs_Timepoint_Treatment.txt', open = 'at')
write.table(Time_Treat_cols, file = TT_f, sep='\t', row.names = FALSE, 
            col.names = FALSE,
			qmethod = NULL, quote = FALSE)
close(TT_f)
###############################################################################
# Generate the Remission colouring file
###############################################################################
Remission = data.frame(SampleID = samps$X.SampleID,
				  TT = samps$Remission)
Remission %>% 
  mutate(TT = case_when(TT == "Slurry" ~ paste("Donor"),
                        TRUE ~ as.character(TT))) -> Remission

Remission %>% left_join(Remission_color, by = "TT") %>%
	select(-TT) -> Remission_cols
TT_head = c('DATASET_COLORSTRIP',	
			'SEPARATOR TAB',
			'DATASET_LABEL\tRemission',
			'COLOR_BRANCHES\t0',	
			'LEGEND_TITLE\tRemission',
			'LEGEND_SHAPES\t1\t1\t1',
		    'LEGEND_COLORS\t#00441b\t#E7298A\t#7570B3',
		    'LEGEND_LABELS\tResponder\tNon-responder\tDonor',
		    'DATA')
TT_f = file('ITOL/AllPTs_Remission.txt')
writeLines(TT_head, TT_f)
close(TT_f)
TT_f = file('ITOL/AllPTs_Remission.txt', open = 'at')
write.table(Remission_cols, file = TT_f, sep='\t', row.names = FALSE, 
            col.names = FALSE,
			qmethod = NULL, quote = FALSE)
close(TT_f)
################################################################################
## Generate the patient colour file
################################################################################
pats = samps[,c('X.SampleID','Donor')]
pats = data.frame(SampleID= samps$X.SampleID,
                  Donor = samps$Donor)
pats %>% left_join(patient_color, by = "Donor") %>%
	select(-Donor) %>%
	mutate(clade = 'clade',normal = 'normal') %>%
	select(SampleID, clade, colour, normal)-> pat
pat_head = c('TREE_COLORS',
			 'SEPARATOR TAB',
			 'DATA')
pat_f = file('ITOL/AllPTs_Donor_cladeColor.txt')
writeLines(pat_head, pat_f)
close(pat_f)
pat_f = file('ITOL/AllPTs_Donor_cladeColor.txt', open = 'at')
write.table(pat, file = pat_f, sep = '\t', row.names = FALSE, col.names = FALSE,
			qmethod = NULL, quote = FALSE)
close(pat_f)
################################################################################
## Generate the taxa bar chart file
################################################################################

phy_col <- data.frame(Phylum= rev(levels(ps_phyl_df$Phylum)),
           colour= rev(colours[1:9]))

ps_phyl_df %>%
	select(X.SampleID, Phylum, Abundance) %>%
	spread(Phylum, Abundance, fill = 0) -> ps_phyl_itol
ps_phyl_itol[,2:ncol(ps_phyl_itol)] = ps_phyl_itol[,ncol(ps_phyl_itol):2]
names(ps_phyl_itol) = names(ps_phyl_itol)[c(1,ncol(ps_phyl_itol):2)]
ps_phy_names = rev(levels(ps_phyl_df$Phylum))
ps_phy_here = as.character(phy_col$colour)
ps_phy_head = c('DATASET_MULTIBAR',
			 'SEPARATOR TAB',
			 'DATASET_LABEL\tPhyla',
			 stri_join(c('FIELD_LABELS', ps_phy_names), collapse = '\t'),
			 stri_join(c('FIELD_COLORS',ps_phy_here), collapse = '\t'),
			 'LEGEND_TITLE\tAbundant Phyla',
			 stri_join(c('LEGEND_SHAPES',rep('1',length(ps_phy_names))), 
			 		  collapse = '\t'),
			 stri_join(c('LEGEND_COLORS',ps_phy_here), 
			 		  collapse = '\t'),
			 stri_join(c('LEGEND_LABELS', ps_phy_names), 
			 		  collapse = '\t'),
			 'DATA')
ps_phy_f = file('ITOL/AllPTs_taxa_bar_phylum.txt')
writeLines(ps_phy_head,ps_phy_f)
close(ps_phy_f)
ps_phy_f = file('ITOL/AllPTs_taxa_bar_phylum.txt', open = 'at')
write.table(ps_phyl_itol, file = ps_phy_f, sep = '\t', 
			row.names = FALSE, col.names = FALSE, qmethod = NULL, quote = FALSE)
close(ps_phy_f)
```

Making iTOL mapfiles for only DonorB samples:

```{r echo=FALSE, message=FALSE}

###############################################################################
# Generate the sample labels file: 'Donor_Timepoint'
###############################################################################
samps_DonB = sample_data(ps_DonB_rel)
# Generate the tip label file
labs = data.frame(SampleID = samps_DonB$X.SampleID, 
				  #Label = paste(samps_DonB$Donor, 
				  #              samps_DonB$Timepoint, sep = '_'))
          Label = samps_DonB$Fig_lab)
ID_lab_head = c('LABELS',
				'SEPARATOR TAB',
				'DATASET_LABEL,IDc_labels',
				'DATA')
id_f = file('ITOL/DonBPTs_ID_labs.txt')
writeLines(ID_lab_head, id_f)
close(id_f)
id_f = file('ITOL/DonBPTs_ID_labs.txt', open = 'at')
write.table(labs, col.names = FALSE, file = id_f,
			sep = '\t', row.names = FALSE, qmethod = NULL, quote = FALSE)
close(id_f)

###############################################################################
# Generate the Timepoint_Treatment colouring file
###############################################################################
Time_Treat = data.frame(SampleID = samps_DonB$X.SampleID,
				  TT = samps_DonB$Timepoint_Treatment)

Time_Treat %>% 
  mutate(TT = case_when(TT == "Slurry" ~ paste("Donor"),
                        TRUE ~ as.character(TT))) -> Time_Treat


Time_Treat %>% left_join(Time_Treat_color, by = "TT") %>%
	select(-TT) -> Time_Treat_cols
TT_head = c('DATASET_COLORSTRIP',	
			'SEPARATOR TAB',
			'DATASET_LABEL\tTimepoint/Treatment',
			'COLOR_BRANCHES\t0',	
			'LEGEND_TITLE\tTimepoint/Treatment',
			'LEGEND_SHAPES\t1\t1\t1\t1\t1',
		    'LEGEND_COLORS\t#B2DF8A\t#33A02C\t#A6CEE3\t#1F78B4\t#7570B3',
		    'LEGEND_LABELS\tFMT0\tFMT6\tPlacebo0\tPlacebo6\tDonor',
		    'DATA')
TT_f = file('ITOL/DonBPTs_Timepoint_Treatment.txt')
writeLines(TT_head, TT_f)
close(TT_f)
TT_f = file('ITOL/DonBPTs_Timepoint_Treatment.txt', open = 'at')
write.table(Time_Treat_cols, file = TT_f, sep='\t', row.names = FALSE, 
            col.names = FALSE,
			qmethod = NULL, quote = FALSE)
close(TT_f)
###############################################################################
# Generate the Remission colouring file
###############################################################################
Remission = data.frame(SampleID = samps_DonB$X.SampleID,
				  TT = samps_DonB$Remission)
Remission %>% 
  mutate(TT = case_when(TT == "Slurry" ~ paste("Donor"),
                        TRUE ~ as.character(TT))) -> Remission

Remission %>% left_join(Remission_color, by = "TT") %>%
	select(-TT) -> Remission_cols
TT_head = c('DATASET_COLORSTRIP',	
			'SEPARATOR TAB',
			'DATASET_LABEL\tRemission',
			'COLOR_BRANCHES\t0',	
			'LEGEND_TITLE\tRemission',
			'LEGEND_SHAPES\t1\t1\t1',
		    'LEGEND_COLORS\t#00441b\t#E7298A\t#7570B3',
		    'LEGEND_LABELS\tResponder\tNon-responder\tDonor',
		    'DATA')
TT_f = file('ITOL/DonBPTs_Remission.txt')
writeLines(TT_head, TT_f)
close(TT_f)
TT_f = file('ITOL/DonBPTs_Remission.txt', open = 'at')
write.table(Remission_cols, file = TT_f, sep='\t', row.names = FALSE, 
            col.names = FALSE,
			qmethod = NULL, quote = FALSE)
close(TT_f)
################################################################################
## Generate the patient colour file
################################################################################
pats_DonB = samps_DonB[,c('X.SampleID','Donor')]
pats_DonB = data.frame(SampleID= samps_DonB$X.SampleID,
                  Donor = samps_DonB$Donor)
pats_DonB %>% left_join(patient_color, by = "Donor") %>%
	select(-Donor) %>%
	mutate(clade = 'clade',normal = 'normal') %>%
	select(SampleID, clade, colour, normal)-> pat_DonB
pat_head = c('TREE_COLORS',
			 'SEPARATOR TAB',
			 'DATA')
pat_f_DonB = file('ITOL/DonBPTs_Donor_cladeColor.txt')
writeLines(pat_head, pat_f_DonB)
close(pat_f_DonB)
pat_f_DonB = file('ITOL/DonBPTs_Donor_cladeColor.txt', open = 'at')
write.table(pat_DonB, file = pat_f_DonB, sep = '\t', row.names = FALSE, 
            col.names = FALSE,
			qmethod = NULL, quote = FALSE)
close(pat_f_DonB)
################################################################################
## Generate the donor B color file for tree clades
################################################################################
pats_DonB = samps_DonB[,c('X.SampleID','Donor')]
pats_DonB = data.frame(SampleID= samps_DonB$X.SampleID,
                  Donor = samps_DonB$Donor)
pats_DonB %>% left_join(patient_color, by = "Donor") %>%
	select(-Donor) %>%
	mutate(clade = 'clade',normal = 'normal') %>%
	select(SampleID, clade, colour, normal)-> pat_DonB
pat_head = c('TREE_COLORS',
			 'SEPARATOR TAB',
			 'DATA')
pat_f_DonB = file('ITOL/DonBPTs_Donor_cladeColor.txt')
writeLines(pat_head, pat_f_DonB)
close(pat_f_DonB)
pat_f_DonB = file('ITOL/DonBPTs_Donor_cladeColor.txt', open = 'at')
write.table(pat_DonB, file = pat_f_DonB, sep = '\t', row.names = FALSE, 
            col.names = FALSE,
			qmethod = NULL, quote = FALSE)
close(pat_f_DonB)

################################################################################
## Generate the taxa bar chart file
################################################################################

DonB_phy_col <- phy_col[-7, ]
DonB_phy_col <- phy_col

ps_DonB_phyl_df %>%
	select(X.SampleID, Phylum, Abundance) %>%
	spread(Phylum, Abundance, fill = 0) -> ps_DonB_phyl_itol
ps_DonB_phyl_itol[,2:ncol(ps_DonB_phyl_itol)] = ps_DonB_phyl_itol[,ncol(ps_DonB_phyl_itol):2]
names(ps_DonB_phyl_itol) = names(ps_DonB_phyl_itol)[c(1,ncol(ps_DonB_phyl_itol):2)]
ps_DonB_phy_names = rev(levels(ps_DonB_phyl_df$Phylum))
ps_DonB_phy_here = as.character(DonB_phy_col$colour)
ps_DonB_phy_head = c('DATASET_MULTIBAR',
			 'SEPARATOR TAB',
			 'DATASET_LABEL\tPhyla',
			 stri_join(c('FIELD_LABELS', ps_DonB_phy_names), collapse = '\t'),
			 stri_join(c('FIELD_COLORS',ps_DonB_phy_here), collapse = '\t'),
			 'LEGEND_TITLE\tAbundant Phyla',
			 stri_join(c('LEGEND_SHAPES',rep('1',length(ps_DonB_phy_names))), 
			 		  collapse = '\t'),
			 stri_join(c('LEGEND_COLORS',ps_DonB_phy_here), 
			 		  collapse = '\t'),
			 stri_join(c('LEGEND_LABELS', ps_DonB_phy_names),
			 		  collapse = '\t'),
			 'DATA')
ps_DonB_phy_f = file('ITOL/DonBPTs_taxa_bar_phylum.txt')
writeLines(ps_DonB_phy_head, ps_DonB_phy_f)
close(ps_DonB_phy_f)
ps_DonB_phy_f = file('ITOL/DonBPTs_taxa_bar_phylum.txt', open = 'at')
write.table(ps_DonB_phyl_itol, file = ps_DonB_phy_f, sep = '\t', 
			row.names = FALSE, col.names = FALSE, qmethod = NULL, quote = FALSE)
close(ps_DonB_phy_f)

```


#temporary for mike's grant
```{r}
mike_ordinate <- function(phyloOBJ, ordinat, METHOD) {

## all ###################################################################
p0 <- plot_ordination(phyloOBJ, ordinat, color = 'TimeTreat') +
   scale_color_manual(values = c("#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4",
                                  "#7570B3", "#252525"),
                       breaks = c("F0", "F6", "P0", "P6", "DonB", "Donor")) +
    geom_line(aes(group= PatientNo), color = "#bdbdbd", linetype = "dashed") +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")
# FMT ###################################################################
p <- plot_ordination(phyloOBJ, ordinat, color = 'Timepoint_Treatment') +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#bdbdbd", "#bdbdbd",
                                  "#bdbdbd"),
                       breaks = c("F0", "F6", "P0", "P6", "Slurry"), 
	                   labels = c("FMT0", "FMT6", 
	                              "Placebo0", "Placebo6", "Donor"))
p1 <- p + stat_ellipse(data= subset(p$data, 
                                    Timepoint_Treatment %in% c("F0", "F6"))) +
  annotate(geom = 'text', size=3, fontface= "bold",
           # adonis result
           label = paste("FMT:",
           adonis_ord(phyloOBJ, Treatment %in% c("FMT"), METHOD = METHOD),
           sep = " "),
           x = -Inf, y = Inf, hjust = 0, vjust = 1) +
     geom_line(aes(group= PatientNo), color = "#bdbdbd", linetype = "dashed") +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")

plot_grid(p0, p1, 
      ncol = 2, align = "hv")
}
```


```{r}
metagenome <- read.csv("map_of_metagenomics.txt", sep = "\t")
metagenome %>% filter(!Fig_lab=="DonB") %>% pull(Fig_lab)-> meta_pts
unique(meta_pts) -> meta_pts

METAps <- subset_samples(ps_DonB, Fig_lab %in% meta_pts)

METAps_rel <- METAps %>%
transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

bray_dis_ps_DonB = phyloseq::distance(METAps_rel, method = 'bray')
bray_ord_ps_DonB = ordinate(METAps_rel, method = 'PCoA', 
                            distance = bray_dis_ps_DonB)
# save_plot
mike_ordinate(METAps_rel, bray_ord_ps_DonB, "bray")
ggsave('~/Drive2/UC_FMT/UCFMT1_ShotgunReads/temp/PcoA_bray_16s.png', 
       heigh = 6, width = 12, units = "cm")

#CLR transform ### ONLY DONOR B SAMPLES ################
METAps_clr <- microbiome::transform(METAps, "clr")

psclr_dist= phyloseq::distance(METAps_clr, method = "euclidean")
psclr_ord= ordinate(METAps_clr, method = "PCoA", distance = psclr_dist)
# save_plot
mike_ordinate(METAps_clr, psclr_ord, "euclidean")
ggsave('~/Drive2/UC_FMT/UCFMT1_ShotgunReads/temp/RDA_Ait_16s.png', 
       heigh = 6, width = 12, units = "cm")

```


```{r}

Metaps_taxa <- subset_samples(ps_DonB, Fig_lab %in% meta_pts)
(rare <- phyloseq::rarefy_even_depth(Metaps_taxa, 
                          sample.size = min(sample_sums(Metaps_taxa)),
                                       rngseed = 123))

adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(rare, measures = "Shannon"),
  "TimeTreat" = phyloseq::sample_data(rare)$Timepoint)
  adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon"))) %>%
  ggplot(aes(x = TimeTreat, y = value)) +
  geom_boxplot(outlier.color = NA) +
  facet_wrap(~metric, scales = "free") +
  geom_jitter(aes(color = TimeTreat), height = 0, width = .2) +
   scale_color_manual(
     values = c("#238b45", "#005824", "#99d8c9", "#66c2a4",
                "#A6CEE3", "#1F78B4",
                "#7570B3"),
     breaks = c("FR0", "FR6", "FN0", "FN6", 
                "PP0", "PP6",
                "DonB")) +
    ggplot2::theme_bw() +
    theme(legend.title = element_blank(),
        text = element_text(size = 20),
        legend.position = "top", axis.title = element_blank())
  ggsave('~/Drive2/UC_FMT/UCFMT1_ShotgunReads/temp/alpha_taxa_16s.png', 
       heigh = 6, width = 12, units = "cm")



```




