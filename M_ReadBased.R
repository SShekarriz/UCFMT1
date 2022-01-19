---
title: "UCFMT1_ShotgunRead_analysis"
author: "Sharok"
date: "27/08/2019"
output: html_document
---

Here are my analysis to visualize the Humann2 and metaphlan2 outputs. I'll 
convert metaphlan2 and pathabundance file into phyloseq object and calculate
bray-curtis and visualize it in a tree. I'll use gene family file to find
engraftment by using DonorB 2013 sample.

To make a phyloseq object the taxonomy/function information must be separated
from the count/normalized data make them both as a matrix and import them
both seprately as OTU table and taxonomy table into phyloseq.

```{r echo=FALSE, message=FALSE}
library(tidyverse)
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

# Loading JAke's package
devtools::load_all('~/Drive2/software/aftersl1p/')

```

Importing all the Humann2/metaphlan2 outputs as well as the mapfile

```{r}

genefamilies <- read.csv("Humann2/UCFMT1_ShotgunRead_genefamilies_norm.txt",  
                         sep = "\t", stringsAsFactors = FALSE)

#pathcoverage <- read.csv("UCFMT1_ShotgunRead_pathcoverage.txt", 
#                         sep = "\t", stringsAsFactors = FALSE)

#pathabundance <- read.csv("UCFMT1_ShotgunRead_pathabundance_norm.txt", 
#                          sep = "\t", stringsAsFactors = FALSE)

metaphlan2 <- read.csv("Humann2/UCFMT1_ShotgunRead_metaphlan2.txt", 
                       sep = "\t", stringsAsFactors = FALSE)

```


```{r}
# edit the sample names and proper variable for merging to mapfile
metaphlan2 %>%
  rename_all(funs(str_replace(., "_clean_metaphlan_bugs_list", ""))
  ) %>%
  rename_all(funs(str_replace(., "X", "s"))
  ) -> metaphlan2
colnames(metaphlan2)[1] <- "SampleID"

# collapse the data at the species-level
metaphlan2 %>%
  filter(grepl("s__", SampleID)) %>%
  filter(!grepl("t__", SampleID)) -> metaphlan2

# making a data.frame for visulaization
metaphlan2 %>%
  filter(grepl("s__", SampleID)) %>%
  filter(!grepl("t__", SampleID)) %>%
  gather(SampleName, Abundance, -SampleID) %>%
  separate(SampleID, c("Kingdom", "Phylum", "Class", "Order", "Family", 
                         "Genus", "Species"),sep = "[|]") %>%
  #select(-Class) %>%
  left_join(mapfile, by = "SampleName")-> metaphlan2_df

metaphlan2 %>%
  rownames_to_column("OTU") %>%
  mutate(OTU = paste("OTU", OTU, sep = "_")) %>%
  column_to_rownames("OTU") %>%
  select(-SampleID) -> metaphlan_otu

metaphlan2 %>%
  rownames_to_column("OTU") %>%
  mutate(OTU = paste("OTU", OTU, sep = "_")) %>%
  column_to_rownames("OTU") %>%
  select(SampleID) %>%
  separate(SampleID, c("Kingdom", "Phylum", "Class", "Order", "Family", 
                         "Genus", "Species"),sep = "[|]") %>%
  select(-Class) -> metaphlan_taxa

# creat a phyloseq object contain OTU, Taxa, and mapfile
metaphlan_TAX = tax_table(as.matrix(metaphlan_taxa))
#metaphlan_OTU <- otu_table(as.matrix(metaphlan_otu), taxa_are_rows = TRUE)
metaphlan_OTU <- otu_table(as.matrix(metaphlan_otu), taxa_are_rows = TRUE)
ps_metaphlan = phyloseq(metaphlan_OTU, metaphlan_TAX)
sample_data(ps_metaphlan) <- mapfile


###############################################################################
###############################################################################

# metaphlan output is relative to 100, here convert this to realtive to 
ps_metaphlan_rel <- ps_metaphlan %>%
transform_sample_counts(function(x) {x / 100} )

ps_prop <- prop_tax_down(ps_metaphlan_rel, TRUE)
ps_phylum_df = make_phy_df(ps_prop,'Phylum', 0.01, TRUE, prop = FALSE)
ps_phyl_df = make_phy_df(ps_prop,'Genus', 0.01, TRUE, prop = FALSE)

Phy_donorB <- ggplot(ps_phyl_df,
                 aes(x = X.SampleID, y = Abundance, fill = Family)) + 
  #facet_grid(~Donor, space = "free_x", scales = "free") +
  facet_grid(~Status+ Fig_lab, scales = "free", space = "free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colours1) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        text = element_text(size=20),
        legend.position="top",
        #legend.position="none",
        axis.title.y = element_blank(), 
        #axis.text.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 90),
        strip.text.x = element_text(size = 20)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1))
pdf('figs/StackBar_Fam.pdf', heigh = 8, width = 15)
Phy_donorB
dev.off()
Phy_donorB


```

## Functions to calculate alpha diversity
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
}

  UCFMT_adiv_stat <- function(phyloOBJ) {
     
#Observed    
table <- data.frame(estimate_richness(phyloOBJ, measures = "Observed"),
                            sample_data(phyloOBJ))
lmer(Observed ~ Timepoint + (1 | Trial_no),
      data = table %>% filter(Timepoint %in% c("WK0", "WK6") &
                                Treatment == "FMT")) -> lm_FMT_obs
FMT1 = paste("FMT, WK0 vs WK6: p=", 
      format(round(anova(lm_FMT_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lmer(Observed ~ Timepoint + (1 | Trial_no),
      data = table %>% filter(TimeTreat %in% c("FR0", "FR6"))) -> lm_Res_obs
Res1= paste("FMT Responder, WK0 vs WK6: p=", 
      format(round(anova(lm_Res_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lmer(Observed ~ Timepoint + (1 | Trial_no),
      data = table %>% filter(TimeTreat %in% c("FN0", "FN6"))) -> lm_NoRes_obs
NoRes1= paste("FMT NoResponder, WK0 vs WK6: p=", 
      format(round(anova(lm_NoRes_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lm(Observed ~ donorB,
      data = table) -> lm_DonB_obs
DonB1= paste("DonB vs All other samples: p=", 
      format(round(anova(lm_DonB_obs)$"Pr(>F)"[1], 3), nsmall = 3), sep = "")

#Shannon
table <- data.frame(estimate_richness(phyloOBJ, measures = "Shannon"),
                            sample_data(phyloOBJ))
FMT2 = paste("FMT, WK0 vs WK6: p=", 
      format(round(anova(lm_FMT_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lmer(Shannon ~ Timepoint + (1 | Trial_no),
      data = table %>% filter(TimeTreat %in% c("FR0", "FR6"))) -> lm_Res_obs
Res2= paste("FMT Responder, WK0 vs WK6: p=", 
      format(round(anova(lm_Res_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lmer(Shannon ~ Timepoint + (1 | Trial_no),
      data = table %>% filter(TimeTreat %in% c("FN0", "FN6"))) -> lm_NoRes_obs
NoRes2= paste("FMT NoResponder, WK0 vs WK6: p=", 
      format(round(anova(lm_NoRes_obs)$"Pr(>F)", 3), nsmall = 3), sep = "")
lm(Shannon ~ donorB,
      data = table) -> lm_DonB_obs
DonB2= paste("DonB vs All other samples: p=", 
      format(round(anova(lm_DonB_obs)$"Pr(>F)"[1], 3), nsmall = 3), sep = "")

data.frame(Observed= c(FMT1, Res1, NoRes1, DonB1),
           Shannon=  c(FMT2, Res2, NoRes2, DonB2))
}

```


# Shannon Measurment
Sample are rarefied to calculate these values
```{r}
##############################################################
####Do rarefication normalization
#Subsample reads
(ps_metaphalan_rare <- phyloseq::rarefy_even_depth(ps_metaphlan, 
                          sample.size = min(sample_sums(ps_metaphlan)),
                                       rngseed = 123))
UCFMT1_adiv(ps_metaphalan_rare)
ggsave('figs/metaphlan_Alpha_all.pdf.pdf', heigh = 5, width = 10)

UCFMT_adiv_stat(ps_metaphalan_rare) -> Alpha_all_stat
Alpha_all_stat

```


```{r}
# Functions for Adonis test on distances (bray-Curtis and Aitchinson)
adonis_ord <- function(phyloOBJ, ..., METHOD) {
  
sub <- subset_samples(phyloOBJ, ...)
sub_DIS <- phyloseq::distance(sub, method = METHOD)
# adonis on timepoint
set.seed(10)
adon_sub = adonis(sub_DIS ~  Trial_no + Timepoint,
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

# all ###################################################################
p0 <- plot_ordination(phyloOBJ, ordinat, color = 'TimeTreat') +
  scale_color_manual(
     values = c("#B2DF8A", "#33A02C", "#B2DF8A", "#33A02C",
                "#A6CEE3", "#1F78B4",
                "#7570B3"),
     breaks = c("FR0", "FR6", "FN0", "FN6", 
                "PP0", "PP6",
                "DonB")) +
    geom_line(aes(group= Initial), color = "#bdbdbd", linetype = "dashed") +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")

# FMT ###################################################################
p <- plot_ordination(phyloOBJ, ordinat, color = 'Timepoint_Treatment') +
  scale_color_manual(
     values = c("#B2DF8A", "#33A02C", "#bdbdbd", "#bdbdbd",
                "#bdbdbd"),
     breaks = c("F0", "F6", "P0", "P6", 
                "DonB")) +
    geom_line(aes(group= Initial), color = "#bdbdbd", linetype = "dashed") 

p1 <- p + stat_ellipse(data= subset(p$data, Timepoint_Treatment %in% c("F0", "F6")))+
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
  

# Donor B samples #########################################################

p <- plot_ordination(phyloOBJ, ordinat, color = 'TimeTreat') +
    scale_color_manual(
     values = c("#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd",
                "#bdbdbd", "#bdbdbd",
                "#7570B3"),
     breaks = c("FR0", "FR6", "FN0", "FN6", 
                "PP0", "PP6",
                "DonB")) +
  geom_line(aes(group= Initial), color = "#bdbdbd", linetype = "dashed")

p5 <- p + stat_ellipse(data= subset(p$data, TimeTreat %in% c("DonB")))+
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

plot_grid(p0, p1, p5, 
      ncol = 3, align = "hv")

}

```


```{r}

# ordination plot function for FMT

UCFMT1_ordinate2 <- function(phyloOBJ, ordinat, METHOD) {
  
# all ###################################################################
p0 <- plot_ordination(phyloOBJ, ordinat, color = 'TimeTreat') +
  scale_color_manual(
     values = c("#B2DF8A", "#33A02C", "#B2DF8A", "#33A02C",
                "#A6CEE3", "#1F78B4",
                "#7570B3"),
     breaks = c("FR0", "FR6", "FN0", "FN6", 
                "PP0", "PP6",
                "DonB")) +
    geom_line(aes(group= Initial), color = "#bdbdbd", linetype = "dashed") +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")  

# Responders #############################################################
p <- plot_ordination(phyloOBJ, ordinat, color = 'TimeTreat') +
  scale_color_manual(
     values = c("#B2DF8A", "#33A02C", "#bdbdbd", "#bdbdbd",
                "#bdbdbd", "#bdbdbd",
                "#bdbdbd"),
     breaks = c("FR0", "FR6", "FN0", "FN6", 
                "PP0", "PP6",
                "DonB")) +
    geom_line(aes(group= Initial), color = "#bdbdbd", linetype = "dashed")

p3 <- p + stat_ellipse(data= subset(p$data, TimeTreat %in% c("FR0", "FR6")))+
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
p <- plot_ordination(phyloOBJ, ordinat, color = 'TimeTreat') +
 scale_color_manual(
     values = c("#bdbdbd", "#bdbdbd", "#B2DF8A", "#33A02C",
                "#bdbdbd", "#bdbdbd",
                "#bdbdbd"),
     breaks = c("FR0", "FR6", "FN0", "FN6", 
                "PP0", "PP6",
                "DonB")) +
  geom_line(aes(group= Initial), color = "#bdbdbd", linetype = "dashed")

p4 <- p + stat_ellipse(data= subset(p$data, TimeTreat %in% c("FN0", "FN6")))+
    annotate(geom = 'text', size=3, fontface= "bold",
           # adonis result
           label = paste("FMT NoRes:",
           adonis_ord(phyloOBJ, Remission %in% c("NoRes") &
                                Treatment == "FMT", METHOD = METHOD),
           sep = " "),
           x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    geom_point(alpha = 0.5) +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.position = "none")

# save_plot
plot_grid(p0, p3, p4, 
      ncol = 3, align = "hv")

}

```

#Aitchison’s distanc
Using raw count data to transform the data to centeriod log ratio and visualize
it in via euclidian distance in PCA space.

```{r}

#CLR transform
ps_metaphlan_clr <- microbiome::transform(ps_metaphlan, "clr")

psclr_dist= phyloseq::distance(ps_metaphlan_clr, method = "euclidean")
psclr_ord= ordinate(ps_metaphlan_clr, method = "RDA", distance = psclr_dist)

# save_plot
UCFMT1_ordinate(ps_metaphlan_clr, psclr_ord, "euclidean")
ggsave('figs/PCOA_Aithchinson_metaphlan_1.pdf', heigh = 6, width = 18, units = "cm")

# save_plot
UCFMT1_ordinate2(ps_metaphlan_clr, psclr_ord, "euclidean")
ggsave('figs/PCOA_Aithchinson_metaphlan_2.pdf', heigh = 6, width = 18, units = "cm")

```


# Bray-Curtis
Using relative abundance information to plot PCOA and UPCGMA Tree for All samples
and only DonorB datasets.

```{r}

bray_dis_ps = phyloseq::distance(ps_metaphlan_rel, method = 'bray')
bray_ord_ps = ordinate(ps_metaphlan_rel, method = 'PCoA', distance = bray_dis_ps)
# save_plot
UCFMT1_ordinate(ps_metaphlan_rel, bray_ord_ps, "bray")
ggsave('figs/PCOA_bray_metaphlan_1.pdf', heigh = 6, width = 18, units = "cm")

UCFMT1_ordinate2(ps_metaphlan_rel, bray_ord_ps, "bray")
ggsave('figs/PCOA_bray_metaphlan_2.pdf', heigh = 6, width = 18, units = "cm")


# UPGMA Tree for all PTS
bray_upgma_ps_meta = as.phylo(hclust(bray_dis_ps, method = 'average'))
write.tree(phy = bray_upgma_ps_meta, file = 'ITOL/bray_upgma_allPT_metaphlan.nwk')

```


#Engraftment Curve for donor B (based on relative abund) using metaphlan2
Extracting donorB species and look at DonorB treated patients and placebo

```{r}

#extracting donorB's ASVs:
subset_samples(ps_metaphlan_rel, donorB %in% c("DonB")) -> onlyDonB
psmelt(onlyDonB) -> df_DonB

df_DonB %>%
  select(OTU, Phylum, Family, Genus, Species, Sample, Abundance) %>%
  spread(Sample, Abundance, fill=0) -> df_DonB_wide

# selecting OTUs that are present at least in one DonorB samples:
#df_DonB_wide[rowSums(df_DonB_wide[,6:39] > 0) >= 1, ] %>%
df_DonB_wide[rowSums(df_DonB_wide[,6:9] > 0) >= 1, ] %>%
  pull(OTU)-> core_DonB

# Finding DonB OTUs in the whole DonorB dataset
subset_taxa(ps_metaphlan_rel, 
            rownames(tax_table(ps_metaphlan_rel)) %in% core_DonB) -> core_DonB_ps
core_DonB_ps <- prop_tax_down(core_DonB_ps, TRUE)


###############################################################################
# A function to get the commonly engrafted ASVs
###############################################################################
#removing DonB ASvs
subset_samples(core_DonB_ps, Fig_lab != "DonB") %>%
# Keeping both the placebo and FMT patients
  psmelt() %>%
  select(Fig_lab, Timepoint, Abundance, OTU, Phylum, Family, Genus, Species) %>%
   spread(key=Timepoint, value=Abundance) %>%
     #mutate(model= case_when(WK0 <= 0 & WK6 > 0 ~ paste("Engraft"),
     mutate(model= case_when(WK0 <= 0 & WK6 > 0.001 ~ paste("Engraft"),
                          WK0 > 0.001 & WK6 <= 0 ~ paste("Lost"),
                          WK0 > 0.001 & WK6 > 0.001~ paste("Present"),
                          WK0 <= 0 & WK6 <= 0 ~ paste("NotPresent"),
                            TRUE ~ paste("MinorChange"))) -> Sp_models
# Classified donor B ASVs in patients:
table(Sp_models$model)
# Engrafte ASV in WK6
 Sp_models %>%
  filter(model == "Engraft") %>%
  select(-WK0) %>%
  spread(Fig_lab, WK6, fill = 0) -> Sp_engraft
# Forloop to count the number of common engraftment across patients
EngraftOUT <- list()
EngraftCount <- list()
for (i in 1:length(unique(Sp_models$Fig_lab))) {
# selecting any ASV > 0 Across all the pts
EngraftOUT[[i]] <- Sp_engraft[rowSums(Sp_engraft[,7:ncol(Sp_engraft)] > 0) >= i, ]
EngraftCount[[i]] <- data_frame(Npts = paste(i),
                                eSpecies = length(EngraftOUT[[i]]$OTU))
}

###############################################################################

do.call(rbind.data.frame, EngraftCount) -> metaphlan_FMT_curve

#Visualize the output
metaphlan_FMT_curve %>%
  mutate(Treatment= paste("FMT")) %>%
  filter(as.numeric(Npts) <= 6) %>%
  ggplot(aes(as.numeric(Npts), eSpecies , colour = Treatment)) +
  geom_point(aes(group=Treatment))+ #, size=3.5) +
  geom_line(aes(group=Treatment))+#, size=3) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  scale_color_manual(
  values = c("#006837"),
  breaks = c("FMT"),
  labels = c("11 patients")) +
  theme_bw() +
  scale_y_sqrt() +
  theme(text = element_text(size = 10),
        #legend.position='none',
        legend.position=c(.6,.85),
        legend.title=element_blank(),
        strip.text = element_blank())+
   xlab("Number of patients") +
   ylab("Strains")-> fig_metaCurve
#ggsave('figs/Metaphlan_curve.pdf', heigh = 5, width = 7)


data.frame(EngraftOUT[3]) -> eSp_3pt
colnames(eSp_3pt)
eSp_3pt %>%
  gather(Fig_lab, Abundance, -OTU, -Phylum, -Family, -Genus, -Species,
         -model, na.rm = TRUE) %>%
  mutate(model= case_when(Abundance >= 0.001 ~ paste("Engraft"),
                          Abundance <= 0.001 ~ paste("Other"))) %>%
  left_join(
  data.frame(sample_data(ps_metaphlan_rel)) %>%
  filter(!Treatment == "Slurry") %>%
  select(Fig_lab, Remission) %>% distinct(), by = "Fig_lab") %>% 
  ggplot(aes(OTU, Fig_lab, fill=model)) +
  facet_grid(Remission+Fig_lab~Family, scales = "free", space = "free") +
  scale_fill_manual(values = c("#006d2c", "NA"),
                    breaks = c("Engraft", "Other")) +
  geom_tile() +
  theme_classic() +
  theme(axis.text= element_blank(),
        strip.text.x = element_text(angle=0), text = element_text(size = 10),
        legend.position = "none") + xlab("Species") + ylab("Donor B treated patients") 
  ggsave('figs/eMetaphlan_3pt.pdf', heigh = 10, width = 8, units = "cm")



```


Visualizing the composition of the genefamily file. It shows Gene families that 
are groups of evolutionarily-related protein-coding sequences that often 
perform similar functions. The normalized, merged genefamilies regrouped and 
then renamed all of these via humann2 code/pipes.


```{r}
#Importing the Uniref90(all humann2 marker) functions and their completename
gene_uniref90 <- read.csv("Humann2/UCFMT1_ShotgunRead_genefamilies_norm_uniref90.txt", 
                     sep = "\t", stringsAsFactors = FALSE)
# Edit the genefamily file
gene_uniref90 %>%
   rename_all(funs(str_replace(., "_clean_Abundance.RPKs", ""))) %>%
   rename_all(funs(str_replace(., "X", "s"))) %>%
   rownames_to_column("OTU") %>%
   mutate(OTU = paste("OTU", OTU, sep = "_")) %>%
   column_to_rownames("OTU") %>%
   dplyr::rename(GeneFamily = s..Gene.Family) %>%
   separate(GeneFamily, c("Genefamily", "Taxonomy"), sep = "[|]") %>%

# Removing taxonomic assignment from the data (Rows are replicated for taxa)
  filter(is.na(Taxonomy)) %>%
  mutate(Type= case_when(!Genefamily == "UniRef90_unknown" &
                         !Genefamily == "UNMAPPED" ~ paste("Uniref90_known"),
                         TRUE ~ as.character(Genefamily))) -> gene_uni_E_onlyfunc

#make the data long just for the sake of visualization
gene_uni_E_onlyfunc %>%
  gather(SampleName, Abundance, -Genefamily, -Taxonomy, -Type) %>%
  select(SampleName, Type, Abundance) %>%
  group_by(SampleName, Type) %>%
  summarize(NewAbundance= sum(Abundance)) %>%
  left_join(mapfile, by="SampleName") %>%
  mutate(Timepoint = case_when(Fig_lab == "DonB" ~ paste(X.SampleID),
                             TRUE ~ as.character(Timepoint))) %>%
  mutate(Timepoint = gsub("DMG_", "", Timepoint)) %>%
ggplot(aes(Timepoint, NewAbundance, fill = Type)) + 
  geom_bar(stat= "identity") +
  facet_grid(~Remission+Fig_lab, scales = "free", space = "free") +
  theme_classic() +
  scale_fill_manual(values = c(
    "#3182bd", "#636363", "#bdbdbd"), breaks = c(
    "Uniref90_known", "UniRef90_unknown", "UNMAPPED"
    )) +
  theme(axis.title.x = element_blank(),
        text = element_text(size=10),
        legend.position="bottom",
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1))
ggsave('figs/barplot_genefamily.pdf', height = 10, width = 14, units = "cm")

######################################################################
# Removing the unmapped and ungroupped data and make a phyloseq object
gene_uni_E_onlyfunc %>%
   filter(!Genefamily == "UNMAPPED" & !Genefamily == "UniRef90_unknown") %>%
   rownames_to_column("OTU") %>%
   mutate(OTU = paste("OTU", OTU, sep = "_")) %>%
   separate(Genefamily, c("Uniref90", "Uniref90_name"), sep = ":") %>%
   #removing those with uknown names
   filter(!Uniref90_name == " NO_NAME") %>%
   mutate(Genefamily = paste(Uniref90, Uniref90_name, sep = ":")) %>%
   select(-Uniref90, -Uniref90_name) %>%
   column_to_rownames("OTU") -> gene_uni_final

# gene info matrix
gene_uni_final %>%
  select(Genefamily, Type) -> gene_info
# Abundace matrix
gene_uni_final %>%
  select(-Genefamily, -Type, -Taxonomy) -> gene_matrix
# creat a phyloseq object contain OTU, Fuctions, and mapfile
gene_INF = tax_table(as.matrix(gene_info))
gene_OTU <- otu_table(as.matrix(gene_matrix), 
                      taxa_are_rows = TRUE)
ps_genefamily = phyloseq(gene_OTU, gene_INF)
sample_data(ps_genefamily) <- mapfile
                                                             
# the orgiinal genefamily file is relAbund but after removing unmappd/uknonw
ps_genefamily_rel <- ps_genefamily %>%
transform_sample_counts(function(x) {x/sum(x) * 100} )# Transform to rel. abundance
# here is sum of relative adbundace after removing unmapped and unknown uniref90
sample_sums(ps_genefamily)[1:10]
# and after re-cal relative abund
sample_sums(ps_genefamily_rel)[1:10]

```


#Aitchison’s distanc
Using relative abudance to transform the data to centeriod log ratio and visualize
it in via euclidian distance in PCA space.

```{r}
#CLR transform
(ps_clr <- microbiome::transform(ps_genefamily_rel, "clr"))  

psclr_dist= phyloseq::distance(ps_clr, method = "euclidean")
psclr_ord= ordinate(ps_clr , method = "RDA", distance = psclr_dist)

# save_plot
UCFMT1_ordinate(ps_clr, psclr_ord, "euclidean") 
ggsave('figs/RDA_Aithchinson_genefamily_1.pdf', heigh = 6, width = 18, units = "cm")

# save_plot
UCFMT1_ordinate2(ps_clr, psclr_ord, "euclidean") 
ggsave('figs/RDA_Aithchinson_genefamily_2.pdf', heigh = 6, width = 18, units = "cm")



```

# Bray-Curtis
Using relative abundance information to plot PCOA and UPCGMA Tree

```{r}

bray_dis_ps = phyloseq::distance(ps_genefamily_rel, method = 'bray')
bray_ord_ps = ordinate(ps_genefamily_rel, method = 'PCoA', 
                       distance = bray_dis_ps)
# save_plot
UCFMT1_ordinate(ps_genefamily_rel, bray_ord_ps, "bray") 
ggsave('figs/PCOA_bray_genefamily_1.pdf', heigh = 6, width = 18, units = "cm")

UCFMT1_ordinate2(ps_genefamily_rel, bray_ord_ps, "bray") 
ggsave('figs/PCOA_bray_genefamily_2.pdf', heigh = 6, width = 18, units = "cm")

# UPGMA Tree for all PTS
bray_upgma_ps_gene = as.phylo(hclust(bray_dis_ps, method = 'average'))
write.tree(phy = bray_upgma_ps_gene, file = 'ITOL/bray_upgma_allPT_genefamily.nwk')


```


#Engraftment Curve for donor B (relative abund) using genefamily (Human2)
Extracting donorB genefamily and look at DonorB treated patients and placebo

```{r}

#extracting donorB's ASVs:
subset_samples(ps_genefamily, donorB %in% c("DonB")) -> onlyDonB
source("psmelt.R") # a function that does the psmelt() but FAST
psmelt(onlyDonB, as = "df") -> df_DonB

df_DonB %>%
  select(OTU, Genefamily, Sample, Abundance) %>%
  spread(Sample, Abundance, fill=0) -> df_DonB_wide


# Genefamily detection threshold:
# metaphlan threshold is
# 0.0001(0.01%) = 0.1 coverage of 5Mbp
# x = 5 coverage 5Mbp=> x=0.005
# given 1000 genes per Mbp, we expect:
# 0.005(5x coverage) = 5 Mbp, x = 0.005 >> x = 1e-07
#cutoff <- 0.0000001
cutoff <- 0.000005

df_DonB_wide[rowSums(df_DonB_wide[,3:6] > cutoff) >= 4, ] %>%
  pull(OTU)-> core_DonB

# Finding DonB OTUs in the whole DonorB dataset
subset_taxa(ps_genefamily_rel, 
            rownames(tax_table(ps_genefamily_rel)) %in% core_DonB) -> core_DonB_ps

#removing DonB ASvs
subset_samples(core_DonB_ps, Fig_lab != "DonB")  -> temp
#subset_samples(!Treatment %in% c("Donor")) -> temp
  psmelt(temp, as = "df") %>%
  select(Fig_lab, Timepoint, Abundance, OTU, Genefamily) %>%
   spread(key=Timepoint, value=Abundance) %>%
     #mutate(model= case_when(WK0 <= 0 & WK6 > 0 ~ paste("Engraft"),
     mutate(model= case_when(WK0 <= 0 & WK6 > cutoff ~ paste("Engraft"),
                          WK0 > cutoff & WK6 <= 0 ~ paste("Lost"),
                          WK0 > cutoff & WK6 > cutoff ~ paste("Present"),
                          WK0 <= cutoff & WK6 <= cutoff~ paste("NotPresent"),
                            TRUE ~ paste("MinorChange"))) -> Gene_models
  format(summary(Gene_models$WK6), scientific = FALSE)

# Classified donor B ASVs in patients:
table(Gene_models$model)
# Engrafte ASV in WK6
Gene_models %>%
  filter(model == "Engraft") %>%
  select(-WK0) %>%
  spread(Fig_lab, WK6, fill = 0) -> Gene_engraft
# Forloop to count the number of common engraftment across patients
EngraftOUT <- list()
EngraftCount <- list()
for (i in 1:length(unique(Gene_models$Fig_lab))) {
# selecting any ASV > 0 Across all the pts
EngraftOUT[[i]] <- Gene_engraft[rowSums(Gene_engraft[,4:ncol(Gene_engraft)] > 0) >= i, ]
EngraftCount[[i]] <- data_frame(Npts = paste(i),
                                eSpecies = length(EngraftOUT[[i]]$OTU))
}
do.call(rbind.data.frame, EngraftCount)
###############################################################################


commonly_engrafted <- data.frame(EngraftOUT[3]) 
commonly_engrafted %>%
  gather(Fig_lab, Abundance, -OTU, -Genefamily, -model) %>%
  mutate(model= case_when(Abundance >= 0.001 ~ paste("Engraft"),
                          Abundance <= 0.001 ~ paste("Other"))) -> comm_fig
  
comm_fig$Fig_lab <- factor(comm_fig$Fig_lab,
  levels=rev(c("pt60","pt75","pt80","pt84",
           "pt25",
           "pt10", "pt4", "pt56", "pt74",  "pt79","pt85")))

  comm_fig %>%
   left_join(
  data.frame(sample_data(ps_metaphlan_rel)) %>%
  filter(!Treatment == "Slurry") %>%
  select(Fig_lab, Remission) %>% distinct(), by = "Fig_lab") %>% 
  ggplot(aes(OTU, Fig_lab, fill=model)) +
  scale_fill_manual(values = c("#006d2c", "NA"),
                    breaks = c("Engraft", "Other")) +
  geom_tile() +
  theme_classic() +
  facet_grid(Remission+Fig_lab~., scales = "free") +
  theme(legend.position = "none", text = element_text(size = 10),
        axis.text = element_blank()) + xlab("Genefamilies_uniref90") +
  ylab("Donor B treated patients") 
  ggsave('figs/commonlyEngraft_genefamily.pdf', heigh = 14, width = 8, units = "cm")

   
do.call(rbind.data.frame, EngraftCount) -> genefamily_FMT_curve

#Visualize the output
genefamily_FMT_curve %>%
  mutate(Treatment= paste("FMT")) %>%
  filter(as.numeric(Npts) <= 6) %>%
  ggplot(aes(as.numeric(Npts), eSpecies , colour = Treatment)) +
  geom_point(aes(group=Treatment))+#, size=3.5) +
  geom_line(aes(group=Treatment))+#, size=3) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  scale_color_manual(
  values = c("#006837"),
  breaks = c("FMT"),
  labels = c("11 patients")) +
  scale_y_sqrt() +
  theme_bw() +
  theme(text = element_text(size = 10),
         legend.position=c(.6,.85),
        legend.title=element_blank(),
        strip.text = element_blank()) +
  xlab("Number of patients") +
  ylab("Genefamilies") -> fig_GeneCurve

plot_grid(fig_metaCurve, fig_GeneCurve, ncol = 2, align = "hv")
ggsave('figs/Engraftment_curve.pdf', heigh = 6, width = 15, units = "cm")


```


ITOL mapfiles:
```{r}
###############################################################################
# Generate the sample labels file: (For both metaphlan/genefamilies files)
###############################################################################
samps_DonB = sample_data(ps_metaphlan_rel)
# Generate the tip label file
labs = data.frame(SampleID = samps_DonB$SampleName, 
				  #Label = paste(samps_DonB$Initial, 
				  #              samps_DonB$Timepoint, sep = '_')
				   Label = samps_DonB$Fig_lab)
ID_lab_head = c('LABELS',
				'SEPARATOR TAB',
				'DATASET_LABEL,IDc_labels',
				'DATA')
id_f = file('ITOL/ID_labs.txt')
writeLines(ID_lab_head, id_f)
close(id_f)
id_f = file('ITOL/ID_labs.txt', open = 'at')
write.table(labs, col.names = FALSE, file = id_f,
			sep = '\t', row.names = FALSE, qmethod = NULL, quote = FALSE)
close(id_f)

###############################################################################
# Generate the Timepoint_Treatment colouring file
###############################################################################
Remission_color <- read.csv("colours/Remission.txt", sep = "\t")
Remission = data.frame(SampleID = samps_DonB$SampleName,
				  TT = samps_DonB$Remission)
Remission %>% left_join(Remission_color, by = "TT") %>%
	select(-TT) -> Remission_cols
RR_head = c('DATASET_COLORSTRIP',	
			'SEPARATOR TAB',
			'DATASET_LABEL\tRemission',
			'COLOR_BRANCHES\t0',	
			'LEGEND_TITLE\tRemission',
			'LEGEND_SHAPES\t1\t1\t1\t1',
		    'LEGEND_COLORS\t#1B9E77\t#E7298A\t#A6761D\t#7570B3',
		    'LEGEND_LABELS\tResponder\tNo-responder\tPlacebo\tDonor',
		    'DATA')
RR_f = file('ITOL/DonBPTs_Remission.txt')
writeLines(RR_head, RR_f)
close(RR_f)
RR_f = file('ITOL/DonBPTs_Remission.txt', open = 'at')
write.table(Remission_cols, file = RR_f, sep='\t', row.names = FALSE, 
            col.names = FALSE,
			qmethod = NULL, quote = FALSE)
close(RR_f)
################################################################################
## Generate the patient colour file
################################################################################
patient_color <- read.csv("colours/metagenome_PTs.txt", sep = "\t")
pats_DonB = samps_DonB[,c('X.SampleID','Initial')]
pats_DonB = data.frame(SampleID= samps_DonB$SampleName,
                  Initial = samps_DonB$Initial)
pats_DonB %>% left_join(patient_color, by = "Initial") %>%
	select(-Initial) %>%
	mutate(clade = 'clade',normal = 'normal') %>%
	select(SampleID, clade, colour, normal)-> pat_DonB
pat_head = c('TREE_COLORS',
			 'SEPARATOR TAB',
			 'DATA')
pat_f_DonB = file('ITOL/DonBPTs_cladeColor.txt')
writeLines(pat_head, pat_f_DonB)
close(pat_f_DonB)
pat_f_DonB = file('ITOL/DonBPTs_cladeColor.txt', open = 'at')
write.table(pat_DonB, file = pat_f_DonB, sep = '\t', row.names = FALSE, 
            col.names = FALSE,
			qmethod = NULL, quote = FALSE)
close(pat_f_DonB)
################################################################################
## Generate the taxa bar chart file for metaphlan
################################################################################
#import the phylum color, compatible with 16S data
phy_col <- read.csv("colours/phylum_colours.txt", sep = "\t")

ps_phylum_df %>%
	select(SampleName, Phylum, Abundance) %>%
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
ps_phy_f = file('ITOL/DonBPTs_taxa_bar_phylum.txt')
writeLines(ps_phy_head, ps_phy_f)
close(ps_phy_f)
ps_phy_f = file('ITOL/DonBPTs_taxa_bar_phylum.txt', open = 'at')
write.table(ps_phyl_itol, file = ps_phy_f, sep = '\t', 
			row.names = FALSE, col.names = FALSE, qmethod = NULL, quote = FALSE)
close(ps_phy_f)

################################################################################
## Generate the taxa bar chart file for genefamily
################################################################################
gene_uni_E_onlyfunc %>%
  gather(SampleName, Abundance, -Genefamily, -Taxonomy, -Type) %>%
  select(SampleName, Type, Abundance) %>%
  group_by(SampleName, Type) %>%
  summarize(NewAbundance= sum(Abundance)) %>%
  left_join(mapfile, by="SampleName") %>%
  filter(!SampleName %in% Not_nec) %>%
  mutate(Timepoint = case_when(Fig_lab == "DonB" ~ paste(X.SampleID),
                             TRUE ~ as.character(Timepoint))) %>%
  select(SampleName, Type, NewAbundance) %>%
  spread(Type, NewAbundance, fill = 0) -> genefamily_itol

genefamily_name = c("Uniref90_known", "UniRef90_unknown", "UNMAPPED")
genefamily_here= c("#3182bd","#636363","#bdbdbd")
ps_gene_head = c('DATASET_MULTIBAR',
			 'SEPARATOR TAB',
			 'DATASET_LABEL\tGenefamily',
			 stri_join(c('FIELD_LABELS', genefamily_name), collapse = '\t'),
			 stri_join(c('FIELD_COLORS',genefamily_here), collapse = '\t'),
			 'LEGEND_TITLE\tGenefamily',
			 stri_join(c('LEGEND_SHAPES',rep('1',length(genefamily_name))), 
			 		  collapse = '\t'),
			 stri_join(c('LEGEND_COLORS',genefamily_here), 
			 		  collapse = '\t'),
			 stri_join(c('LEGEND_LABELS', genefamily_name),
			 		  collapse = '\t'),
			 'DATA')
ps_gene_f = file('ITOL/DonBPTs_taxa_bar_genefamily.txt')
writeLines(ps_gene_head, ps_gene_f)
close(ps_gene_f)
ps_gene_f = file('ITOL/DonBPTs_taxa_bar_genefamily.txt', open = 'at')
write.table(genefamily_itol, file = ps_gene_f, sep = '\t', 
			row.names = FALSE, col.names = FALSE, qmethod = NULL, quote = FALSE)
close(ps_gene_f)


```

