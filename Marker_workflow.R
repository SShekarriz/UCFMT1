---
title: "Marker_analysis"
author: "Sharok"
date: "7/1/2021"
output: html_document
---

We built species-specific marker using pangenome alignments. In order to test
the accuracy of CEGs as strain marker and our specices-specific marker a 
collection of 1.2K genomes from surette lab were mapped to these markers with
stringent cutoff: bwa mem -O 60 -E 10 -L 100

```{r}
library(tidyverse)
library(tidytext)
```

# Surette Genome collection taxonomy
```{r}
#importing GC taxonomy
GC_taxa <- read.csv("GC_gtdbtk_taxonomy_2021.csv")
GC_taxa%>% 
  separate(classification, c("Kingdom", "Phylum", "Class",
                                   "Order", "Family", "Genus",
                                   "Species"), ";.__") %>%
  select(user_genome, Kingdom:Species) -> GC_taxa
```

# genes from rsCEGs in the phylogeny
```{r, echo=FALSE, message=FALSE, warning=FALSE}
path="genes_from_rsCEGs"
patt="_CEGs.blastout"
COLS <- c("Genome","qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
"qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp")
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = F))) %>%
         unnest() %>%
  mutate(sample.id = gsub(patt, "", sample.id)) -> blastOuts
colnames(blastOuts) <- COLS

blastOuts %>%
    filter(pident >= 90, qcovhsp >= 90) %>%
    mutate(temp=sseqid) %>%
    separate(temp, c("Ref_genomes", "Ref_contig"), sep = ";") %>%
    filter(Ref_genomes %in% c("GCF_902497355.1_P9094",
                              "GC313_hybrid",
                              "GC568")) %>%
   rename(gene_id= qseqid) %>%
   select(gene_id, Genome, Ref_genomes, Ref_contig) %>%
   mutate(gene_id= paste("g", gene_id, sep = "_"))-> gene_in_rsCEGs

# Number of CEGs within each representative
table(gene_in_rsCEGs$Genome)
```


# Species markers- genomes
```{r}
path="perfect_align_freeMismatch/species_markers"
patt="_coverage.txt"

data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(Contig= `#rname` ) %>%
         # make sure to have "__" for all contigs
         separate(Contig, c("Contig", "Cluster"), sep = "__") %>%
         # make a new variable for genome names  
         mutate(Genome = case_when(grepl("Isolate_19", Contig) ~ paste("Fus"),
                                   grepl("RefSeq_89", Contig) ~ paste("Fpr"),
                                   grepl("GC568", Contig) ~ paste("Dor"),
                                   TRUE ~ as.character("Other"))) %>%
         mutate(length= endpos) %>%
         mutate(Type=paste("Conserved"))-> cons
```


# Strain marker -commonly engrafted genes from phylogeny-genomes
```{r, echo=FALSE, message=FALSE, warning=FALSE}

path="perfect_align_freeMismatch/CEGenes"
patt="_coverage.txt"
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(gene_id= `#rname`) %>%
         mutate(gene_id= paste("g", gene_id, sep = "_")) %>%
         left_join(gene_in_rsCEGs, by = "gene_id") %>%
         filter(!is.na(Genome)) %>% 
         # selecting gene with minimum 100 bps
         filter(endpos >= 100) %>%
  # all the genes within a genome together
  group_by(sample.id, Genome) %>%
  summarise(length= sum(endpos),
            numreads=sum(numreads),
            covbases=sum(covbases),
            coverage=covbases/length * 100) %>%
  mutate(Type= paste("CEGs"))-> cegs

```


# visualize genomic data for marker (test the accuracy)
```{r}

#separating the color of Fus strains
Fus_of_int <- c("GC313", "GC313_hybrid", "GC828", "GC917")

# merge the two tables
cons %>%
  select(colnames(cegs)) %>%
  bind_rows(cegs) %>%
  left_join(GC_taxa %>% rename(sample.id= user_genome),
                               by = "sample.id") %>%
  mutate(taxa= case_when(
    sample.id %in% Fus_of_int ~ paste("Fusicatenibacter saccharivorans_A"),
    grepl("Dorea", Genus) ~ paste(Species),
    grepl("Fusicatenibacter", Genus) ~ paste(Species),
    grepl("Faecalibacterium", Genus) ~ paste(Species),
    coverage > 80 ~ paste(Species),
    TRUE ~ paste("Other"))
    ) %>%
  mutate(Cluster= paste(Type, paste(length, "bp", sep = ""),
                        sep = ":"))-> tbl

sp_col <- c("Other" = "black",
            "NA" = "grey",
            "Dorea_A longicatena" = "#d95f02",
            "Dorea_A longicatena_B" =  "#1b9e77",
            "Dorea formicigenerans" = "#7570b3",
            "Fusicatenibacter saccharivorans" = "#e7298a",
            "Fusicatenibacter saccharivorans_A" = "#66a61e",
            "Faecalibacterium prausnitzii" =  "#e6ab02")

ggplot(tbl, aes(sample.id, coverage, color= taxa)) +
  geom_point(size=2.5, alpha=0.6) +
  geom_hline(yintercept=80, linetype="dashed", 
                color = "black") +
  facet_wrap(~Genome~Cluster, scales = "free_x", ncol = 2) +
  scale_color_manual(values = sp_col)+
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=10),
        legend.text = element_text(size=6),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank())
ggsave("marker_accuracy1.png", width = 12, height = 18, units = "cm")

```



```{r}

tbl %>%
select(sample.id, Genome, Type, coverage) %>%
spread(Type, coverage) %>%
left_join(GC_taxa %>% rename(sample.id= user_genome),
                               by = "sample.id") %>%
  mutate(taxa= case_when(
    sample.id %in% Fus_of_int ~ paste("Fusicatenibacter saccharivorans_A"),
    grepl("Dorea", Genus) ~ paste(Species),
    grepl("Fusicatenibacter", Genus) ~ paste(Species),
    grepl("Faecalibacterium", Genus) ~ paste(Species),
    Conserved > 80 ~ paste(Species),
    CEGs > 80 ~ paste(Species),
    TRUE ~ paste("Other"))
    ) -> tbl2

ggplot(tbl2, aes(CEGs, Conserved, color=taxa, fill=taxa)) +
    # setting the background colors
   annotate("rect", xmin=80,xmax=Inf,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#008837") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#a6dba0") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#bababa") +
   annotate("rect", xmin=80,xmax=Inf,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#404040") +
  # main plot component
  geom_point(alpha=0.6,size=2.5,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=80, linetype="dashed", 
                color = "black") +
  geom_vline(xintercept=80, linetype="dashed", 
                color = "black") +
  facet_wrap(~Genome, ncol = 1) +
  scale_color_manual(values = sp_col) +
  scale_fill_manual(values = sp_col) +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank())
ggsave("marker_accuracy2.png", width = 6, height = 18, units = "cm")


```



# UCFMT mapfile
```{r, echo=FALSE, message=FALSE, warning=FALSE}
pt_Fate <- data.frame(Fig_lab= c(
"pt60", "pt75", "pt80", "pt84",
"pt25",
"pt10", "pt4", "pt56", "pt74","pt79", "pt85",
"DonB_O17", "DonB_M17", "DonB_13", "DonB_16"),
Fate=c(
 "N", "N", "N", "N",
 "P",
 "R", "R", "R", "R", "R", "R",
 "H", "H", "H", "H"
))
mapfile <- read.csv("map_healthy_vs_UC.txt", sep="\t", header = F)
head <- c("sample.id", "Fig_lab", "Timepoint", "Status")
colnames(mapfile) <- head
#depth of seq
d_UC <- read.csv("UC_seq_depth.txt", header = FALSE)
d_B <- read.csv("DonorB_seq_depth.txt", header = FALSE)
d_UC %>%
  bind_rows(d_B) %>%
  separate(V1, c("sample.id", "ReadN"), sep = "__") %>%
  mutate(sample.id = gsub("_001_inter.fastq", "", sample.id))-> depth
mapfile %>% left_join(depth) %>% left_join(pt_Fate)-> mapfile

```



# UCFMT dataset
```{r, echo=FALSE, message=FALSE, warning=FALSE}

path="perfect_align_freeMismatch/species_markers_for_FMTdata"
patt="_coverage.txt"
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(Contig= `#rname` ) %>%
         # make sure to have "__" for all contigs
         separate(Contig, c("Contig", "Cluster"), sep = "__") %>%
         # make a new variable for genome names  
         mutate(Genome = case_when(grepl("Isolate_19", Contig) ~ paste("Fus"),
                                   grepl("RefSeq_89", Contig) ~ paste("Fpr"),
                                   grepl("GC568", Contig) ~ paste("Dor"),
                                   TRUE ~ as.character("Other"))) %>%
         mutate(length= endpos) %>%
         mutate(Type=paste("Conserved"))-> cons


path="perfect_align_freeMismatch/CEGenes_for_FMTdata"
patt="_coverage.txt"
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(gene_id= `#rname`) %>%
         mutate(gene_id= paste("g", gene_id, sep = "_")) %>%
         left_join(gene_in_rsCEGs, by = "gene_id") %>%
         filter(!is.na(Genome)) %>% 
         # selecting gene with minimum 100 bps
         filter(endpos >= 100) %>%
  # all the genes within a genome together
  group_by(sample.id, Genome) %>%
  summarise(length= sum(endpos),
            numreads=sum(numreads),
            covbases=sum(covbases),
            coverage=covbases/length * 100) %>%
  mutate(Type= paste("CEGs"))-> cegs


# merge the two tables
cons %>%
  select(colnames(cegs)) %>%
  bind_rows(cegs) %>%
  filter(!sample.id %in% c("721_CCGTCC", "797_GTGGCC")) %>%
  mutate(Cluster= paste(Type, paste(length, "bp", sep = ""),
                        sep = ":")) %>%
  left_join(mapfile) %>%
  filter(!is.na(Fig_lab)) %>%
  filter(Fate != "P") %>% #remove Placebo
  mutate(Time= case_when(Timepoint == "H" ~ paste("DonorB(n=4)"),
                         Timepoint == "0" ~ paste("UC_base(n=10)"),
                         Timepoint == "6" ~ paste("UC_FMT(n=10)"))) -> tbl_fmt
  

```



# metagenomics dataset: PRJNA400072, number of individuals:
```{r, echo=FALSE, message=FALSE, warning=FALSE}
#original study
#https://www.nature.com/articles/s41564-018-0306-4
map_sra <- read.csv("PRJNA400072_SraRunTable.txt")
Diagnos <- read.csv("PRJNA400072_Diagnosis.csv")
Diagnos %>%
  mutate(Diagnosis= gsub("CD", "CD(n=88)", Diagnosis)) %>%
  mutate(Diagnosis= gsub("Control", "Control(n=56)", Diagnosis)) %>%
  mutate(Diagnosis= gsub("UC", "UC(n=76)", Diagnosis)) %>%
  mutate(Sample.Name= case_when(grepl("^[0-9]", 
                                      local_sample_id) ~ paste("PRISM",
                                                               local_sample_id,
                                                               sep = "_"),
                                TRUE ~ as.character(local_sample_id))) %>%
  select(Sample.Name, Diagnosis) %>%
  left_join(map_sra, by = "Sample.Name")-> map_sra
table(map_sra$Diagnosis)

```


# PRJNA400072 dataset
```{r}

path="perfect_align_freeMismatch/species_markers_for_PRJNA400072"
patt="_coverage.txt"
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(Contig= `#rname` ) %>%
         # make sure to have "__" for all contigs
         separate(Contig, c("Contig", "Cluster"), sep = "__") %>%
         # make a new variable for genome names  
         mutate(Genome = case_when(grepl("Isolate_19", Contig) ~ paste("Fus"),
                                   grepl("RefSeq_89", Contig) ~ paste("Fpr"),
                                   grepl("GC568", Contig) ~ paste("Dor"),
                                   TRUE ~ as.character("Other"))) %>%
         mutate(length= endpos) %>%
         mutate(Type=paste("Conserved"))-> cons


path="perfect_align_freeMismatch/CEGenes_for_PRJNA400072"
patt="_coverage.txt"
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(gene_id= `#rname`) %>%
         mutate(gene_id= paste("g", gene_id, sep = "_")) %>%
         left_join(gene_in_rsCEGs, by = "gene_id") %>%
         filter(!is.na(Genome)) %>% 
         # selecting gene with minimum 100 bps
         filter(endpos >= 100) %>%
  # all the genes within a genome together
  group_by(sample.id, Genome) %>%
  summarise(length= sum(endpos),
            numreads=sum(numreads),
            covbases=sum(covbases),
            coverage=covbases/length * 100) %>%
  mutate(Type= paste("CEGs"))-> cegs


# merge the two tables
cons %>%
  select(colnames(cegs)) %>%
  bind_rows(cegs) %>%
  left_join(map_sra %>% select(Run, Diagnosis, AvgSpotLen, Bases) %>%
              rename(sample.id=Run), 
            by="sample.id") %>%
  mutate(Cluster= paste(Type, paste(length, "bp", sep = ""),
                        sep = ":")) %>%
  mutate(Depth = Bases / AvgSpotLen) -> tbl_PRJNA400072


```




```{r}

strain_col <- c("Fus" = "#e41a1c",
                "Dor" = "#377eb8",
                "Fpr" = "#984ea3")

#Surette dataset:
tbl_fmt %>%
  rename(Diagnosis = Time) %>%
  mutate(Dataset=paste("Surette")) %>%
  select(sample.id:Cluster, Diagnosis, Dataset) %>%
#Add SRA dataset
  bind_rows(tbl_PRJNA400072 %>%
            select(sample.id:Diagnosis, Cluster) %>%
            mutate(Dataset=paste("SRA"))) -> tbl

tbl %>%
  select(sample.id, Genome, Type, coverage) %>%
  spread(Type, coverage) %>%
  left_join(
  tbl %>% select(sample.id, Diagnosis, Dataset) %>%
  distinct(), by = "sample.id") -> tbl_viz

tbl_viz$Diagnosis <- factor(tbl_viz$Diagnosis,
                            levels = c("DonorB(n=4)", "UC_base(n=10)", 
                                       "UC_FMT(n=10)",
                                       "Control(n=56)", "UC(n=76)", 
                                       "CD(n=88)"))

ggplot(tbl_viz, aes(CEGs, Conserved, color=Genome, fill=Genome)) +
    # setting the background colors
   annotate("rect", xmin=80,xmax=Inf,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#008837") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#a6dba0") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#bababa") +
   annotate("rect", xmin=80,xmax=Inf,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#404040") +
  # main plot component
  geom_point(alpha=0.6,position="dodge", size=2.5) +
  geom_hline(yintercept=80, linetype="dashed", 
                color = "black") +
  geom_vline(xintercept=80, linetype="dashed", 
                color = "black") +
  facet_wrap(~Diagnosis) +
  scale_color_manual(values = strain_col) +
  scale_fill_manual(values = strain_col) +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank())
ggsave("marker_in_metagenome.png", width = 15, height = 10, units = "cm")

```





```{r}

ggplot(tbl_viz %>% filter(Genome=="Dor"), 
       aes(CEGs, Conserved, color=Genome, fill=Genome)) +
    # setting the background colors
   annotate("rect", xmin=80,xmax=Inf,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#008837") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#a6dba0") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#bababa") +
   annotate("rect", xmin=80,xmax=Inf,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#404040") +
  # main plot component
  geom_point(alpha=0.6,position="dodge", size=2.5) +
  geom_hline(yintercept=80, linetype="dashed", 
                color = "black") +
  geom_vline(xintercept=80, linetype="dashed", 
                color = "black") +
  facet_wrap(~Diagnosis) +
  scale_color_manual(values = strain_col) +
  scale_fill_manual(values = strain_col) +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_blank())
ggsave("Dor_in_metagenome.png", width = 15, height = 10, units = "cm")


ggplot(tbl_viz %>% filter(Genome=="Fpr"), 
       aes(CEGs, Conserved, color=Genome, fill=Genome)) +
    # setting the background colors
   annotate("rect", xmin=80,xmax=Inf,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#008837") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#a6dba0") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#bababa") +
   annotate("rect", xmin=80,xmax=Inf,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#404040") +
  # main plot component
  geom_point(alpha=0.6,position="dodge", size=2.5) +
  geom_hline(yintercept=80, linetype="dashed", 
                color = "black") +
  geom_vline(xintercept=80, linetype="dashed", 
                color = "black") +
  facet_wrap(~Diagnosis) +
  scale_color_manual(values = strain_col) +
  scale_fill_manual(values = strain_col) +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_blank())
ggsave("Fpr_in_metagenome.png", width = 15, height = 10, units = "cm")


ggplot(tbl_viz %>% filter(Genome=="Fus"), 
       aes(CEGs, Conserved, color=Genome, fill=Genome)) +
    # setting the background colors
   annotate("rect", xmin=80,xmax=Inf,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#008837") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#a6dba0") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#bababa") +
   annotate("rect", xmin=80,xmax=Inf,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#404040") +
  # main plot component
  geom_point(alpha=0.6,position="dodge", size=2.5) +
  geom_hline(yintercept=80, linetype="dashed", 
                color = "black") +
  geom_vline(xintercept=80, linetype="dashed", 
                color = "black") +
  facet_wrap(~Diagnosis) +
  scale_color_manual(values = strain_col) +
  scale_fill_manual(values = strain_col) +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_blank())
ggsave("Fus_in_metagenome.png", width = 15, height = 10, units = "cm")

```




# SHCM-DMG dataset
```{r, echo=FALSE, message=FALSE, warning=FALSE}

path="perfect_align_freeMismatch/species_markers_for_SHCM-DMG"
patt="_coverage.txt"
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(Contig= `#rname` ) %>%
         # make sure to have "__" for all contigs
         separate(Contig, c("Contig", "Cluster"), sep = "__") %>%
         # make a new variable for genome names  
         mutate(Genome = case_when(grepl("Isolate_19", Contig) ~ paste("Fus"),
                                   grepl("RefSeq_89", Contig) ~ paste("Fpr"),
                                   grepl("GC568", Contig) ~ paste("Dor"),
                                   TRUE ~ as.character("Other"))) %>%
         mutate(length= endpos) %>%
         mutate(Type=paste("Conserved"))-> cons


path="perfect_align_freeMismatch/CEGenes_for_SHCM-DMG"
patt="_coverage.txt"
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = T))) %>%
         unnest() %>%
         mutate(sample.id = gsub(patt, "", sample.id)) %>%
         dplyr::rename(gene_id= `#rname`) %>%
         mutate(gene_id= paste("g", gene_id, sep = "_")) %>%
         left_join(gene_in_rsCEGs, by = "gene_id") %>%
         filter(!is.na(Genome)) %>% 
         # selecting gene with minimum 100 bps
         filter(endpos >= 100) %>%
  # all the genes within a genome together
  group_by(sample.id, Genome) %>%
  summarise(length= sum(endpos),
            numreads=sum(numreads),
            covbases=sum(covbases),
            coverage=covbases/length * 100) %>%
  mutate(Type= paste("CEGs"))-> cegs


# merge the two tables
cons %>%
  select(colnames(cegs)) %>%
  bind_rows(cegs) %>%
  mutate(Cluster= paste(Type, paste(length, "bp", sep = ""),
                        sep = ":")) %>%
  mutate(Donor= gsub("_Stool.*", "", sample.id)) %>%
  mutate(Donor= gsub("-stool.*", "", Donor)) %>%
  mutate(Donor= gsub("Stool.*", "", Donor)) -> tbl_SHCM




strain_col <- c("Fus" = "#e41a1c",
                "Dor" = "#377eb8",
                "Fpr" = "#984ea3")


tbl_SHCM %>%
  select(sample.id, Genome, Type, coverage, Donor) %>%
  spread(Type, coverage) -> tbl_viz

ggplot(tbl_viz, aes(CEGs, Conserved, color=Genome, fill=Genome)) +
    # setting the background colors
   annotate("rect", xmin=80,xmax=Inf,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#008837") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=80,ymax=Inf, 
            alpha=0.5, fill="#a6dba0") +
   annotate("rect", xmin=-Inf,xmax=80,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#bababa") +
   annotate("rect", xmin=80,xmax=Inf,ymin=-Inf,ymax=80, 
            alpha=0.5, fill="#404040") +
  # main plot component
  geom_point(alpha=0.6,position="dodge", size=2.5) +
  geom_hline(yintercept=80, linetype="dashed", 
                color = "black") +
  geom_vline(xintercept=80, linetype="dashed", 
                color = "black") +
  facet_wrap(~Donor) +
  scale_y_continuous(limits=c(0,100)) +
  scale_x_continuous(limits=c(0,100)) +
  scale_color_manual(values = strain_col) +
  scale_fill_manual(values = strain_col) +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_blank())
ggsave("marker_in_SHCM-DMG.png", width = 15, height = 10, units = "cm")

```




```{r}

q_col <- c("B1" = "#008837",
           "B2" = "#a6dba0",
           "B3" = "#bababa",
           "B4" = "#404040",
           "Dor" = "#b9cedf",
           "Fpr" = "#d5bfd9",
           "Fus" = "#f4a4a6")

tbl_viz %>%
mutate(Quarter = case_when(Conserved > 80 & CEGs > 80 ~ paste("B1"),
                           Conserved > 80 & CEGs < 80 ~ paste("B2"),
                           Conserved < 80 & CEGs < 80 ~ paste("B3"),
                           Conserved < 80 & CEGs > 80 ~ paste("B4"),
                             TRUE ~ paste("Other"))) -> tbl_spec


tbl_spec %>%
  group_by(Diagnosis, Genome) %>% summarise(total=n()) -> total_G

tbl_spec %>%
 group_by(Diagnosis, Quarter, Genome) %>% tally() %>%
 left_join(total_G, by = c("Diagnosis", "Genome")) %>%
 mutate(percent= n/total * 100) %>%
 ggplot(aes(Quarter, percent, fill=Quarter)) +
 geom_rect(aes(fill=Genome), xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf) +
  geom_bar(stat = "identity") +
  facet_wrap(~Diagnosis+Genome, ncol = 3) +
  scale_fill_manual(values = q_col) +
  theme_classic() +
  geom_text(aes(label=round(percent, digits = 0)))+
  theme(text = element_text(size=10),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) 
  ggsave("marker_in_metagenome2.png", width = 6, height = 12, units = "cm")

```




```{r}
library(gtools)

tbl_spec %>%
 group_by(Diagnosis, Quarter, Genome) %>% tally() %>%
 left_join(total_G, by = c("Diagnosis", "Genome")) %>%
 mutate(percent= n/total * 100) -> propor_tbl

propor_tbl %>%
  filter(Diagnosis %in% c("UC_base(n=10)", "UC_FMT(n=10)")) %>%
  select(-n, -total) %>%
  ungroup() %>%
  spread(Diagnosis, percent, fill = 0) %>%
  mutate(Type= paste("FMT"))-> propor_FMT
propor_FMT$foldchange <- foldchange(as.numeric(propor_FMT$`UC_FMT(n=10)`), 
                                    as.numeric(propor_FMT$`UC_base(n=10)`))

propor_tbl %>%
  filter(Diagnosis %in% c("UC(n=76)", "Control(n=56)")) %>%
  select(-n, -total) %>%
  ungroup() %>%
  spread(Diagnosis, percent, fill = 0) %>%
  mutate(Type= paste("UC-Healthy"))-> propor_UC
propor_UC$foldchange <- foldchange(as.numeric(propor_UC$`Control(n=56)`), 
                                    as.numeric(propor_UC$`UC(n=76)`))

propor_tbl %>%
  filter(Diagnosis %in% c("CD(n=88)", "Control(n=56)")) %>%
  select(-n, -total) %>%
  ungroup() %>%
  spread(Diagnosis, percent, fill = 0) %>%
  mutate(Type= paste("CD-Healthy")) -> propor_CD
propor_CD$foldchange <- foldchange(as.numeric(propor_CD$`Control(n=56)`), 
                                    as.numeric(propor_CD$`CD(n=88)`))

propor_FMT %>%
  rename("0" = "UC_base(n=10)",
         "1" = "UC_FMT(n=10)") %>%
  bind_rows(propor_UC %>%
            rename("0" = "UC(n=76)",
                   "1" = "Control(n=56)")) %>%
  bind_rows(propor_CD %>%
            rename("0" = "CD(n=88)",
                   "1" = "Control(n=56)")) -> propor
propor$foldchange2logratio <- foldchange2logratio(propor$foldchange)
    
  
ggplot(propor, aes(Quarter, foldchange2logratio, color=Type)) +
  geom_point() +
  facet_grid(~Genome) +
  theme_bw()

```

