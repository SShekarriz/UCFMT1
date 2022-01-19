---
title: "Phylogenomics"
author: "Sharok"
date: "2/12/2021"
output: html_document
---

#Dorea Genomes and CEGs
```{r}
library(tidyverse)

res <- read.csv("gtdbtk.bac120.summary.tsv", sep = "\t")
res %>% separate(classification, c("Kingdom", "Phylum", "Class",
                                   "Order", "Family", "Genus",
                                   "Species"), ";.__") %>%
  select(user_genome, Kingdom:Species) %>%
  mutate(Species = gsub("_B", "", Species)) %>%
  mutate(Species = gsub("sp.*", "sp.", Species)) %>%
  select(user_genome, Species)-> sp_labels

###############################################################################
# selecting only Dorea longicatena genome for further core-genome analysis
###############################################################################
write.table(
sp_labels %>%
  filter(grepl("longicatena", Species)) %>%
  select(user_genome), "Only_Dlongicatena.txt",
sep = "\t", col.names = F, row.names = F, quote = F
)
###############################################################################
# Generate the genome labels file for iTOL
###############################################################################
ID_lab_head = c('LABELS',
				'SEPARATOR TAB',
				'DATASET_LABEL,IDc_labels',
				'DATA')
id_f = file('AllGenome_ID_labs.txt')
writeLines(ID_lab_head, id_f)
close(id_f)
id_f = file('AllGenome_ID_labs.txt', open = 'at')
write.table(sp_labels, col.names = FALSE, file = id_f,
			sep = '\t', row.names = FALSE, qmethod = NULL, quote = FALSE)
close(id_f)

###############################################################################
# Generate a tree bars for the engrafted genes
###############################################################################

blastOuts <- read.table("allfna_Common-engraft.blastout", sep = "\t",
                      header = FALSE)
COLS <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
"qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp")
colnames(blastOuts) <- COLS
  
blastOuts %>%
    filter(pident >= 90, qcovhsp >= 90) %>%
    mutate(temp=sseqid) %>%
    separate(temp, c("Ref_genomes", "Ref_contig"), sep = ";") %>%
    group_by(qseqid, Ref_genomes) %>%
    summarise(Gene_counts_within_genome=n()) %>%
    distinct()  %>%
    group_by(Ref_genomes) %>%
    summarise(Gene_counts_within_genome=n()) %>%
    filter(Gene_counts_within_genome >= 1) -> Genes_per_genomes

ID_lab_head = c('DATASET_SIMPLEBAR',
				'SEPARATOR COMMA',
				'DATASET_LABEL,CEGs',
				'COLOR,#006d2c',
				'DATA')
id_f = file('Genome_CEGs_bars.txt')
writeLines(ID_lab_head, id_f)
close(id_f)
id_f = file('Genome_CEGs_bars.txt', open = 'at')
write.table(Genes_per_genomes, col.names = FALSE, file = id_f,
			sep = ',', row.names = FALSE, qmethod = NULL, quote = FALSE)
close(id_f)

```


#Phylogenetic tree1 (marker based)
```{r}

library(ggtree)
library(ape)
library(ggtreeExtra)
library(ggnewscale)
library(grid)
library(gtable)

res %>% separate(classification, c("Kingdom", "Phylum", "Class",
                                   "Order", "Family", "Genus",
                                   "Species"), ";.__") %>%
  select(user_genome, Kingdom:Species) %>%
  mutate(Lab = gsub("_B", "", Species)) %>%
  mutate(Lab = gsub("sp.*", "sp.", Lab)) %>%
  mutate(Source = case_when(grepl("GCF_", user_genome) ~ paste("RefSeq"),
                            TRUE ~ as.character("Isolate"))) %>%
  group_by(Source) %>% mutate(gCount = 1:n()) %>%
  ungroup() %>%
  mutate(Lab.org = paste(Source, gCount, sep="_")) %>%
  left_join(Genes_per_genomes %>%
              rename(user_genome="Ref_genomes",
                     CEGs="Gene_counts_within_genome"),
                     by="user_genome") %>%
  mutate_each(funs(replace(., which(is.na(.)), 0))) %>%
  select(user_genome, Species, Lab, Lab.org, CEGs)-> info

write.table(info, "Dorea_complete_info.txt", sep = "\t", 
            row.names = F, quote = F)

tree <- read.tree("gtdbtk.bac120.user_msa.tree")

cols <- c("Dorea formicigenerans"='#377eb8', 
          "Dorea longicatena"='#006d2c', 
          "Dorea sp."='#984ea3')

info %>% select(user_genome,Lab,Lab.org) -> info2
bar_tbl <- info %>% select(user_genome, CEGs)

###############################################################################
# No annotation
###############################################################################

p <- ggtree(tree) 
## attach the sampling information data set 
## and add symbols colored by location
p <- p %<+% info2 + geom_tippoint(aes(color=Lab), size=1) +
                    geom_tiplab(aes(label=Lab.org, color=Lab), 
                                size=1.5, align = T, fontface="bold") +
                    scale_color_manual(values=cols) +
                    scale_fill_manual(values =cols)
  
## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
p2 <- p + geom_facet(panel = "CEGs", data = bar_tbl, 
                     geom = ggstance::geom_barh, 
             aes(x = CEGs, color = Lab, fill = Lab), 
             stat = "identity", width = .6) + xlim_tree(0.25) +
   theme_tree2(legend.position=c(0.2, 0.85)) +
   theme(strip.text = element_blank(), legend.title = element_blank(),
     strip.background.x=element_rect(fill=NA))


gt = ggplot_gtable(ggplot_build(p2))
#gtable_show_layout(gt) # will show you the layout - very handy function
gt # see plot layout in table format
gt$layout$l[grep('panel-1-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[7] = 0.5*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
#grid.draw(gt) # doesn't work in Rmarkdown, so below function fix the problem
print.ggExtraPlot <- function(x, newpage = is.null(vp), vp = NULL,...) {
  if (newpage) grid::grid.newpage()
  grid::grid.draw(x)
}

pdf("Dorea.phylo1.pdf", width = 3.54, height = 6)
print.ggExtraPlot(gt)
dev.off()
print.ggExtraPlot(gt)

```


#Phylogenetic tree2 (marker based)
```{r}

p <- ggtree(tree) 

## attach the sampling information data set 
## and add symbols colored by location
p <- p %<+% info2 + geom_tippoint(aes(color=Lab), size=1) +
                    geom_tiplab(aes(label=Lab.org, color=Lab), 
                                size=1.5, align = T, fontface="bold") +
                    scale_color_manual(values=cols) +
                    scale_fill_manual(values =cols) +
                    #D.formicigenerans
                    geom_cladelabel(node=95, label="D.formicigenerans", 
                    size=5, offset=-0.02, barsize = NA, 
                    fontsize=5, fontface="italic", angle = 90,
                    color="#377eb8", align=FALSE) +
                    #Dorea.sp
                    geom_cladelabel(node=132, label="Dorea sp.", 
                    size=5, offset=-0.2, barsize = NA, 
                    fontsize=5, fontface="italic", angle = 0,
                    color="#984ea3", align=FALSE) +
                    #D.longicatena
                    geom_cladelabel(node=119, label="D.longicatena", 
                    size=5, offset=-0.05, barsize = NA, 
                    fontsize=5, fontface="italic", angle = 90,
                    color="#006d2c", align=FALSE) #+
                    # A
                    #geom_cladelabel(node=134, label="D.longicatena A", 
                    #size=5, offset=-0.22, barsize = NA, 
                    #fontsize=5, fontface="italic", angle = 0,
                    #color="#006d2c", align=FALSE) +
                    # B
                    #geom_cladelabel(node=153, label="D.longicatena B", 
                    #size=5, offset=-0.22, barsize = NA, 
                    #fontsize=5, fontface="italic", angle = 0,
                    #color="#006d2c", align=FALSE)
  
## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
p2 <- p + geom_facet(panel = "CEGs", data = bar_tbl, 
                     geom = ggstance::geom_barh, 
             aes(x = CEGs, color = Lab, fill = Lab), 
             stat = "identity", width = .6) + xlim_tree(0.25) +
   theme_tree2(legend.position="none") +
   theme(strip.text = element_blank(),
     strip.background.x=element_rect(fill=NA))


gt = ggplot_gtable(ggplot_build(p2))
#gtable_show_layout(gt) # will show you the layout - very handy function
gt # see plot layout in table format
gt$layout$l[grep('panel-1-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[7] = 0.5*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
#grid.draw(gt) # doesn't work in Rmarkdown, so below function fix the problem
print.ggExtraPlot <- function(x, newpage = is.null(vp), vp = NULL,...) {
  if (newpage) grid::grid.newpage()
  grid::grid.draw(x)
}

pdf("Dorea.phylo2.pdf", width = 3.54, height = 6)
print.ggExtraPlot(gt)
dev.off()
print.ggExtraPlot(gt)

##############################################################################
# Just to find parent node information
ggtree(tree, branch.length = "none")+ geom_text(aes(label=node), hjust=-.3)

```


# Dorea.longicatena (core-genome)

```{r}

tree <- read.tree("Dorea_longicatena/core_gene_alignment.newick")

cols <- c("Dorea longicatena"='black', 
          "Dorea longicatena_B"='#006d2c')

info %>% 
  filter(Lab== "Dorea longicatena")%>%
  select(user_genome,Species, Lab,Lab.org) %>%
  mutate(Lab.org= case_when(user_genome=="GC568" ~ paste(Lab.org,
                                                         "***",
                                                         sep = ""),
                            TRUE ~ as.character(Lab.org)))-> info2
bar_tbl <- info %>% select(user_genome, CEGs)


p <- ggtree(tree) 
## attach the sampling information data set 
## and add symbols colored by location
p <- p %<+% info2 + geom_tippoint(aes(color=Species), size=1) +
                    geom_tiplab(aes(label=Lab.org, color=Species), 
                                size=2, align = T)+ #,fontface="bold") +
                    scale_color_manual(values=cols) +
                    scale_fill_manual(values =cols) +
                    #D.longicatena
                    geom_cladelabel(node=89, label="Dorea longicatena A", 
                    size=5, offset=-0.02, barsize = NA, 
                    fontsize=5, fontface="italic", angle = 90,
                    color="black", align=FALSE)+
                    #D.longicatena B
                    geom_cladelabel(node=27, label="D.longicatena B", 
                    size=5, offset=0.03, barsize = NA, 
                    fontsize=5, fontface="italic", angle = 90,
                    color="#006d2c", align=FALSE)+
                    geom_hilight(node=59, fill="#006d2c", 
                                 alpha=0.2, extendto=0.02)
                    
  
  
## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
p2 <- p + geom_facet(panel = "CEGs", data = bar_tbl, 
                     geom = ggstance::geom_barh, 
             aes(x = CEGs, color = Species, fill = Species), 
             stat = "identity", width = .6) + xlim_tree(0.17) +
   theme_tree2(legend.position="none") +
   theme(strip.text = element_blank(), legend.title = element_blank(),
     strip.background.x=element_rect(fill=NA))


gt = ggplot_gtable(ggplot_build(p2))
#gtable_show_layout(gt) # will show you the layout - very handy function
gt # see plot layout in table format
gt$layout$l[grep('panel-1-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[7] = 0.5*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
#grid.draw(gt) # doesn't work in Rmarkdown, so below function fix the problem
print.ggExtraPlot <- function(x, newpage = is.null(vp), vp = NULL,...) {
  if (newpage) grid::grid.newpage()
  grid::grid.draw(x)
}

pdf("Dorea.phylo3.pdf", width = 3.54, height = 4.7)
print.ggExtraPlot(gt)
dev.off()
print.ggExtraPlot(gt)

ggtree(tree, branch.length = "none")+ geom_text(aes(label=node), hjust=-.3)

```













