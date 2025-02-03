library(tidyverse)
library(AfterSl1p)

tax_16s = read.csv('./data/featTaxonomy/engrafted_16s_lineage.csv', 
                   row.names = 1)
tax_sp = read.csv('./results/engrafted_sp_taxa.csv', row.names = 1)
tax_st = read.csv('./data/featTaxonomy/engrafted_st_lineage.csv',
                  row.names = 1)
tax_mg = read.csv('./data/featTaxonomy/engrafted_mg_lineage.csv',
                  row.names = 1)

all_tax = (tax_16s
           %>% mutate(Strain = NA,
                      FeatureType = '16S')
           %>% select(FeatureType, everything()))

tax_sp_sep = (tax_sp
              %>% mutate(Kingdom = TaxInfo)
              %>% separate_wider_delim(Kingdom, delim = '|', 
                                       names = c('Kingdom', 'Phylum', 'Class',
                                                 'Order', 'Family', 'Genus',
                                                 'Species', 'Strain'),
                                       too_few = 'align_start',
                                       too_many = 'error')
              %>% mutate(FeatureType = 'SpeciesMarkers')
              %>% select(FeatureType, everything()))
tax_st_sep = (tax_st
              %>% unite(col = 'TaxInfo', SGB, st_marker)
              %>% separate_wider_delim(lineage, delim = '|',
                                       names = c('Kingdom', 'Phylum', 'Class',
                                                 'Order', 'Family', 'Genus',
                                                 'Species', 'Strain'),
                                       too_few = 'align_start',
                                       too_many = 'error')
              %>% mutate(FeatureType = 'StrainMarkers')
              %>% select(FeatureType, everything()))
tax_mg_sep = (tax_mg
              %>% separate_wider_delim(classification, delim = ';',
                                       names = c('Kingdom', 'Phylum', 'Class',
                                                 'Order', 'Family', 'Genus',
                                                 'Species', 'Strain'),
                                       too_few = 'align_start',
                                       too_many = 'error')
              %>% mutate(FeatureType = 'MAGs')
              %>% select(FeatureType, everything()))

all_tax = rbind(all_tax, tax_sp_sep, tax_st_sep, tax_mg_sep)

fix_tax = function(x){
  new_x = str_replace(x, '^[a-z]__', '')
  new_x[new_x == ''] = NA
  return(new_x)
}

all_tax_clean = (all_tax
           %>% select(-Strain, -TaxInfo)
           %>% mutate(across(Kingdom:Species, fix_tax))
           %>% filter(!is.na(Kingdom)))

tax_propped = AfterSl1p:::prop_tax_tab(select(all_tax_clean, Kingdom:Species),
                                       indic = FALSE)
all_tax_clean[,3:9] = tax_propped

all_tax_gl_sp = (all_tax_clean
                   %>% group_by(FeatureType, Kingdom, Phylum, Class, Order,
                                Family, Genus, Species)
                   %>% summarize(nPat = sum(nPat),
                                 .groups = 'drop')
                   %>% group_by(FeatureType)
                   %>% mutate(relPat = nPat/sum(nPat)))

all_tax_gl_ph = (all_tax_gl_sp
                 %>% group_by(FeatureType, Phylum)
                 %>% summarize(nPat = sum(nPat), .groups = 'drop')
                 %>% group_by(FeatureType)
                 %>% mutate(relPat = nPat/sum(nPat)))

# Make a data frame for plotting at family level

top_fam_df = (all_tax_gl_sp
              %>% group_by(FeatureType, Family)
              %>% summarize(nPat = sum(nPat), .groups = 'drop')
              %>% group_by(FeatureType)
              %>% mutate(relPat = nPat/sum(nPat))
              %>% group_by(FeatureType)
              %>% arrange(desc(relPat), .by_group = TRUE)
              %>% group_by(FeatureType)
              %>% slice_head(n = 10))


fam_plt = ggplot(top_fam_df, aes(x = Family, y = relPat)) +
  geom_bar(stat = 'identity', aes(colour = FeatureType, fill = FeatureType),
           position = position_dodge(width = 0.6),
           width = 0.5) +
  scale_colour_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  ylab('Proportion of Feature Type') +
  theme_classic() +
  rotate_ticks()
fam_plt
ggsave('plots_paper/suppl_taxonomy.png', width = 12, height = 8)





