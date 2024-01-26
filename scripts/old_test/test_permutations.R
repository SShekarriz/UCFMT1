load('../permut_data/intermed/perms_ge_rem.RData')
names(perms_ge_rem)
dim(perms_ge_rem[['cts_array']])
cts_array = perms_ge_rem[['cts_array']]
table(cts_array[,1,1])
table(cts_array[,2,1])

tab_mat = matrix(nrow = 13, ncol = 2000)
rownames(tab_mat) = as.character(0:12)
colnames(tab_mat) = paste('it', 1:2000, sep = '_')
tab_mat[1:10,1:10]
tst = table(cts_array[,1,1])
names(tst)
tab_mat[,1][names(tst)] = tst

tab_tx = matrix(nrow = 13, ncol = 2000)
tab_bl = matrix(nrow = 13, ncol = 2000)
rownames(tab_tx) = rownames(tab_bl) = as.character(0:12)
colnames(tab_tx) = colnames(tab_bl) = paste('it', 1:2000, sep = '_')

for (i in 1:2000){
    tx = table(cts_array[,1,i])
    bl = table(cts_array[,2,i])
    tab_tx[names(tx),i] = tx
    tab_bl[names(bl),i] = bl
}

tab_tx[1:10,1:10]

library(tidyverse)
tab_df_tx = (tab_tx
          %>% data.frame()
          %>% rownames_to_column('Npats')
          %>% pivot_longer(names_to = 'Iterations', values_to = 'Counts',
                           -Npats)
          %>% mutate(Treatment = 'FMT'))
head(tab_df_tx)

tab_df_bl = (tab_bl
          %>% data.frame()
          %>% rownames_to_column('Npats')
          %>% pivot_longer(names_to = 'Iterations', values_to = 'Counts',
                           -Npats)
          %>% mutate(Treatment = 'Placebo'))

tab_df = rbind(tab_df_tx, tab_df_bl)
tab_df 
head(tab_df)
tail(tab_df)

tab_df = filter(tab_df, !is.na(Counts), Npats != '0')

ggplot(filter(tab_df, Iterations %in% paste('it', 1, sep = '_')),
       aes(x = as.numeric(Npats), y = Counts, colour = Treatment)) +
    geom_point() +
    geom_line() +
    facet_wrap(~Iterations) +
    scale_colour_brewer(palette = 'Dark2') +
    scale_x_continuous(breaks = 1:6) +
    theme_bw()
