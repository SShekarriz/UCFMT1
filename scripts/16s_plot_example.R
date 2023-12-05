library(tidyverse)

load('../permut_data/intermed/perms_16s.RData')

# For the 16s data, we should probably plot the averaged counts from the
# subsampled data, and not the full counts from the full data.

cts_16s = perms_16s[['cts_array']][,,1]
head(cts_16s)

# These averages are not integers, so I'm going to bin them by rounding. I can't
# think of a better way to do that.

cts_16s = round(cts_16s)
head(cts_16s)
(fmt_16s = table(cts_16s[,1]))
(pl_16s = table(cts_16s[,2]))

df_16s = data.frame('Npat' = 1:8,
                    'FMT' = as.vector(fmt_16s[-1]),
                    'Placebo' = c(pl_16s[2:6],0,pl_16s[7],0))
df_16s = (df_16s
          %>% pivot_longer(names_to = 'Treatment', values_to = 'Count', -Npat))
ggplot(df_16s, aes(x = Npat, y = Count, colour = Treatment)) +
    geom_point() +
    geom_line() +
    scale_colour_brewer(palette = 'Dark2') +
    theme_bw()

# Do it all again for responders:

load('../permut_data/intermed/perms_16s_rem.RData')
cts_16s_rem = perms_16s_rem[['cts_array']][,,1]

cts_16s_rem = round(cts_16s_rem)

(fmt_16s_rem = table(cts_16s_rem[,1]))
(pl_16s_rem = table(cts_16s_rem[,2]))

df_16s_rem = data.frame('Npat' = 1:3,
                    'FMT' = as.vector(fmt_16s_rem[-1]),
                    'Placebo' = as.vector(pl_16s_rem[-1]))
df_16s_rem = (df_16s_rem
          %>% pivot_longer(names_to = 'Treatment', values_to = 'Count', -Npat))
ggplot(df_16s_rem, aes(x = Npat, y = Count, colour = Treatment)) +
    geom_point() +
    geom_line() +
    scale_colour_brewer(palette = 'Dark2') +
    theme_bw()
