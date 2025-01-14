library(tidyverse)

map <- "data/UCFMT1_METAmap.txt"
precom_thresh = "data/VallesColomerM_2022_Jan21_thresholds.tsv"
path = "strainphlan"
patt = "_nGD.tsv"

pwdis = (data.frame(marker.id = paste(dir(path, pattern = patt))) 
         %>% mutate(file_contents = map(marker.id, ~ 
                                            read_tsv(file.path(path, .),
                                                     col_names = F))) 
         %>% unnest() 
         %>% mutate(marker.id = gsub(patt, "", marker.id)))

length(unique(pwdis$marker.id))

donB_pwdis = (pwdis
              # only samples compared to donorB:
              %>% filter(str_detect(X1, "DonorB_D")) 
              # not interested in donorB to donorB comparison for now:
              %>% filter(!str_detect(X2, "DonorB")) 
              # I don't think the pairwise column was doing anything
              %>% mutate(distance= X3,
                         sample= gsub(".mp", "", X2),
                         SGB= marker.id,
                         Donor_S= gsub(".mp", "", X1))
              # Reduce it to one row per SGB/Sample combo by taking the smallest
              # distance value
              %>% summarise(distance = min(distance, na.rm = TRUE),
                            .by = c(SGB, sample))
              %>% ungroup()
              # Now we don't need the donor sample column anymore. These are all
              # comparisons to _some_ donor sample, which is all that matters.
              %>% select(SGB, sample, distance))

length(unique(donB_pwdis$SGB))

thresholds <- read.csv(precom_thresh, sep = "\t")
mapfile <- read.csv(map, sep = "\t")

events = (donB_pwdis 
       %>% mutate(Study_ID=gsub("_.*", "", sample)) 
       %>% left_join(mapfile %>% select(Study_ID, Fig_lab, Timepoint)) 
       # Need to remove Study_ID from select or it fills in a bunch of
       # extraneous 1s in the spread
       %>% select(SGB, Fig_lab, Timepoint, distance)
       %>% spread(Timepoint, distance, fill = 1)
       %>% left_join(thresholds) 
       # This way there is only one row per SGB/Fig Lab combo now
       %>% mutate(event= case_when(WK0 > used_nGD_score_percentile &
                          WK6 < used_nGD_score_percentile ~ paste("Engraft"),
                          WK0 > used_nGD_score_percentile &
                          WK6 > used_nGD_score_percentile ~ paste("Donor"),
                          WK0 < used_nGD_score_percentile &
                          WK6 < used_nGD_score_percentile ~ paste("Shared"),
                          is.na(used_nGD_score_percentile) ~ paste("no-cutoff"),
                          TRUE ~ paste("other"))))
eng_events = (events 
              %>% filter(event == "Engraft")
              %>% left_join((mapfile 
                             %>% select(Fig_lab, Treatment, Remission)
                             # Each line is there twice because of the timepoints
                             %>% unique()))
              %>% select(SGB, Fig_lab, Treatment, Remission))

# total number of unique engrafted SGB with strainphlan
length(unique(eng_events$SGB))

results = (eng_events 
           %>% group_by(SGB, Treatment)
           # I think what we're doing here is counting how much engraftment
           # there is for each SGB within each treatment. So we have the number
           # of patients in each treatment group for which that SGB was
           # engrafted.
           %>% summarise(patients = n_distinct(Fig_lab), .groups = "drop")
           %>% group_by(Treatment, patients)
           # Now we've counted how many SGBs there are within group for each
           # number of engrafted patients
           %>% summarise(Engraftment = n(), .groups = "drop"))

engr_sgbs = (eng_events
             %>% count(SGB, Treatment, Remission))
head(engr_sgbs)
write.csv(engr_sgbs, file = 'results/engrafted_sgbs.csv')


fig1 = ggplot(results, aes(patients, Engraftment, color=Treatment)) +
  geom_line() + geom_point() +
  scale_color_manual(values = c("FMT" = "#1a9641",
                               "Placebo" = "#bababa")) +
  theme_bw() +
  theme()

results_rem = (eng_events 
               %>% filter(Treatment == "FMT") 
               %>% group_by(SGB, Remission) 
               %>% summarise(patients = n_distinct(Fig_lab), .groups = "drop") 
               %>% group_by(Remission, patients) 
               %>% summarise(Engraftment = n(), .groups = "drop"))

fig3 = ggplot(results_rem, aes(patients, Engraftment, color=Remission)) +
  geom_line() + geom_point() +
  scale_color_manual(values = c("Res" = "#2b83ba",
                               "NoRes" = "#d7191c")) +
  theme_bw() +
  theme()

head(eng_events)

tst = rbind((eng_events
             %>% select(SGB, Treatment)
             %>% unique()
             %>% filter(Treatment == 'FMT')),
             (eng_events
              %>% select(SGB, Treatment)
              %>% unique()
              %>% filter(Treatment == 'Placebo')))
both = (tst
        %>% count(SGB)
        %>% filter(n > 1))

eng_by_cat = (tst
              %>% mutate(Category = case_when(SGB %in% both$SGB ~ 'Both',
                                               TRUE ~ Treatment))
              %>% select(-Treatment)
              %>% unique())
fig2_cols <- c("FMT" = "#1b7837",
                   "Placebo" = "#999999",
                   "Both" = "#fc8d59")
tst1 = (eng_events
        %>% dplyr::count(SGB)
        %>% left_join(eng_events)
        %>% left_join(eng_by_cat))
fig2 = ggplot(tst1, aes(n, colour = Category, fill = Category)) +
    geom_bar(stat = 'count') +
    scale_fill_manual(values = fig2_cols) +
    scale_colour_manual(values = fig2_cols) +
    theme_classic()
