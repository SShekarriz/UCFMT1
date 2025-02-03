library(tidyverse)

c_inBin <- read.csv("../data/Assembly/contigs_inBins.txt", sep = "\t", header = F)
c_inMAG <- read.csv("../data/Assembly/contigs_inMAGs.txt", sep = "\t", header = F)
colnames(c_inBin) <- c("contig", "bin_id")
colnames(c_inMAG) <- c("contig", "bin_id")

# NB that c_inMAG is a proper subset of c_inBin, so I'm not sure we need both.

# Get a vector of the unique MAG ids.
MAGs = (c_inMAG 
        %>% select(bin_id) 
        %>% distinct() 
        %>% pull())


##### Read in the MAG coverage ####
# define colnames for marker_lvl_mg
heads_cover <- c("sample.id",
                 "contig", "startpos", "endpos", 
                 "numreads", "covbases", "coverage", "meandepth",
                 "meanbaseq", "meanmapq",
                 "ReadSampling", "Quality")

path="../../UCFMT1-old/data/perfect_mapping/DBD/"
patt="_perfect_.*.cover" # these are sample coverage files. Rows are contigs.
marker_lvl_mg = (data.frame(sample.id = paste(dir(path,
                                                  pattern = patt))) 
                 # read each file from the directory (current path) and then unnest
                 %>% mutate(file_contents = map(sample.id, ~ 
                                                    read.csv(file.path(path, .),
                                                             sep = "\t", header = F,
                                                             comment.char = "#"))) 
                 %>% unnest() 
                 %>%  mutate(sample.id = gsub("_perfect_BDB.cover", 
                                              "", sample.id)))
colnames(marker_lvl_mg) <- heads_cover

# Join the coverage info with the bin info, filter to only keep mags, calculate
# normalized coverage
marker_lvl_mg_all = (marker_lvl_mg 
                     %>% left_join(c_inBin) 
                     %>% mutate(bin_q= case_when(bin_id %in% MAGs ~ paste("MAG"),
                                                 TRUE ~ paste("Bin"))) 
                     %>% filter(bin_q == "MAG") 
                     %>% select(-bin_q) 
                     %>% dplyr::rename(sample=sample.id) 
                     %>% group_by(sample, bin_id) 
                     %>% summarise(bin_len= sum(endpos),
                                   bin_covbases=sum(covbases),
                                   bin_numreads=sum(numreads),
                                   bin_coverage=bin_covbases/bin_len*100))
save(MAGs, file = '../data/MAGs.RData')
write.csv(marker_lvl_mg_all, file = '../data/marker_lvl_mg.csv')
