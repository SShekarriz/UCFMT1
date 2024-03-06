#### benchmarking functions ####

# a function to filter marker above cutoff based on relative abundance
benchmark1 <- function(table, cutoff){
    # total number of unqiue species in each sample
    table %>%
        group_by(sample) %>%
        summarise(total=n()) -> total.marker
    
    # total number of unique marker above cutoff
    cutoff_tab = (table %>%
                      mutate(relative_abund= abundance) %>%
                      filter(relative_abund >= cutoff) %>%
                      group_by(sample) %>%
                      summarise(pos_n=n()) %>%
                      left_join(total.marker) %>%
                      mutate(pos_p=pos_n/total * 100) %>%
                      mutate(cutoff=paste(cutoff)))
    
    return(cutoff_tab)
}

# comparison of detection cut-offs
benchmark <- function(table, vector, f_type, bs_4 = NULL) {
    
    vector <- as.character(vector)  # Convert vector to character for proper factor levels
    
    cols <- c("2013" = "#bae4b3",
              "2016" = "#74c476", 
              "May17" = "#31a354",
              "Oct17" = "#006d2c")
    
    # comparison of detection cut-offs
    res <- list()
    for (i in vector) {
        res[[i]] <- benchmark1(table, as.numeric(i))
    }
    
    tbl <- do.call(rbind.data.frame, res)
    tbl$cutoff <- factor(tbl$cutoff, levels = vector)
    
    f1 <- ggplot(tbl, aes(cutoff, pos_p)) 
    if(f_type == "16s"){
        f1 = f1 +
            geom_line(aes(group = sample), colour = "#006d2c")
    } else {
        f1 = f1 +
            geom_line(aes(group = sample, colour = sample)) +
            scale_colour_manual(values = cols)
    }
    f1 = f1 +
        theme_bw() +
        theme(legend.position = c(0.8,0.7),
              axis.title = element_text(size=7),
              legend.title = element_blank()) +
        xlab("relative abundance cutoff (0 - 1)") +
        ylab("proportion of unique markers kept")
    
    
    table %>%
        ungroup() %>%
        group_by(sample) %>%
        summarise(percent=as.numeric(abundance)
                  /sum(as.numeric(abundance))) -> hist_table
    
    # get the 20% quantile of the data
    q05 <- quantile(hist_table$percent, 0.20)
    q05 = formatC(q05, format = 'e', digits = 1)
    
    if(!is.null(bs_4)){
        hist_table = (hist_table
                      %>% filter(sample %in% bs_4))
    }
    
    f2 <- ggplot(hist_table, aes(x=percent)) +
        geom_histogram(position="identity", colour="grey40", fill="#31a354",
                       alpha=0.7, bins = 10) +
        scale_x_log10() +
        facet_grid(. ~ sample) +
        geomtextpath::geom_textvline(label = q05, color = "black", 
                                     xintercept = as.numeric(q05), vjust = 1.3) +
        theme_classic() +
        theme(legend.position = "none",
              axis.title = element_text(size=7),
              axis.text.x = element_text(angle = 90)) +
        xlab("relative abundance (0 - 1)") +
        ylab("# of markers present at a given relative abundance")
    
    fig = cowplot::plot_grid(f1, f2, ncol = 1)
    
    return(fig)
    
}

# a function to filter marker above cutoff based on genomic coverage
benchmark2 <- function(table, cutoff){
    # total number of unqiue species in each sample
    table %>%
        group_by(sample) %>%
        summarise(total=n()) -> total.marker
    
    # total number of unique marker above cutoff
    table = (table %>%
                 filter(coverage >= cutoff) %>%
                 group_by(sample) %>%
                 summarise(pos_n=n()) %>%
                 left_join(total.marker) %>%
                 mutate(pos_p=pos_n/total * 100) %>%
                 mutate(cutoff=paste(cutoff)))
    
    return(table)
}

# comparison of detection cut-offs
Genome_benchmark <- function(table, vector, feature = 'UHGG') {
    
    cols <- c("2013" = "#bae4b3",
              "2016" = "#74c476", 
              "May17" = "#31a354",
              "Oct17" = "#006d2c")
    
    
    vector <- as.character(vector)  # Convert vector to character for proper factor levels
    
    res <- list()
    for (i in vector) {
        res[[i]] <- benchmark2(table, as.numeric(i))
    }
    
    tbl <- do.call(rbind.data.frame, res)
    tbl$cutoff <- factor(tbl$cutoff, levels = vector)
    
    f1 <- ggplot(tbl, aes(cutoff, pos_p, color = sample, fill=sample)) +
        geom_line(aes(group = sample)) +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        theme_bw() +
        theme(legend.position = c(0.8,0.7),
              axis.title = element_text(size=7),
              legend.title = element_blank()) +
        xlab("1x genome coverage cutoff (%)") +
        ylab("proportion of unique genome kept")
    
    table %>%
        select(sample, coverage) -> hist_table
    
    # get the 20% quantile of the data
    # q_top <- quantile(hist_table$coverage, 0.5)
    # q_top = formatC(q_top, format = 'e', digits = 1)
    # q_bottom <- quantile(hist_table$coverage, 0.27)
    # q_bottom = formatC(q_bottom, format = 'e', digits = 2)
    
    
    f2 <- ggplot(hist_table, aes(x=coverage, fill = sample, color= sample)) +
        geom_histogram(position="identity", colour="grey40", alpha=0.7, bins = 10) +
        #scale_x_log10() +
        facet_grid(. ~ sample) +
        geomtextpath::geom_textvline(label = '0.75', color = "black", 
                                     xintercept = 0.75, vjust = 1.3) +
        scale_fill_manual(values = cols) +
        theme_classic() +
        theme(legend.position = "none",
              axis.title = element_text(size=7),
              axis.text.x = element_text(angle = 90)) +
        xlab("1x genome coverage (%)") +
        ylab("# of genome present at a given coverage")
    
    if(feature == 'MAG'){
        f2 = f2 +
            geomtextpath::geom_textvline(label = '0.20', color = "black",
                                         xintercept = 0.20, 
                                         vjust = 1.3)
    }
    
    
    fig = cowplot::plot_grid(f1, f2, ncol = 1)
    return(fig)
}

#### f_paper functions ####

f_paper <- function(f1, f2, f3, tx_pv, rs_pv, pvals) {
    
    # f1, f2, and f3 are the ggplot object 
    # generated by plot_fig13() and plot_fig2()
    # tx_pv and rs_pv named vectors of patient IDs/groups
    #pvals a datafame contains labels for f1 and f3 figures
    
    #pull out ggplot data
    data1 <- ggplot_build(f1)$plot$data
    
    fig1_cols <- c("FMT" = "#1b7837",
                   "Placebo" = "#999999")
    
    p1 = f_paper_p1(data1, tx_pv, f1, fig1_cols, pvals)
    
    fig2_cols <- c("FMT" = "#1b7837",
                   "Placebo" = "#999999",
                   "Both" = "#fc8d59")
    p2 = f_paper_p2(data1, f2, tx_pv, fig2_cols)
    
    fig3_cols <- c("Res" = "#2166ac",
                   "NoRes" = "#b2182b")
    
    p3 = f_paper_p3(rs_pv, f3, fig3_cols, pvals)
    
    fig = cowplot::plot_grid(p1, p2, p3, ncol = 1)
    return(fig)
}

f_paper_p1 = function(data1, tx_pv, f1, fig1_cols, pvals){
    
    data.frame(Treatment = tx_pv) %>%
        group_by(Treatment) %>%
        tally() %>%
        mutate(annot= paste(Treatment, "(n=", n, ")", sep = "")) %>%
        select(Treatment, annot)-> i_tx
    # re-label the color variable that include n
    labels <- setNames(i_tx$annot, i_tx$Treatment)
    # pvalue-label
    p_label = pvals %>% filter(plot == "f1") %>% pull(label)
    p1 = f1 +
        scale_color_manual(values = fig1_cols,
                           labels = labels) +
        annotate("text", 
                 x = max(data1$Npats) * 0.85, y = max(data1$Count) * 0.4, 
                 label = p_label, size = 2.5,
                 colour = "black", parse = FALSE) +
        theme_classic() +
        theme(legend.position = c(0.75, 0.8),
              legend.title = element_blank(),
              axis.title = element_blank())
    
    return(p1)
}

f_paper_p2 = function(data1, f2, tx_pv, fig2_cols){
    #pull out ggplot data
    data2 <- ggplot_build(f2)$plot$data
    # a label for second figure, showing n
    labels <- sprintf("FMT+Placebo\n(n=%d)", length(tx_pv))
    
    p2 = f2 +
        scale_fill_manual(values = fig2_cols) +
        scale_colour_manual(values = fig2_cols) +
        annotate("text", 
                 x = max(data2$Total) * 0.8, y = max(data1$Count) * 0.8, 
                 label = labels, size = 3,
                 colour = "black", parse = FALSE) +
        theme_classic() +
        theme(legend.position = "none",
              legend.title = element_blank(),
              axis.title = element_blank())
    
    return(p2)
}

f_paper_p3 = function(rs_pv, f3, fig3_cols, pvals){
    
    data.frame(Remission = rs_pv) %>%
        group_by(Remission) %>%
        tally() %>%
        mutate(annot= paste(Remission, "(n=", n, ")", sep = "")) %>%
        select(Remission, annot)-> i_rs
    #pull out ggplot data
    data3 <- ggplot_build(f3)$plot$data
    # re-label the color variable that include n
    labels <- setNames(i_rs$annot, i_rs$Remission)
    # pvalue-label
    p_label = pvals %>% filter(plot == "f3") %>% pull(label)
    p3 = f3 +
        scale_colour_manual(values = fig3_cols,
                            labels = labels) +
        annotate("text", 
                 x = max(data3$Npats) * 0.85, y = max(data3$Count) * 0.4, 
                 label = p_label, size = 2.5,
                 colour = "black", parse = FALSE) +
        theme_classic() +
        theme(legend.position = c(0.75, 0.8),
              legend.title = element_blank(),
              axis.title = element_blank())
    
    return(p3)
}

#### p_paper functions ####


p_paper <- function(profile_pat, profile_don, mapfile) {
    
    cols <- c("Engraft" = "#1a9641",
              "Other" = "#a6d96a",
              "Shared" = "#2b83ba",
              "Donor" = "#abd9e9")
    
    # profile number of markers in each Treatment
    p1 = p_paper_p1(profile_pat, mapfile)
    
    f1 = p_paper_f1(p1, cols)
    
    p2 = p_paper_p2(p1, profile_pat)
    
    f2 <- p_paper_f2(p2, cols)
    f3 <- p_paper_f3(profile_don, cols)
    fig = cowplot::plot_grid(f1, f2, f3, nrow = 1)
    return(fig)
}

p_paper_p1 = function(profile_pat, mapfile){
    
    mapfile %>%
        select(Fig_lab, Treatment, Remission) %>%
        distinct()-> map
    
    profile_pat %>%
        rownames_to_column("Marker") %>%
        gather(Fig_lab, engrafted, -Marker) %>%
        group_by(Fig_lab, engrafted) %>%
        summarise(Marker_n= n()) %>%
        #mutate(proportion= Marker_n / donor_markers * 100) %>%
        left_join(map) %>%
        mutate(Fig_lab = case_when(Remission == "Res" ~ paste("*", Fig_lab,
                                                              sep = ""),
                                   TRUE ~ paste(Fig_lab)),
               engrafted = factor(engrafted, levels = rev(c('Donor',
                                                            'Shared',
                                                            'Engraft',
                                                            'Other')))) -> p1
    return(p1)
}

p_paper_p2 = function(p1, profile_pat){
    
    donor_markers = nrow(profile_pat)
    p1 %>%
        ungroup() %>%
        group_by(Treatment, engrafted) %>%
        summarise(Mean_marker= mean(Marker_n) / donor_markers ) %>%
        mutate(engrafted = factor(engrafted, levels = rev(c('Donor',
                                                            'Shared',
                                                            'Engraft',
                                                            'Other')))) -> p2
    return(p2)
}

p_paper_f1 = function(p1, cols){
    f1 <- ggplot(p1, aes(Marker_n, Fig_lab, fill=engrafted)) +
        geom_bar(stat = "identity") +
        facet_grid(Treatment~., scales = "free", space = "free") +
        scale_fill_manual(values = cols) +
        theme_classic() +
        theme(legend.position = "none") +
        xlab("# of donor features") +
        ylab("Patients")
    
    return(f1)
}

p_paper_f2 = function(p2, cols){
    
    f2 <- ggplot(p2, aes(Treatment, Mean_marker, fill=engrafted)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        geom_text(aes(label = format(Mean_marker, format = 'e', digits = 1)),
                  position = position_stack(vjust = .5)) +
        theme_classic() +
        theme(legend.position = "rig") +
        ylab("Mean proportion across all patients (%)")
    return(f2)
}

p_paper_f3 = function(profile_don, cols){
    
    f3 <- ggplot(profile_don, aes(x = status)) +
        geom_bar(aes(fill = status), position = position_dodge()) +
        scale_fill_manual(values = cols) +
        theme_classic() +
        theme(legend.position = "none") +
        facet_grid(Treatment~.) +
        ylab('# of donor features')
    return(f3)
}

#### profile functions ####

## Profile the type of donor Marker in patients:
# 1. Engraft, 2. UniqueToDonor, 3. SharedWithDonor, 4. Other events
get_profile_pat = function(long_markerlvl, mapfile, cutoff_abs, cutoff_pres,  ...){
    
    tst = (long_markerlvl 
           %>% prep_profile(mapfile, cutoff_abs, cutoff_pres, ...)
           %>% select(Marker, Fig_lab,  engrafted)
           %>% pivot_wider(names_from = Fig_lab, values_from = engrafted,
                           values_fill = NA)
           # Make the Marker column the rownames
           %>% column_to_rownames('Marker'))
    return(tst)
}

get_profile_don = function(long_markerlvl, mapfile, cutoff_abs, cutoff_pres,
                           ...){
    tst = (long_markerlvl
           %>% prep_profile(mapfile, cutoff_abs, cutoff_pres, 
                            Treatment %in% c('FMT', 'Placebo'))
           %>% group_by(Marker, Treatment)
           %>% summarize(status = get_status(engrafted))
    )
    shr = (tst
           %>% filter(status == 'Both')
           %>% mutate(status = 'Shared'))
    tst = (tst
           %>% mutate(status = case_when(status == 'Both' ~ 'Engraft',
                                         TRUE ~ status))
           %>% rbind(shr)
           %>% mutate(status = factor(status, levels = c('Donor',
                                                         'Shared',
                                                         'Engraft','Other'))))
    return(tst)
}

get_status = function(engrafted){
    status = case_when(all(engrafted == 'Donor') ~ 'Donor',
                       any(engrafted == 'Engraft') & 
                           any(engrafted == 'Shared') ~ 'Both',
                       any(engrafted == 'Engraft') ~ 'Engraft',
                       any(engrafted == 'Shared') ~ 'Shared',
                       TRUE ~ 'Other')
    return(status)
}

prep_profile = function(long_markerlvl, mapfile, cutoff_abs, cutoff_pres, ...){
    tst = (long_markerlvl
           %>% left_join(mapfile, by = c('sample' = 'Study_ID'))
           %>% filter(Treatment %in% c('FMT', 'Placebo'))
           %>% select(Marker, Timepoint, abundance, Fig_lab, Treatment)
           %>% spread(Timepoint, abundance, fill = 0)
           %>% mutate(engrafted = 
                          case_when(WK0 <= cutoff_abs & WK6 >= cutoff_pres ~ 
                                        "Engraft",
                                    WK0 <= cutoff_abs & WK6 <= cutoff_abs  ~ 
                                        "Donor",
                                    WK0 >= cutoff_pres ~ "Shared",
                                    TRUE ~ "Other")))
    return(tst)
}