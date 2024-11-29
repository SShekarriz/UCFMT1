#### Functions for sensitivity analysis. Requires permutation_functions.R.

sens_fig_13_df = function(arr_cts, txrm){
    
    # Check the inputs
    if (!(txrm %in% c('Treatment', 'Remission'))){
        stop('txrm must be one of "Treatment" or "Remission".')
    }
    
    # Set up to iterate over all the presence/absence cutoff values
    n = dim(arr_cts)[1]
    
    # Extract the absence/presence cutoff values
    reps = dimnames(arr_cts)$replicate
    tmp = get_reps(reps)
    abs = tmp[['abs']]
    pres = tmp[['pres']]
    
    # Create the base of the dataframe
    
    plt_all = data.frame(Npats = NA, Treatment = NA, Count = NA, Replicate = NA,
                         AbsCutoff = NA, PresCutoff = NA)
    for (i in 1:n){
        cts = arr_cts[i,,]
        cts = cts[rowSums(cts)>0,] 
        cts_tx = table(cts[,1])
        cts_bl = table(cts[,2])
        
        mtx = as.numeric(names(cts_tx)[length(cts_tx)])
        mbl = as.numeric(names(cts_bl)[length(cts_bl)])
        maxlen = max(mtx, mbl) + 1
        plt_mat = matrix(nrow = maxlen, ncol = 2)
        rownames(plt_mat) = as.character(0:(maxlen - 1))
        if (txrm == 'Treatment'){
            colnames(plt_mat) = c('FMT', 'Placebo')
        } else {
            colnames(plt_mat) = c('Res', 'NoRes')
        }
        plt_mat[names(cts_tx),1] = cts_tx
        plt_mat[names(cts_bl),2] = cts_bl
        plt_mat = plt_mat[-1,]
        plt_mat[is.na(plt_mat)] = 0
        
        plt_df = (plt_mat
                  %>% data.frame()
                  %>% rownames_to_column('Npats')
                  %>% mutate(Npats = as.numeric(Npats))
                  %>% pivot_longer(names_to = txrm, values_to = 'Count',
                                   -Npats)
                  %>% mutate(Replicate = dimnames(arr_cts)$replicate[i],
                             AbsCutoff = abs[i],
                             PresCutoff = pres[i]))
        plt_all = rbind(plt_all, plt_df)
    }
    plt_all = plt_all[-1,]
    plt_all = (plt_all
               %>% mutate(AbsCutoff = factor(AbsCutoff, levels = abs[1:9],
                                             labels = paste('abs:',abs[1:9])),
                          PresCutoff = factor(PresCutoff, 
                                              levels = unique(pres),
                                              labels = paste('pres:',
                                                             unique(pres)))))
    return(plt_all)
}


sens_plot_fig_13 = function(plt_all, txrm){
    
    f1 = ggplot(plt_all, aes(Npats, Count, colour = .data[[txrm]])) +
        geom_point() +
        geom_line() +
        facet_grid(PresCutoff ~ AbsCutoff) +
        scale_x_continuous(breaks = 1:max(plt_all$Npats))
    return(f1)
}

sens_fig_2_df = function(arr_cts){
    
    n = dim(arr_cts)[1]
    reps = dimnames(arr_cts)$replicate
    tmp = get_reps(reps)
    abs = tmp[['abs']]
    pres = tmp[['pres']]
    
    cts_all = data.frame(Total = NA, Uniqueness = NA, names = NA, Count = NA, 
                         Replicate = NA, AbsCutoff = NA, PresCutoff = NA)
    for (i in 1:n){
        cts = arr_cts[i,,]
        cts_df = (data.frame(cts)
                  %>% mutate(Total = rowSums(.),
                             Uniqueness = case_when(tx_tot > 0 & bl_tot > 0 ~
                                                        'Both',
                                                    tx_tot == 0 & bl_tot > 0 ~
                                                        'Placebo',
                                                    tx_tot > 0 & bl_tot == 0 ~
                                                        'FMT',
                                                    tx_tot == 0 & bl_tot == 0 ~
                                                        'None'))
                  %>% filter(Uniqueness != 'None')
                  %>% pivot_longer(tx_tot:bl_tot, names_to = 'names',
                                   values_to = 'Count')
                  %>% filter(Count > 0)
                  %>% mutate(Replicate = reps[i],
                             AbsCutoff = abs[i],
                             PresCutoff = pres[i]))
        # this still counts the Both ones twice ^
        # Get each both once:
        cts_both = (cts_df
                    %>% filter(Uniqueness == 'Both', names == 'tx_tot'))
        # Remove Both from the original
        cts_df = (cts_df
                  %>% filter(Uniqueness != 'Both'))
        cts_df = rbind(cts_df, cts_both)
        
        cts_all = rbind(cts_all, cts_df)
        
    }
    cts_all = cts_all[-1,]
    return(cts_all)
}
    

sens_plot_fig2 = function(cts_all){    
    
    f2 = ggplot(cts_all, aes(x = Total, fill = Uniqueness,
                                        colour = Uniqueness)) +
        geom_bar(stat = 'count') +
        scale_x_continuous(breaks = 1:max(cts_all$Total)) +
        facet_grid(PresCutoff ~ AbsCutoff)
    return(f2)
}

get_reps = function(reps){
    tmp = unlist(str_split(reps, '_'))
    abs = tmp[seq(2, length(tmp), by = 4)]
    pres = tmp[seq(4, length(tmp), by = 4)]
    return(list(abs = abs, pres = pres))
}

rep_perms = function(engr_arr, arr_cts, nperm,
                     tx_pv_mgm, txrm){
    n = dim(engr_arr)[1]
    
    stat_arr = array(dim = c(n, nperm, 3), 
                     dimnames = list(replicate = dimnames(arr_cts)$replicate,
                                     permutations = 1:nperm,
                                     test_stat = c('fx','xfx','x2fx')))
    pval_arr = matrix(nrow = n, ncol = 3)
    for (i in 1:n){
        print(i)
        engr = engr_arr[i,,]
        cts = arr_cts[i,,]
        obs = get_stats(cts)
        perms = do_permute(engr, cts, obs, tx_pv_mgm, nperm, txrm = txrm)
        stat_arr[i,,] = perms[['stat_mat']]
        pval_arr[i,] = round(get_pvals(perms[['stat_mat']]), 5)
    }
    return(list(stat_arr = stat_arr, pval_arr = pval_arr))   
}

sens_perm_df = function(stat_arr){
    
    n = dim(stat_arr)[1]
    reps = dimnames(arr_cts)$replicate
    tmp = get_reps(reps)
    abs = tmp[['abs']]
    pres = tmp[['pres']]
    
    sm_all = data.frame(test_stat = NA, values = NA, obs_val = NA,
                        Replicate = NA, AbsCutoff = NA, PresCutoff = NA)
    for (i in 1:n){
        sm_df = (as.data.frame(stat_mat)
                 %>% pivot_longer(everything(),names_to = 'test_stat',
                                  values_to = 'values')
                 %>% mutate(test_stat = factor(test_stat, 
                                               levels = c('fx','xfx','x2fx')))
                 %>% arrange(test_stat)
                 %>% mutate(obs_val = rep(stat_mat[1,],
                                          each = nrow(stat_mat)),
                            Replicate = reps[i],
                            AbsCutoff = abs[i],
                            PresCutoff = pres[i]))
        sm_all = rbind(sm_all, sm_df)
    }
    sm_all = sm_all[-1,]
    
    return(sm_all)
}

sens_perm_plot = function(sm_all, test){
    if (!(test %in% c('fx','xfx','x2fx'))){
        msg = "'test' must be 'fx', 'xfx', or 'x2fx'."
        stop(msg)
    }
    
    sm_use = filter(sm_all, test_stat == test)
    pl_pv = ggplot(sm_use, aes(x = values)) +
        geom_histogram(fill = 'grey59') +
        geom_vline(aes(xintercept = obs_val),
                   colour = 'violetred4',
                   linetype = 2) +
        facet_grid(PresCutoff~AbsCutoff)
    return(pl_pv)
}
