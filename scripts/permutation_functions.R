# New Permutation Functions

## Functions
# Make the marker data frame long
mark_to_long = function(marker_lvl, donB = NULL) {
    
    # For some feature types there is not a resticted list of donor B features
    # because all features are from donor B
    if(!is.null(donB)){
        # But if there is such a list, remove any features that are absent from
        # it
        marker_lvl = filter(marker_lvl, Marker %in% donB)
    } 
    
    # Gather the marker data frame to a long format and fix the naming problems
    # in the sample names
    long_markerlvl = (marker_lvl
                      %>% gather(sample_orig, abundance, -Marker)
                 %>% mutate(sample = case_when(str_detect(sample_orig,"PMCL") ~ 
                                             paste(gsub("_.*", "",sample_orig)),
                                                    TRUE ~ paste(sample_orig))))
    return(long_markerlvl)
}

## Identify engrafted features
get_engraft = function(long_markerlvl, mapfile, cutoff_abs, cutoff_pres,  ...){
    engraft = (long_markerlvl
               # Join the feature count/coverage data with the metadata and
               # filter to the desired groups
               %>% left_join(mapfile, by = c('sample' = 'Study_ID'))
               %>% filter(...)
               # Reshape the table so we can determine what has been engrafted
               %>% select(Marker, Timepoint, abundance, Fig_lab)
               %>% spread(Timepoint, abundance, fill = 0)
               # Set engrafted column to 1 if engraftment criteria are met, 0
               # otherwise
               %>% mutate(engrafted = case_when(WK0 <= cutoff_abs & WK6 >= cutoff_pres ~ 1,
                                                TRUE ~ 0))
               # Select the columns and reshape the data frame to create a wide
               # matrix of 1/0 engraftment values
               %>% select(Marker, Fig_lab,  engrafted)
               %>% pivot_wider(names_from = Fig_lab, values_from = engrafted,
                               values_fill = 0)
               # Make the Marker column the rownames
               %>% column_to_rownames('Marker'))
    
    # Turn it into a matrix and remove all rows that are 0 everywhere. This
    # removal has no effect on the test statistic, it's just tidier
    engraft = as.matrix(engraft)
    engraft = engraft[rowSums(engraft) > 0,]
    return(engraft)
}

## Make a vector of patient IDs/treatment group correspondences
pat_to_vect = function(mapfile, patcol, catcol){
    # Create a two-column data frame of unique patient IDs and the grouping column
    # (probably Treatment or Remission)
    to_perm = mapfile[,c(patcol,catcol)]
    to_perm = unique(to_perm)
    
    # Create a vector of the grouping variable and name it with the patient IDs
    # to make it easy to permute
    perm = to_perm[,catcol]
    names(perm) = to_perm$Fig_lab
    
    # Remove any NAs
    perm = na.omit(perm)
    return(perm)
}

## Count the number of engrafted features in baseline and non-baseline groups
count_engraft = function(engraft, pat_vect, txrm){
    # Is this treatment or remission?
    if (txrm == 'Treatment'){
        bl = 'Placebo'
        tx = 'FMT'
    } else if (txrm == 'Remission'){
        bl = 'NoRes'
        tx = 'Res'
    }
    
    # Subset the engraftment matrix to just the columns identified as FMT or Res
    # (see above) in the patient vector and take the row sums to get the number
    # of times each feature was engrafted in the FMT (Res) group
    engr_tx = engraft[,names(pat_vect[pat_vect == tx])]
    tx_tot = rowSums(engr_tx)
    
    # As above, but with the baseline (Placebo/NoRes) group
    engr_bl = engraft[,names(pat_vect[pat_vect == bl])]
    bl_tot = rowSums(engr_bl)
    
    # create a data frame with the two sets of counts as columns
    engr_ct = cbind(tx_tot, bl_tot)
    return(engr_ct)
}

get_stats = function(cts_mat){
    
    # Get the difference in number of engrafted features (features where
    # engraftment is > 0)
	stat_v = c('fx' = sum(cts_mat[,'tx_tot'] > 0) - sum(cts_mat[,'bl_tot'] > 0),
	# Get the difference in the number of engraftment events           
	           'xfx' = sum(cts_mat[,'tx_tot']) - sum(cts_mat[,'bl_tot']),
	# Get the difference in square-weighted number of engraftment events (weighted
	# by the number of people they are engrafted in)
	           'x2fx' = sum(cts_mat[,'tx_tot']^2) - sum(cts_mat[,'bl_tot']^2))
	return(stat_v)
}

do_permute = function(engr, cts, obs, pv, nperm = 2000, subsample = 0,
                      txrm = NULL){
    
    # Check the inputs. Make sure that Treatment vs. Remission is specified if
    # subsample is set. Set the variables to use lower down.
    if (is.null(txrm) | !(txrm %in% c('Treatment','Remission'))){
        stop(paste('txrm must be set to',
                   '"Treatment" or "Remission"'))
    }
    if (subsample > 0){
        if (txrm == 'Treatment'){
            bl = 'Placebo'
            tx = 'FMT'
        } else if (txrm == 'Remission'){
            bl = 'NoRes'
            tx = 'Res'
        }
    }
    
    # Create an array to put the permuted engraftment counts in
    cts_array = array(dim = c(dim(cts), nperm))
    
    # Add the observed values as the first matrix in the array
    cts_array[,,1] = cts
    
    # Create a matrix to store the permuted test statistic values in
    stat_mat = matrix(nrow = nperm, ncol = 3)
    colnames(stat_mat) = c('fx','xfx','x2fx')
    
    # Add the observed values as the first row in the matrix
    stat_mat[1,] = obs
    
    # Create a matrix to store the permuted patient ID vectors in
    perm_pvs = matrix(nrow = nperm, ncol = length(pv))
    
    # Put the actual patient vector in the first row
    perm_pvs[1,] = pv
    
    # Permute the remaining nperm-1 times
    for (i in 2:nperm){
        
        # Copy the patient vector for permuting
        perm_pv = pv
        
        # We actually permute the names (patient IDs) while keeping the values
        # (group membership) in place. It doesn't matter.
        names(perm_pv) = names(pv)[sample(1:length(pv), length(pv))]
        
        # If subsample >0, that means the groups are not equal sizes. We need to
        # subsample down to the smaller group size and take an average to get
        # unbiased count values
        if (subsample > 0){
            
            # Create an array to hold the subsampled count values
            cts_sub = array(dim = c(dim(cts), 100), 
                            dimnames = list(rownames(cts),
                                            colnames(cts),
                                            NULL))
            # Create a matrix to hold the subsampled test statistics
            stat_sub = matrix(ncol = 3, nrow = 100)
            # subsample 100 times
            for (j in 1:100){
                # Sample from both groups the same depth. This should, in
                # principle, be the depth of one of the groups so sampling will
                # have a null effect. But since I don't know a priori which
                # group that is, it doesn't hurt to do it to both.
                perm_pv_sub = c(sample(perm_pv[perm_pv == tx], subsample),
                                sample(perm_pv[perm_pv == bl], subsample))
                
                # Add the counts from this subsampling to the subsampling array
                cts_sub[,,j] = count_engraft(engr, perm_pv_sub, txrm)
                stat_sub[j,] = get_stats(cts_sub[,,j])
            }
            
            
            # Take the mean of all the subsamplings just to have something
            engr_ct = apply(cts_sub, c(1,2), mean)
            engr_sd = apply(cts_sub, c(1,2), sd)
            
            # Take the mean of the test statistics
            stat_mat[i,] = colMeans(stat_sub)
        } else{
            
            # If we're not subsampling, just get the count matrix and stats once
            engr_ct = count_engraft(engr, perm_pv, txrm)
            stat_mat[i,] = get_stats(engr_ct)
        }
        # Add the count matrix to its array and the test statistics to their
        # matrix
        cts_array[,,i] = engr_ct
        perm_pvs[i,] = perm_pv
    }
    
    # Return a list with the full permutation information and output so it can
    # easily be reused
    ret_l = list('cts_array' = cts_array, 'stat_mat' = stat_mat, 
                'perm_pv' = perm_pvs)
    return(ret_l)
}

get_pvals = function(stat_mat){
	
    # Create a vector for the three p-values
	pvals = c()

	# For each column in the matrix of test-statistics, create a function to
	# identify quantiles and then find the quantile of the first value in the
	# matrix, which is always the observed value
	for (i in colnames(stat_mat)){
		quant = ecdf(stat_mat[,i])
		p = 1 - quant(stat_mat[1,i])
		pvals = c(pvals, p)
	}
    
	# Name the p-values to correspond with their respective test statistics
	names(pvals) = colnames(stat_mat)
	
	# Correct them so they are two-tailed. NEGATIVE VALUES MEANS THE WRONG
	# DIRECTION
	
	pvals = ifelse(pvals < 0.5,
	               2*pvals,
	               -2*(1-pvals))
	return(pvals)
}


plot_fig_13 = function(cts, txrm){
    
    # Check Inputs
    if (!(txrm %in% c('Treatment', 'Remission'))){
        stop('txrm must be one of "Treatment" or "Remission".')
    }
    
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
                                   -Npats))
    f1 = ggplot(plt_df, aes(Npats, Count, colour = .data[[txrm]])) +
        geom_point() +
        geom_line() +
        scale_x_continuous(breaks = 1:max(plt_df$Npats))
    return(f1)
}

plot_fig_2 = function(cts){
    
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
              %>% select(-names)
              %>% filter(Count > 0))
    
    f2 = ggplot(cts_df, aes(x = Total, y = Count, fill = Uniqueness,
                                        colour = Uniqueness)) +
        geom_bar(stat = 'identity') +
        scale_x_continuous(breaks = 1:max(cts_df$Total))
    return(f2)
}
