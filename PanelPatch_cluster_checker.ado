program define PanelPatch_cluster_checker, rclass // program to check if clusters from k-means wrapper (below) have more or less cases
    syntax varlist [if] [in], ///varlist is the cluster var, totally marker
        marker(varlist) /// marker variable for which we need to keep some numbers of tagged values
        threshold(numlist) [thresholdfactor(real 3.0)] // number of cases to maintain in each cluster 
    marksample touse
    tempname N m r 
    levelsof `varlist' if `touse', local(k_list) //get values of cluster var
    local number_k = wordcount("`k_list'") //number k
    **count number donors r in each cluster
    local too_low = 0 //test starts at 0
    local too_high = 0 //test starts at 0
    
    *return values of too-high groups in local macro
    local high_k = "" //start list of high ones
    local high_k_values = ""
    local n_high_k = 0
    foreach j in `k_list' { //go through each cluster
        count if `touse' & `varlist' == `j' & `marker' == 1
        local r = r(N)
        if  `r' <= `threshold' {
            local ++too_low //mark test if is less than or equal to threshold
        }
        if  `r' > `threshold'*`thresholdfactor' {
            local ++too_high //mark test if r is 3 times threshold
            local high_k `high_k' `j'
            local n_high_k = `n_high_k' + `r'
            local high_k_values `high_k_values' `r'
        } 
    }
    *return value of 0 if there are between r and 3*r cases (inclusive)
    *return value of -1 if there are r or less in any one cluster
    *return value of 1 if there are more than 3*r in all clusters
    
    if `too_low' > 0 {
        return scalar ok = -1 //one is too low
    }
    else {
        if `too_high' == `number_k' { // twice more in each 
            return scalar ok = 1 //all too high
        }
        else {
            return scalar ok = 0 //just right
        }
    }
    local high_k_values : list sort high_k_values
    return local high_k "`high_k'"
    return local high_k_values "`high_k_values'"
    return scalar n_high_k = `n_high_k'
end

