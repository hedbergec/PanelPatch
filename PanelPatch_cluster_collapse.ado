program define PanelPatch_cluster_collapse
    syntax varlist [if] , marker(varlist) cluster(varlist) ///
        [thresholdratio(real .5)] /// rate of reponders
        [thresholdcases(integer 25)] //number of cases to maintain in each cluster 

    marksample touse, novarlist
    local total_j = wordcount("`varlist'")

    tempvar newk
    quietly : gen `newk' = `cluster' + 1000 if `touse'
    local ok = 0
    local displayrun = 1
    while `ok' == 0 {
        quietly : levelsof `newk' if `touse', local(k_list)
        local total_k = wordcount("`k_list'")
        local total_k2 = `total_k'^2
        local test_counter = 0
        local small_k_list ""
        foreach k in `k_list' {
            quietly : count if `newk' == `k' & `marker' == 1 & `touse'
            local t_`k' = r(N)
            quietly : count if `newk' == `k' & `touse'
            local p_`k' = `t_`k''/r(N)
            if `p_`k'' >= `thresholdratio' & `t_`k'' >= `thresholdcases' {
                local ++test_counter
            }
            else {
                local small_k_list `small_k_list' `k'
            }
        }
        if `test_counter' == `total_k' {
            local ok = 1 //done
            di _newline _newline "Clusters OK"
            quietly : replace `cluster' = `newk' if `touse' //move new values to old variable
        }
        else {
            if `displayrun' == 1 { //if first run, set up dots
                nois _dots 0, title("Adjusting Clusters")
                local run = 1
                local displayrun = 0 //set of 0 to prevent new titles
            }
            quietly : levelsof `newk' if `touse', local(k_list)
            local small_k_list ""
            foreach k in `k_list' {
                quietly : count if `newk' == `k' & `marker' == 1 & `touse'
                local t_`k' = r(N)
                quietly : count if `newk' == `k' & `touse'
                local p_`k' = `t_`k''/r(N)
                if `p_`k'' >= `thresholdratio' & `t_`k'' >= `thresholdcases' {
                    local ++test_counter
                }
                else {
                    local small_k_list `small_k_list' `k'
                }
            }

            if "`small_k_list'" != "" {
                foreach k in `k_list' { //get means of each cluster
                    foreach j in `varlist' { //of each variable 
                        quietly : sum `v' if `touse' & `newk' == `k'
                        local mu_`j'_`k' = r(mean)
                    }
                }

                foreach k in `k_list' { //compute differences across varlist
                    foreach kp in `k_list' {
                        if `k' != `kp' {
                            local d_`k'_`kp' = 0
                            foreach j in `varlist' {
                                local d_`k'_`kp' = `d_`k'_`kp'' + (`mu_`j'_`k''-`mu_`j'_`kp'')^2
                            }
                            local d_`k'_`kp' = sqrt(`d_`k'_`kp'')
                        }
                    }
                }

                local running_dist = .
                local recode_pair ""
                foreach k in `small_k_list' { //find the smallest distance among the smaller ks and all others 
                    foreach kp in `k_list' {
                        if `k' != `kp' {
                            if `d_`k'_`kp'' < `running_dist' & `d_`k'_`kp'' < . {
                                local running_dist = `d_`k'_`kp''
                                local recode_pair "(`k' = `kp')"
                            }
                        }
                    }
                }
                quietly : recode `newk' `recode_pair' //recode just the one
                nois _dots `run' 0
                local ++run
            }
            else {
                local ok = 1 //done
                di _newline _newline "Clusters OK"
                quietly : replace `cluster' = `newk' if `touse' //move new values to old variable
            }

            
            
        }

        
    }
    tempvar stringnewk newnewk
    quietly : tostring `newk', gen(`stringnewk') 
    drop `newk'
    quietly : gen `newnewk' = `stringnewk' if `touse'
    quietly : encode `newnewk', gen(`newk')
    quietly : replace `cluster' = `newk' if `touse'
end