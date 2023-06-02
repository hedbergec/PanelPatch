program define PanelPatch_cluster_algorithm, rclass //runs cluster_wrapper, checks sizes of donors, and reruns if necessary
    syntax varlist [if] [in], /// outcomes, row selection
        marker(varlist) GENerate(namelist) /// marker variable for which we need to keep some numbers of tagged values
        [threshold(integer 5)] ///number of cases to maintain in each cluster 
        [starter_k(integer 5)] 
    *https://acarril.github.io/posts/progess-bar
    marksample touse
    capture confirm variable `generate'
    if _rc == 0 {
        drop `generate'
    }
    cluster drop _all 
    tempvar leaf1 leaf2 templeaf2 templeaf1 abandon
    quietly : gen `leaf1' = 1 if `touse'
    quietly : PanelPatch_cluster_checker `leaf1' if `touse', marker(`marker') threshold(`threshold')
    quietly : gen `leaf2' = 0
    quietly : gen `abandon' = 0
    label var `abandon' "Finished cluster"
    local high_k = r(high_k)
    local l1_rok = r(ok)
    local run = 1
    local old_high_k "0"
    else {
        local runlimit = .
    }
    di as result "Clustering, key for runs:" ///
       _newline _col(10) ". is a set of clusters or subclusters with none that have less than `threshold'" ///
       _newline _col(10) "s is a set of clusters or subclusters with too few donors each" ///
       _newline _col(10) "x is a set of clusters or subclusters with too many donors each" 
       
    _dots 0, title("Runs")
    while `l1_rok' >= 0 & "`high_k'" != "." & "`high_k'" != "`old_high_k'" & `run' < `runlimit' {
        local old_high_k `high_k'
        quietly : replace `leaf2' = 0
        foreach l in `high_k' {
            *di "run `run' value `l', starter is `starter_k'"
            quietly : PanelPatch_cluster_wrapper `varlist' if `touse' & `leaf1' == `l', k(`starter_k') gen(`templeaf2') 
            quietly : PanelPatch_cluster_checker `templeaf2' if `touse', marker(`marker') threshold(`threshold')
            local l2_rok = r(ok)
            _dots `run' `l2_rok'
            local ++run
            *di "run `run' value `l' result is `l2_rok'"
            if `l2_rok' == -1 {
                local running_k = `starter_k'
                while `l2_rok' == -1 & `running_k' > 2 {
                    local --running_k
                    *di "run `run' value `l' starter is `running_k'"
                    drop `templeaf2'
                    quietly : PanelPatch_cluster_wrapper `varlist' if `touse' & `leaf1' == `l', k(`running_k') gen(`templeaf2') 
                    quietly : PanelPatch_cluster_checker `templeaf2' if `touse', marker(`marker') threshold(`threshold')
                    local l2_rok = r(ok)
                    *di "run `run' value `l' result is `l2_rok'"
                    _dots `run' `l2_rok'
                    local ++run
                }
                if `l2_rok' != -1 {
                    *di "run `run' value `l' done, transferring values"
                    quietly : replace `leaf2' = `templeaf2' if `touse' & `leaf1' == `l'
                    drop `templeaf2'
                }
                else {
                    *di "run `run' value `l' finished"
                    quietly : replace `leaf2' = `l' if `touse' & `leaf1' == `l'
                    quietly : replace `abandon' = 1 if `touse' & `leaf1' == `l'
                    drop `templeaf2'
                }
            }
            else if `l2_rok' == 1 {
                local running_k = `starter_k' 
                while `l2_rok' == 1 {
                    local ++running_k
                    *di "run `run' value `l' starter is `running_k'"
                    tempvar old_templeaf2 
                    capture drop `old_templeaf2'
                    quietly : gen `old_templeaf2' = `templeaf2'
                    drop `templeaf2'
                    quietly : PanelPatch_cluster_wrapper `varlist' if `touse' & `leaf1' == `l', k(`running_k') gen(`templeaf2') 
                    quietly : PanelPatch_cluster_checker `templeaf2' if `touse', marker(`marker') threshold(`threshold')
                    local l2_rok = r(ok)
                    *di "run `run' value `l' result is `l2_rok'"
                    _dots `run' `l2_rok'
                    local ++run
                }
                if `l2_rok' == -1 {
                    *di "run `run' value `l' finished"
                    quietly : replace `leaf2' = `old_templeaf2' if `touse' & `leaf1' == `l'
                    drop `templeaf2' `old_templeaf2' 
                }
                else if `l2_rok' == 0 {
                    *di "run `run' value `l' done, transferring values"
                    quietly : replace `leaf2' = `templeaf2' if `touse' & `leaf1' == `l'
                    drop `templeaf2'
                }
            }
            else {
                *di "run `run' value `l' done, transferring values"
                quietly : replace `leaf2' = `templeaf2' if `touse' & `leaf1' == `l'
                drop `templeaf2'
            }
        }
        quietly : recode `leaf2' (.=0) if `touse' 
        quietly : egen `templeaf1' = group(`leaf1' `leaf2') if `touse' 
        quietly : replace `leaf1' = `templeaf1' if `touse' 
        *table `leaf1'  `abandon' `marker' if `touse' 
        quietly : replace `leaf2' = 0
        drop `templeaf1' 
        quietly : PanelPatch_cluster_checker `leaf1' if `touse' & `abandon' == 0, marker(`marker') threshold(`threshold')
        local high_k = r(high_k)
        local l1_rok = r(ok)
        _dots `run' `l1_rok'
        local ++run
        *di "run `run' done, result is `l1_rok'" 
    }
    di _newline "Solution found"
    quietly : gen `generate' = `leaf1'
    table `generate' `marker'


end

