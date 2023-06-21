program define PanelPatch_weight_cluster, rclass //runs cluster_wrapper, checks sizes of donors, and reruns if necessary
    syntax anything (equalok) [if] [in], /// outcomes, row selection
        GENerate(namelist) wave(varlist) i(varlist) j(varlist) response(varlist) ///
        [lassofolds(integer 10)] ///
        [thresholdratio(real .5)] /// rate of reponders
        [thresholdcases(integer 25)] //number of cases to maintain in each cluster 
    *https://acarril.github.io/posts/progess-bar
    
    marksample touse, novarlist

    cluster drop _all 

    gettoken kvars xvars : anything, parse("=") //split the varaible list by equal sign
    gettoken eq xvars : xvars, parse("=") //get the predictors

    quietly : sum `wave' if `touse' //get details about the wave variable and store local macros
    local firstwave = r(min)
    quietly : sum `wave' if `touse' & `wave' != `firstwave'
    local secondwave = r(min)
    quietly : levelsof `wave' if `touse', local(wavelist)
    local numberwaves = wordcount("`wavelist'")
    tokenize "`wavelist'"
    forvalues w = 1/`numberwaves' {
        local w_`w' = ``w''
    }

    

    local total_runs = wordcount("`kvars'") + 1
    nois _dots 0, reps(`total_runs') title("Computing lags of unstable variables and imputations")
    local run = 1
    local lagvars = ""
    foreach v in `response' `kvars' {
        quietly : PanelPatch_mk_compvars `v' if `touse', ///
            id(`i' `j') wave(`wave') 
        local lags = r(lags)
        local lagvars `lagvars' `lags'
        nois _dots `run' 0
        local ++run
    }

    capture drop select_`generate'
    quietly : PanelPatch_running_response `response' if ///
        `touse', i(`i') j(`j') wave(`wave') gen(select_`generate')
    
    label var select_`generate' "Case by Wave Selected for Adjustment"

    di as text _newline _newline  "Cases Requiring Weight Adjustment by Wave"
    tab `wave' `response' if `touse' & select_`generate' == 1

    capture confirm variable `generate'
    if _rc == 0 {
        drop `generate'
    }
    quietly : gen `generate' = .

    local nkvars = wordcount("`kvars'")
    local total_runs = `nkvars' + 1
    forvalues w = 2/`numberwaves'  {
        quietly : sum `response' if `touse' & /// count of folks who are non-responders but have responded so far
            `wave' == `w_`w'' & select_`generate' == 1
        if r(mean) == 1 { //if none-left, just carry over prev wave
            quietly : replace `generate' = 1 if `touse' & /// constant
            `wave' == `w_`w'' & select_`generate' == 1
            di _newline _newline "No more additional attrition at wave `w', cluster is constant"
        }
        else { //if some, get to work
            local toclustervars ""
            di _newline
            nois _dots 0, reps(`total_runs') title("Running linear LASSO Models for response and values for wave `w'. Predictions which do not vary are discarded")
            local run = 1
            foreach v of varlist `response' `kvars' {
                local xset = "i.`j'"
                foreach x of varlist `kvars' `lagvars' `xvars' { //removes related categories of same variable from predictors
                    local xtest = regexr("`x'","_cat_[0-9]*","")
                    if regexm("`x'", "^`xtest'") == 0 {
                        local xset `xset' `x'
                    }
                }
                quietly : lasso linear `v' `xset' if `touse' & ///
                    `wave' == `w_`w'' & select_`generate' == 1 , folds(`lassofolds') 
                tempvar `v'hat`w' 
                quietly : predict ``v'hat`w''  if `touse' & ///
                    `wave' == `w_`w'' & select_`generate' == 1 
                quietly : sum ``v'hat`w'' if `touse' & ///
                    `wave' == `w_`w'' & select_`generate' == 1
                if r(sd) > 0 {
                    quietly : replace ``v'hat`w'' = (``v'hat`w''-r(mean))/r(sd) if `touse' & ///
                        `wave' == `w_`w'' & select_`generate' == 1
                    local toclustervars `toclustervars' ``v'hat`w''
                    nois _dots `run' 0
                }
                else {
                    nois _dots `run' 1
                }
                local ++run
            }
            quietly : replace ``response'hat`w'' = ``response'hat`w'' * sqrt(`nkvars')
            tempvar clustervar`w'
            quietly : count if `touse' & /// 
                `wave' == `w_`w'' & select_`generate' == 1

            local starter_k = floor(r(N)/`thresholdcases')
            di as text _newline "Fitting wave `w' data to `starter_k' clusters"

            quietly : cluster kmeans `toclustervars' if `touse' & /// kmeans on varlist if
                `wave' == `w_`w'' & select_`generate' == 1, /// 
                k(`starter_k') measure(L2) gen(`clustervar`w'') 
            
            *set trace on 
            
            PanelPatch_cluster_collapse `toclustervars' if `touse' & /// 
                `wave' == `w_`w'' & select_`generate' == 1, ///
                marker(`response') cluster(`clustervar`w'')
            
            di _newline "Final clusters for wave `w'"

            tab `clustervar`w'' `response' if `touse' & /// 
                `wave' == `w_`w'' & select_`generate' == 1

            quietly : replace `generate' = `clustervar`w'' if `touse' & /// 
                `wave' == `w_`w'' & select_`generate' == 1
        }
        

    }

end

