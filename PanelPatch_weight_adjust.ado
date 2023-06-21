
program define PanelPatch_weight_adjust, rclass
    version 17 //written for stata 17
    syntax varlist [if], /// varlist is the weight
        j(varlist) i(varlist) Wave(varlist) GENerate(namelist) /// gen is a stem for the factor and cluster
        saving(string) ///
        [interimfile(string)] ///
        [Xvars(varlist)] ///
        [SNumericvars(varlist)] ///
        [SUnorderedvars(varlist)] ///
        [SOrderedvars(varlist)]  ///
        [VNumericvars(varlist)] ///
        [VUnorderedvars(varlist)] ///
        [VOrderedvars(varlist)] ///
        [lassofolds(real 10)] ///
        [Burnin(integer 10)] 
    
    ****This program assumes the data in memory have been imputed using stata's MI IMPUTE, 
    ****but some missing values remain This program imputes the remaining items and finds clusters 
    ****to create weight adjustment factors

    marksample touse, novarlist

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

**?? JW 06/01/2023: added clear option
    quietly : mi export ice, clear //make long mi export

    **bring in interim imputations

    if "" != "`interimfile'" {
        quietly : merge 1:1 `i' `j' `wave' _mj using "`interimfile'", nogen update 
    }

    quietly : drop if _mj == 0
    *gather lists of values to use, if imputations exist, use them

    local total_runs = wordcount("`sunorderedvars' `vunorderedvars' `snumericvars' `sorderedvars' `vnumericvars' `vorderedvars'")
    display _newline _newline
    nois _dots 0, reps(`total_runs') title("Preparing `total_runs' stable and unstable variables for item imputation")
    local run = 1

    foreach v in `snumericvars' `sorderedvars'  `vnumericvars' `vorderedvars'  {
        quietly : PanelPatch_mk_compvars `v' if `touse', ///
                id(`i' `j' _mj) wave(`wave')   
        local intrps = r(intrps) 
        if "`interps'" != "." {
            local yintrps `yintrps' `intrps'
        }
        local y_use_`v' `v'
        local yvars `yvars' `v'
        nois _dots `run' 0
        local ++run
    }
 

    foreach v in `sunorderedvars'  `vunorderedvars'  {
        quietly : PanelPatch_mk_compvars `v' if `touse', ///
                id(`i' `j' _mj) wave(`wave') categorical
        local intrps = r(intrps)
        local cats = r(cats)
        if "`cats'" == "." {
            local y_use_`v' `v'
            local yvars `yvars' `v'
        }
        else {
            local y_use_`v' `cats'
            local yvars `yvars' `cats'
        }
        local yintrps `yintrps' `intrps'
        nois _dots `run' 0
        local ++run
        
    }


    *collapse over the imputations 

    quietly : collapse `yvars' `yintrps' `varlist' if `touse', by(`i' `j' `wave' `touse') fast 

    *item impute stable 
    local nr_s_vars = ""
    local chainedspecs = ""
    local sregisterlist = ""
    
    foreach v in `snumericvars' `sorderedvars' `sunorderedvars' {
        foreach y in `y_use_`v'' {
            local omitlist = subinstr("`y_use_`v''","`y'","",1)
            if "`omitlist'" != "" {
                local omitcall = "omit(`omitlist')"
            }
            else {
                local omitcall = ""
            }
            local chainedspecs `chainedspecs' (pmm, knn(5) `omitcall') `y'
            local sregisterlist `sregisterlist' `y'
            quietly : gen nr_`y' = `y' >= .
            local nr_s_vars `nr_s_vars' nr_`y'
        }
        
    }

    quietly : mi set wide
    quietly : mi register imputed `sregisterlist'
    di _newline "Imputing Remaining Stable Values In First Wave"
    mi impute chained ///
        `chainedspecs' = i.`j' `xvars' if ///
        `touse' & `wave' == `firstwave', add(1) burnin(`burnin')
    
    if r(M) > 0 {
        quietly : mi update

        foreach v in `sregisterlist' {
            tempvar the_constant
            quietly : bysort `i' `j' : egen `the_constant' = mean(_1_`v')
            quietly : replace `v' = `the_constant' if `v' >= .
            drop `the_constant'
        }

        quietly : mi unset 

        capture drop *_1_
    }

    local stablex `sregisterlist'

    *item impute unstable 
    local nr_v_vars = ""
    local chainedspecs = ""
    local vregisterlist = ""
    foreach v in `vnumericvars' `vorderedvars' `vunorderedvars'   {
        foreach y in `y_use_`v'' {
            local omitlist = subinstr("`y_use_`v''","`y'","",1)
            if "`omitlist'" != "" {
                local omitcall = "omit(`omitlist')"
            }
            else {
                local omitcall = ""
            }
            local chainedspecs `chainedspecs' (pmm, knn(5) `omitcall') `y'
            local vregisterlist `vregisterlist' `y'
            quietly : gen nr_`y' = `y' >= .
            local nr_v_vars `nr_v_vars' nr_`y'
        }
        
        
    } 

    quietly : mi set wide
    quietly : mi register imputed `vregisterlist'
    di _newline "Imputing Remaining Unstable Values In All Waves"
    mi impute chained ///
        `chainedspecs' = i.`j' c.`wave'##c.`wave' `yintrps' `stablex' `xvars' if ///
        `touse', add(1) burnin(`burnin')
    
    quietly : mi update

    foreach v in `vregisterlist' {
        quietly : replace `v' = _1_`v' if `v' >= .
    }

    quietly : mi unset 

    capture drop *_1_

    tempvar nr_s nr_v keepweight

    quietly : egen `nr_s' = rowmin(`nr_s_vars') if `touse'
    quietly : egen `nr_v' = rowmin(`nr_v_vars') if `touse'

    *quietly : sum `nr_s'
    *assert r(mean) == 0 //check all stable are imputed

    capture drop ris_`generate'

    quietly : recode `nr_v' (0 = 1) (1 = 0), gen(ris_`generate') //keepweight is the marker for clustering
    label def keepweight 1 "Responder or Imputed" 0 "Nonresponder without Imputations"
    label var ris_`generate' "Responder or Imputed Status"
    label val ris_`generate' keepweight

    
    quietly : keep if `touse'
    keep `i' `j' ris_`generate'  `wave' `varlist' `vregisterlist' `yintrps' `stablex' `xvars' `touse'

    quietly : save "`saving'", replace

    PanelPatch_weight_cluster ///
        `vregisterlist'  = `xvars' `stablex' if `touse', ///
        response(ris_`generate') gen(clst_`generate') ///
        wave(`wave') j(`j') i(`i')

    foreach wgtvar in `varlist' {
        capture drop `wgtvar'_`generate'
        capture drop f_`wgtvar'_`generate'
        quietly : gen f_`wgtvar'_`generate' = .
        quietly : gen `wgtvar'_`generate' = `wgtvar' if `wave' == `w_1' & `touse' & select_clst_`generate' == 1
        forvalues w = 2/`numberwaves' {
            local lastwave = `w' - 1
            tempvar carryforward constantw totallyingw numerw denomw  
            quietly : gen `carryforward' = `wgtvar'_`generate' if ///
                `wave' == `w_`lastwave'' & `touse' 
            quietly : bysort `i' `j' : egen `constantw' = mean(`carryforward') if `touse'
            quietly : gen `totallyingw' = `constantw' if `touse' ///
                & `wave' == `w_`w'' & select_clst_`generate' == 1
            quietly : bysort clst_`generate' : egen `numerw' = total(`totallyingw') if ///
                `wave' == `w_`w'' & select_clst_`generate' == 1 & `touse' 
            quietly : replace `totallyingw' = 0 if `touse' ///
                & `wave' == `w_`w'' & select_clst_`generate' == 1 & ris_`generate' == 0
            quietly : bysort clst_`generate' : egen `denomw' = total(`totallyingw') if ///
                `wave' == `w_`w'' & select_clst_`generate' == 1 & `touse' 
            
            quietly : replace f_`wgtvar'_`generate' = `numerw'/`denomw' if ///
                `wave' == `w_`w'' & select_clst_`generate' == 1 & `touse' 

            quietly : replace `wgtvar'_`generate' = `totallyingw'*f_`wgtvar'_`generate' if ///
                `wave' == `w_`w'' & select_clst_`generate' == 1 & `touse'
            capture drop `carryforward' 
            capture drop `constantw' 
            capture drop `totallyingw' 
            capture drop `numerw' 
            capture drop `denomw' 
        }

        quietly : replace `varlist'_`generate' = 0 if `varlist'_`generate' == .
    }

    keep `i' `j' `wave' *_`generate' 

    quietly : describe *_`generate' 

    local newvars = r(varlist)

    return local newvars = "`newvars'"

    quietly : save "`saving'", replace


end


