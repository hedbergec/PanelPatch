//conductor program to perform all steps
// do  5 interim

program define PanelPatch , rclass
    syntax varlist [if] [in] , ///
        add(integer) j(varlist) i(varlist) Wave(varlist) ///
        WAVEResponseflag(varlist) ///
        [weightvar(varlist)] /// weightvar or set of vars
        [Burnin(integer 10)] ///
        [lassofolds(integer 10)] ///
        [Xvars(varlist)] ///
        [SNumericvars(varlist)] ///
        [SUnorderedvars(varlist)] ///
        [SOrderedvars(varlist)]  ///
        [VNumericvars(varlist)] ///
        [VUnorderedvars(varlist)] ///
        [VOrderedvars(varlist)] ///
        [diagnosticdata(string)] ///
        [minwave(integer 3)] [RUNDIAGnostic] [useold]

    
    /*
    Logic of this program:
        This program makes multiple calls to Stata's native MI package and 
        performs serveral other data modeling and engineering tasks. Except 
        for the last call, which uses our code, every MI run is converted to 
        ICE MI formats and MI settings are cleared. Data produced at each
        step is stored in temporary data sets, and the user can specify that
        each is saved in a real location for inspection. Lots of "preserve" 
        and "restore" commands are used, so if the user "breaks" the routine
        the data they had when executing the command is restored. 

        useold is not documented and is used for testing purposes. it assumes
        the files have all the variables imputed. this program does not check 
        for that. 
    */


    marksample touse, novarlist //mark which rows to use, ignoring missing status

    quietly : gen _touse = `touse'

    ****check user input for errors****

    local ivars "`varlist'"
    local vars_sorted : list sort ivars
    local option_vars "`snumericvars' `sorderedvars' `sunorderedvars' `vnumericvars' `vorderedvars' `vunorderedvars'"
    local option_vars_uniq : list uniq option_vars
    local option_vars_sorted : list sort option_vars_uniq
    local options_vars_sizeof : list sizeof options_vars
    local options_vars_unique_sizeof : list sizeof options_vars_uniq

    *check that variables are not repeated between variable types in options
    if `options_vars_sizeof' != `options_vars_unique_sizeof' {
        display as error "more variables listed in options (sorderedvars, " ///
            "sunorderedvars, vorderedvars, and vunorderedvars) than unique variables across them"
        error 100 
    }
    *check that ivars and uniq options var are the same, it not, find the problem
    if "`vars_sorted'" != "`option_vars_sorted'" { //these two strings should be the same, if not
        foreach v in `vars_sorted' { //loop over the MI delared variables
            local test_result = regexm("`option_vars'", "`v'") //do the options vars mention the MI declared var?
            if `test_result' == 0 { //if not, error out and describe issue
                display as error "`v' is listed to be imputed, " ///
                "but `v' is not in one of the options"
                error 100
            }
        }
        foreach v in `option_vars_sorted' { //loop over option variable
            local test_result = regexm("`ivars'", "`v'") //do the declared vars mention the option var?
            if `test_result' == 0 { //if not, error out and describe issue
                display as error "`v' is listed in the options, " ///
                "but not in variables to be imputed: `ivars'"
                error 100
            }
        }
    }

    if "`snumericvars' `sorderedvars' `sunorderedvars'" == "" {
        display as error "at least one stable variable noted in one or more options snumericvars, sorderedvars, or sunorderedvars, required."
        error 100
    }

    if "`vnumericvars' `vorderedvars' `vunorderedvars'" == "" {
        display as error "at least one unstable variable noted in one or more options vnumericvars, vorderedvars, or vunorderedvars, required."
        error 100
    }

    **** make sure i and j uniquely id cases ****

    quietly : duplicates report `i' `j' `wave' if `touse'
    capture assert r(unique_value) == r(N)
    if _rc != 0 {
        di as error "`i', `j', and `wave' `if' do not uniquely identify observations"
        error 9
    }

    **** sort out waves, get wave values ****
    quietly : sum `wave' if `touse'
    local firstwave = r(min)
    local lastwave = r(max)
    quietly : levelsof `wave' if `touse', local(wavelist)
    local nwaves = wordcount("`wavelist'")
    tokenize "`wavelist'"
    forvalues w = 1/`nwaves' {
        local w_`w'_value = ``w''
    }

    **** clear any mi settings, display note is there were any

    capture mi unset 
    if _rc == 0 {
        display as result "NOTE: Data were MI SET prior to this command, and these settings have been cleared"   
    }

    **** mark non-missing waves of variables ****
    tempvar number_missing
    quietly : egen `number_missing' = rowmiss(`ivars') if `touse'
    label var `number_missing' "PanelPatch number ivars missing"
    *compute "by hand" complete wave status
    quietly : describe `ivars'
    local number_ivars = r(k)
    tempvar nonmissing_wave 
    quietly : gen `nonmissing_wave' = `number_missing' < `number_ivars' if `touse'
    label var `nonmissing_wave' "PanelPatch non-missing wave"

    quietly : replace `nonmissing_wave' = `nonmissing_wave' * `waveresponseflag' if `touse'

    **** save data for which interim imputations will be stored

    if "`useold'" != "" & "`diagnosticdata'" != "" {
        preserve
        foreach f in _interim _predictions _donorcluster _donoriddata _weightadjust {
            capture confirm file "`diagnosticdata'`f'.dta"
            if _rc == 0 {
                di "`diagnosticdata'`f'.dta found"
                quietly : use "`diagnosticdata'`f'.dta", clear
                local `f'_file "`diagnosticdata'`f'.dta"
                local `f'_time = c(filedate)
                local `f'_save_time = clock("``f'_time'","DMYhm")
            }
            else {
                di "`diagnosticdata'`f'.dta not found"
                local `f'_save_time = .
            }
        }
        restore
    }
    else {
        foreach f in _interim _predictions _donorcluster _donoriddata _weightadjust {
            local `f'_save_time = .
        }
    }

    if "`_interim_file'" == "" { //start interim
        preserve
        quietly : keep if `touse' 

        *create comp varaibles 

        local total_runs = wordcount("`vnumericvars' `vorderedvars' `vunorderedvars'")
        display _newline _newline
        nois _dots 0, reps(`total_runs') title("Computing companion variables for `total_runs' unstable variables before imputations")
        local run = 1

        local lag_lead_vars ""
        local lead_vars ""
        local lag_vars ""

        foreach v in `vnumericvars' {
            quietly : PanelPatch_mk_compvars `v' if `touse', id(`i' `j') wave(`wave') 
            local `v'_lags = r(lags)
            local `v'_leads = r(leads)
            local lag_lead_vars `lag_lead_vars' ``v'_lags' ``v'_leads'
            local lag_vars `lag_vars' ``v'_lags' 
            local lead_vars `lead_vars' ``v'_leads'
            nois _dots `run' 0
            local ++run
        }

        foreach v in `vorderedvars' `vunorderedvars' {
            quietly : PanelPatch_mk_compvars `v' if `touse', id(`i' `j') wave(`wave') categorical
            local `v'_lags = r(lags)
            local `v'_leads = r(leads)
            local lag_lead_vars `lag_lead_vars' ``v'_lags' ``v'_leads'
            local lag_vars `lag_vars' ``v'_lags' 
            local lead_vars `lead_vars' ``v'_leads'
            nois _dots `run' 0
            local ++run
        }
    
        tempfile interim_imputations
        quietly : save `interim_imputations' , replace
        restore //bring back user state

        **** stable variable interim imputations ****

        preserve //preserve user state 

        quietly : keep if `touse' & `wave' == `w_1_value' & `nonmissing_wave' == 1 //keep only applicable rows of first wave

        display as text _newline _newline "Imputing first wave values of stable variables"

        mi set flong
        mi register imputed `snumericvars' `sorderedvars' `sunorderedvars' 

        local chainedspecs ""

        if "`snumericvars'" != "" {
            local chainedspecs "`chainedspecs' (pmm, knn(5)) `snumericvars'"
        } 
        if "`sorderedvars'" != "" {
            local chainedspecs "`chainedspecs' (ologit, augment ) `sorderedvars'" //augment is for perfect predictons
        }
        if "`sunorderedvars'" != "" {
            local chainedspecs "`chainedspecs' (mlogit, augment ) `sunorderedvars'" //augment is for perfect predictons
        }

        mi impute chained `chainedspecs' = `xvars', add(`add') burnin(`burnin')

        quietly : mi export ice, clear
        tempfile stable_imputations
        drop `vnumericvars' `vorderedvars' `vunorderedvars' `wave' *_0* _mi
        quietly : save `stable_imputations', replace

        quietly : use `stable_imputations', clear
        quietly : keep if _mj == 1 
        drop _m* 

        tempfile stable_m1_values
        quietly : save `stable_m1_values'

        quietly : use `interim_imputations', clear 
        quietly : merge m:1 `i' `j' using `stable_m1_values', update nogen 
        quietly : save `interim_imputations', replace 

        restore //bring back user state

        **** unstable variable interim imputations ****

        preserve //preserve user state

        quietly : merge m:1 `i' `j' using `stable_m1_values', update nogen 

        *put the i. in front of stable categorical
        local stable_predictors `snumericvars' 
        foreach v in `sorderedvars' `sunorderedvars' {
            local stable_predictors `stable_predictors' i.`v'
        }

        *** specs differ by first, last, or other wave

        local chainedspecs_first ""
        foreach v in `vnumericvars' {
            local chainedspecs_first "`chainedspecs_first' (pmm, knn(5) include(`lead_vars') ) `v'"
        } 
        foreach v in `vorderedvars' {
            local chainedspecs_first "`chainedspecs_first' (ologit, augment include(`lead_vars') ) `v'" //augment is for perfect predictons
        }
        foreach v in `vunorderedvars' {
            local chainedspecs_first "`chainedspecs_first' (mlogit, augment include(``v'_leads') ) `v'" //augment is for perfect predictons
        }

        local chainedspecs_between ""
        foreach v in `vnumericvars' {
            local chainedspecs_between "`chainedspecs_between' (pmm, knn(5) include(`lag_lead_vars')) `v'"
        } 
        foreach v in `vorderedvars' {
            local chainedspecs_between "`chainedspecs_between' (ologit, augment include(`lag_lead_vars')) `v'" //augment is for perfect predictons
        }
        foreach v in `vunorderedvars' {
            local chainedspecs_between "`chainedspecs_between' (mlogit, augment include(``v'_lags' ``v'_leads')) `v'" //augment is for perfect predictons
        }

        local chainedspecs_last ""
        foreach v in `vnumericvars' {
            local chainedspecs_last "`chainedspecs_last' (pmm, knn(5) include(`lag_vars')) `v'"
        } 
        foreach v in `vorderedvars' {
            local chainedspecs_last "`chainedspecs_last' (ologit, augment include(`lag_vars')) `v'" //augment is for perfect predictons
        }
        foreach v in `vunorderedvars' {
            local chainedspecs_last "`chainedspecs_last' (mlogit, augment include(``v'_lags')) `v'" //augment is for perfect predictons
        }

        **** First wave ****

        display as text _newline _newline "Imputing unstable variables in first wave"

        quietly : use `interim_imputations' if ///
            `wave' == `firstwave' & `nonmissing_wave' == 1 ///
            , clear 

        mi set flong
        mi register imputed `vnumericvars' `vorderedvars' `vunorderedvars' 

        mi impute chained `chainedspecs_first' = `xvars' `stable_predictors', ///
            add(`add') burnin(`burnin')

        quietly : mi export ice, clear
        tempfile unstable_imputations
        quietly : save `unstable_imputations', replace

        **** Between waves ****

        display as text _newline _newline "Imputing unstable variables in waves between first and last"

        quietly : use `interim_imputations' if ///
            `wave' != `lastwave' & `wave' != `firstwave' & `nonmissing_wave' == 1 ///
            , clear 

        mi set flong
        mi register imputed `vnumericvars' `vorderedvars' `vunorderedvars' 

        mi impute chained `chainedspecs_between' = `xvars' `stable_predictors' ///
            c.`wave'##c.`wave' , ///
            add(`add') burnin(`burnin')

        quietly : mi export ice, clear
        quietly : append using `unstable_imputations'
        quietly : save `unstable_imputations', replace

        **** Last wave ****

        display as text _newline _newline "Imputing unstable variables in last wave"

        quietly : use `interim_imputations' if ///
            `wave' == `lastwave' & `nonmissing_wave' == 1 ///
            , clear 

        mi set flong
        mi register imputed `vnumericvars' `vorderedvars' `vunorderedvars' 

        mi impute chained `chainedspecs_last' = `xvars' `stable_predictors', ///
            add(`add') burnin(`burnin')

        quietly : mi export ice, clear
        quietly : append using `unstable_imputations'
        drop `snumericvars' `sorderedvars' `sunorderedvars'  *_0* _mi
        quietly : save `unstable_imputations', replace

        quietly : use `unstable_imputations', clear
        quietly : keep if _mj == 1 
        drop _m* 

        tempfile unstable_m1_values
        quietly : save `unstable_m1_values'

        **** fill in values in each wave of interim imputation file

        quietly : use `interim_imputations', clear 
        quietly : merge 1:1 `i' `j' `wave' using `unstable_m1_values', update nogen 
        
        quietly : save `interim_imputations', replace 

        restore //restore user state

        *** add all imputations up to this point ***

        preserve //perserve user state

        local expand = `add' + 1

        quietly : expand `expand'

        quietly : bysort `i' `j' `wave' : gen _mj = _n - 1

        quietly : merge m:1 `i' `j' _mj using `stable_imputations', update nogen 
        quietly : merge 1:1 `i' `j' `wave' _mj using `unstable_imputations', update nogen 
        drop *_lag *_lead *_intrp
        quietly : save `interim_imputations', replace
        
        if "`diagnosticdata'" != ""  {
            save "`diagnosticdata'_interim.dta", replace 
            
        }

        restore //restore user state
    } //end interim
    else {
        preserve 
        use "`diagnosticdata'_interim.dta", clear 
        tempfile interim_imputations
        quietly : save `interim_imputations', replace
        
        restore 
    }
    ***predict missing waves
    if  "`_predictions_file'" == "" | ///
        `_interim_save_time' >= `_predictions_save_time' { //start predictions
        
        preserve

        quietly : use `interim_imputations' if _mj != 0, clear

        display _newline _newline "Computing missing-wave predictions"

        *make _intrp comp vars

        local total_runs = wordcount("`vnumericvars' `vorderedvars' `vunorderedvars'")

        nois _dots 0, reps(`total_runs') title("Computing companion variables for `total_runs' unstable variables after imputations")
        local run = 1

**?? JW 06/02/2023: fixed hardcoded i & j variables
**?? JW 06/02/2023: PanelPatch breaks if user doesn't include numeric repeating variable type. "foreach v of varlist   {"
        foreach v in `vorderedvars' `vunorderedvars'  {
            quietly : PanelPatch_mk_compvars `v' if `touse', id(`i' `j' _mj) wave(`wave') categorical
            nois _dots `run' 0
            local ++run
        }

**?? JW 06/02/2023: fixed hardcoded i & j variables
        foreach v in `vnumericvars' {
            quietly : PanelPatch_mk_compvars `v' if `touse', id(`i' `j' _mj) wave(`wave') 
            nois _dots `run' 0
            local ++run
        }
        di _newline _newline
        nois _dots 0, reps(`total_runs') title("Running LASSO variable selection and prediction models for `total_runs' unstable variables")
        local run = 1

        local vars_to_trend ""
        
        foreach v in `vnumericvars' {
            quietly : PanelPatch_predict_wrapper `v' i.(`sorderedvars' `sunorderedvars') ///
                c.(`snumericvars') *_intrp if `touse', ///
                xalways(i.`j') model(regress) wave(`wave')
            local out_vars = r(out_vars)
            local vars_to_trend `vars_to_trend' `out_vars'
            nois _dots `run' 0
            local ++run
        }

        foreach v in `vorderedvars'  {
            quietly : PanelPatch_predict_wrapper `v' i.(`sorderedvars' `sunorderedvars') ///
                c.(`snumericvars') *_intrp if `touse', ///
                xalways(i.`j') model(ologit) wave(`wave')
            local out_vars = r(out_vars)
            local vars_to_trend `vars_to_trend' `out_vars'
            nois _dots `run' 0
            local ++run
        }

        foreach v in `vunorderedvars'  {
            quietly : PanelPatch_predict_wrapper `v' i.(`sorderedvars' `sunorderedvars') ///
                c.(`snumericvars') *_intrp if `touse', ///
                xalways(i.`j') model(mlogit) wave(`wave')
            local out_vars = r(out_vars)
            local vars_to_trend `vars_to_trend' `out_vars'
            nois _dots `run' 0
            local ++run
        }
        char _dta[vars_to_trend] `vars_to_trend'

        quietly : keep if `touse'

        tempfile predictions
        quietly : save `predictions', replace

        if "`diagnosticdata'" != ""  {
            save "`diagnosticdata'_predictions.dta", replace 
            
        }
        
        restore 

    } //end predictions
    else {
        preserve
        use "`diagnosticdata'_predictions.dta", clear
        local vars_to_trend : char _dta[vars_to_trend]
        tempfile predictions
        quietly : save `predictions', replace
        
        restore
    }

    **trending and clustering

    if  "`_donorcluster_file'" == "" | ///
        `_interim_save_time' >= `_predictions_save_time' | ///
        `_predictions_save_time' >= `_donorcluster_save_time' | ///
        `_interim_save_time' >= `_donorcluster_save_time' { //start trending and clustering

        preserve

        use `predictions', clear 
        local vars_to_trend : char _dta[vars_to_trend]

        quietly : egen _PanelPatchallvarscheck = rowmiss(`ivars') 
        quietly : recode _PanelPatchallvarscheck (0 = 1) (1/max = 0)
        quietly : bysort `i' `j' `wave' : egen _minPPcheck = min(_PanelPatchallvarscheck)
        quietly : replace _PanelPatchallvarscheck = _minPPcheck
        quietly : replace _PanelPatchallvarscheck = . if _mj != 1
        drop _minPPcheck

        local trend_collapse_list ""
        local varreps = wordcount("`vars_to_trend'") 
        di _newline
        nois _dots 0, reps(`varreps') title("Computing OLS trends of `varreps' unstable values and categories for each `j' and `i'")
        local run = 1
        foreach v in `vars_to_trend' { 
            quietly : PanelPatch_poly_trend_by `v' , ///
                time(`wave') j(`i' `j') polymax(2) ///
                gen(trend_`v'_0 trend_`v'_1 trend_`v'_2)
            local return_list = r(return_list) //get make vars
            local trend_collapse_list `trend_collapse_list' `return_list' //add them to collapse list
            nois _dots `run' 0
            local ++run
        }

        foreach v in `sunorderedvars' {
            quietly : tab `v', gen(`v'_cat_)
            drop `v'_cat_1 `v'
            quietly : describe `v'_cat_*, varlist
            local dummies_to_add = r(varlist)
            local trend_collapse_list "`trend_collapse_list' `dummies_to_add'"
        }

        foreach v in `sorderedvars' {
            quietly : bysort `i' `j' : sum `v' if _n == 1, meanonly
            quietly : replace `v' = `v'-r(mean)
            quietly : gen `v'_poly_2 = `v'^2
            local trend_collapse_list "`trend_collapse_list' `v' `v'_poly_2"
        }

        foreach v in  `snumericvars' {
            local trend_collapse_list "`trend_collapse_list' `v'"
        }

        quietly : replace `waveresponseflag' = . if _mj != 1
        
        quietly : collapse (mean) `trend_collapse_list' (sum) `waveresponseflag' _PanelPatchallvarscheck, by(`i' `j') fast

        foreach v in `trend_collapse_list' {
            quietly : sum `v', meanonly
            quietly : replace `v' = `v' - r(mean)
        }

        quietly : gen _donormarker = `waveresponseflag' == `nwaves' & _PanelPatchallvarscheck == `nwaves' 
        quietly : gen _imputemarker = `waveresponseflag' >= `minwave' & `waveresponseflag' < `nwaves' & _donormarker == 0
        tab _donormarker _imputemarker

        drop _PanelPatchallvarscheck
        PanelPatch_cluster_algorithm `trend_collapse_list' if _donormarker == 1 | _imputemarker == 1, /// outcomes, row selection
            marker(_donormarker ) /// marker variable for which we need to keep some numbers of tagged values
            gen(_donorcluster) 

        drop `trend_collapse_list' 

        tempfile cluster_data
        quietly : save `cluster_data', replace

        if "`diagnosticdata'" != ""  {
            save "`diagnosticdata'_donorcluster.dta", replace 
        }

        restore 
    
    } //end clustering and trending
    else {
        preserve 
        use "`diagnosticdata'_donorcluster.dta", replace 
        tempfile cluster_data
        quietly : save `cluster_data', replace
        restore 
    }

    *** making ids

    if  "`_donoriddata_file'" == "" | ///
        `_interim_save_time' >= `_predictions_save_time' | ///
        `_predictions_save_time' >= `_donorcluster_save_time' | ///
        `_interim_save_time' >= `_donorcluster_save_time' | ///
        `_donorcluster_save_time' >= `_donoriddata_save_time' | ///
        `_predictions_save_time' >= `_donoriddata_save_time' | ///
        `_interim_save_time' >= `_donoriddata_save_time' { //start making ids

        preserve 

        di _newline "Making assignments"

        use `cluster_data', clear 

        quietly : gen _donorid_master = .
        forvalues m = 1/`add' {
            quietly : gen _donorid`m' = . 
        }
        tempvar tempdonorid
        quietly : levelsof _donorcluster, local(cluster_list)
        foreach k in `cluster_list' {
            quietly : egen `tempdonorid' = group(`i' `j') if _donorcluster == `k' & _donormarker == 1 & _imputemarker != 1 
            quietly : sum `tempdonorid' if _donorcluster == `k' & _donormarker == 1 & _imputemarker != 1 
            local min = r(min)
            local max = r(max)
            quietly : replace _donorid_master = `tempdonorid' if ///
                    _donorcluster == `k' & _donormarker == 1 & _imputemarker != 1 
            forvalues m = 1/`add' {
                quietly : replace _donorid`m' = `tempdonorid' if ///
                    _donorcluster == `k' & _donormarker == 1 & _imputemarker != 1
                quietly : replace _donorid`m' = runiformint(`min',`max') if ///
                    _donorcluster == `k' & _donormarker != 1 & _imputemarker == 1
            }
            capture drop `tempdonorid'
        } 
        
        quietly : reshape long _donorid, i(`i' `j' _donorcluster _donormarker) j(_mj)
         
        tempfile donorid_data
        quietly : save `donorid_data', replace

        quietly : use `interim_imputations', clear //this file is OG state, expanded

        quietly : merge m:1 `i' `j' _mj using `donorid_data', nogen

        quietly : duplicates report `i' `j' `wave' _mj _donorid _donorcluster if _donormarker == 1
        assert r(unique_value) == r(N)

        quietly : save `donorid_data', replace

        if "`diagnosticdata'" != ""  {
            save "`diagnosticdata'_donoriddata.dta", replace 
            local diagdata `"`diagdata'"' `"`diagnosticdata'_donoriddata.dta"'

        }

        restore 
    
    
    } //end making ids
    else {
        preserve
        use "`diagnosticdata'_donoriddata.dta", clear
        tempfile donorid_data
        quietly : save `donorid_data', replace
        local diagdata `"`diagdata'"' `"`diagnosticdata'_donoriddata.dta"'
        restore 
    }

    di _newline "Starting Value Transfer"
    *typical stata MI set up

    quietly : bysort `i' `j' : egen _total_complete_waves = total(`nonmissing_wave')

    mi set wide
    mi register imputed `snumericvars' `sorderedvars' `sunorderedvars' ///
        `vnumericvars' `vorderedvars' `vunorderedvars'
    
    capture drop _PanelPatchImputed

    mi impute PanelPatch  ///
        `snumericvars' `sorderedvars' `sunorderedvars' ///
        `vnumericvars' `vorderedvars' `vunorderedvars' ///
        if `touse' == 1 & _total_complete_waves >= `minwave' & _total_complete_waves < ., /// 
        add(`add') j(`j') i(`i') wave(`wave') ///
        donordata(`donorid_data') diagnosticdata(`diagnosticdata')
    
    quietly : gen _PanelPatchImputed = `touse' == 1 & _total_complete_waves >= `minwave' & _total_complete_waves < .
    label var _PanelPatchImputed "Observation imputed or complete re PanelPatch"
    quietly : mi update 

    foreach v in `vunorderedvars' `sunorderedvars' {
        capture drop `v'_cat_*
    }

    if "`diagnosticdata'" != ""  {
        save "`diagnosticdata'_finalwaveimpute.dta", replace 
        local diagdata `"`diagdata'"' `"`diagnosticdata'_finalwaveimpute.dta"'
    }
    
    if "`weightvar'" != "" {
        if  "`_weightadjust_file'" == "" | ///
        `_interim_save_time' >= `_predictions_save_time' | ///
        `_predictions_save_time' >= `_donorcluster_save_time' | ///
        `_interim_save_time' >= `_donorcluster_save_time' | ///
        `_donorcluster_save_time' >= `_donoriddata_save_time' | ///
        `_predictions_save_time' >= `_donoriddata_save_time' | ///
        `_interim_save_time' >= `_donoriddata_save_time' | ///
        `_donorcluster_save_time' >= `_weightadjust_save_time' | ///
        `_predictions_save_time' >= `_weightadjust_save_time' | ///
        `_interim_save_time' >= `_weightadjust_save_time' | ///
        `_donoriddata_save_time' >= `_weightadjust_save_time' {
            
            preserve

            if "`showtimes'" != "" di _newline _newline c(current_date) _col(15) c(current_time) _col(25) "Weight Adjustment Starts"

            if "`xvars'" != "" {
                local waxvars "xvars(`xvars')"
            }

            tempfile weightadjust_data            
            

            di _newline _newline "Starting Weight Adjustment"
   
   **?? JW 06/02/2023: fixed hardcoded Baseweight variable
            PanelPatch_weight_adjust `weightvar', j(`j') i(`i') wave(`wave') ///
                gen(wgtadj) `waxvars' ///
                interimfile("`interim_imputations'") ///
                snumericvars(`snumericvars') ///
                sunorderedvars(`sunorderedvars') ///
                sorderedvars(`sorderedvars')  ///
                vnumericvars(`vnumericvars') ///
                vunorderedvars(`vunorderedvars') ///
                vorderedvars(`vorderedvars') ///
                lassofolds(`lassofolds') ///
                burnin(`burnin') ///
                saving("`weightadjust_data'")
            
            if "`diagnosticdata'" != ""  {
                quietly : use `weightadjust_data', clear
                save "`diagnosticdata'_weightadjust.dta", replace 
                local diagdata `"`diagdata'"' `"`diagnosticdata'_weightadjust.dta"'

            }

            restore
        
        }
        else {
            preserve
            use "`diagnosticdata'_weightadjust.dta", clear
            tempfile weightadjust_data
            quietly : save `weightadjust_data', replace
            local diagdata `"`diagdata'"' `"`diagnosticdata'_weightadjust.dta"'
            restore 
        }

        quietly : merge 1:1 `i' `j' `wave' using `weightadjust_data', nogen


        drop _total_complete_waves __* _touse
        
        foreach v in `vunorderedvars' `sunorderedvars' {
            capture drop `v'_cat_*
        }

        if "`diagnosticdata'" != ""  {
            save "`diagnosticdata'_final_weightadjust.dta", replace 
        }
        

    }
    capture drop _total_complete_waves __* _touse
    foreach v in `vunorderedvars' `sunorderedvars' {
        capture drop `v'_cat_*
    }
    
    if "`rundiagnostic'" != "" & "`weightvar'" != "" {
        foreach opt in weightvar snumericvars sunorderedvars vnumericvars ///
        vunorderedvars vorderedvars {
            if "" != "`opt'" {
                local diagcode `diagcode' `opt'(``opt'')
            }
        }
        PanelPatchDiag ///
            , i(`i') j(`j') wave(`wave') ///
            waveresponseflag(`waveresponseflag') `diagcode'
    }

end

