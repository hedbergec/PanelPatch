program define PanelPatch_predict_wrapper, rclass
    version 17 //written for stata 17
    syntax varlist (fv) [if], model(string) [xalways(string)] [wave(varlist)] ///
        [lassofolds(real 10)] [onlylasso] [replace] [onlypredictions]
    marksample touse, novarlist
    if "`onlylasso'" != "" {
        capture assert regexm("`model'", "mlogit") == 0
        if _rc != 0 {
            di as error "PanelPatch_predict_wrapper: onlylasso option is for linear only"
        }
    }
    
    if "" != "`replace'" {
        local outvar_drop_cmd "*"
    }
    else {
        local outvar_drop_cmd "drop"
    }
    if "`xalways'" != "" {
        local xalwaysmod "(`xalways')"
    }

    tokenize `varlist'
    local y `1'
    macro shift
    local x `*'
    local x_use ""
    if "`wave'" != "" {
        local x_use `x_use'  i.`wave' c.`wave'
        local wavemodopen "c.`wave'##("
        local wavemodclose ")"
    }
    local x_maybe ""
    foreach v in `x' {
        if regexm("`v'","\.") == 0 {
            local x_maybe `x_maybe' c.`v'
        }
        else {
            local x_maybe `x_maybe' `v'
        }
    }
    if regexm("`model'","reg") == 1 {
        capture `outvar_drop_cmd' `y'_1
        local y_list = 1
        capture confirm variable `y'_1
        if _rc == 0 {
            replace `y'_1 = `y' if `touse' & `y' < .
        }
        else {
            gen `y'_1 = `y' if `touse' & `y' < .
        }
        lasso linear `y'_1 `xalwaysmod' `wavemodopen'`x_maybe'`wavemodclose' if `touse', folds(`lassofolds')
        estimates store model
        local x_use = e(allvars_sel)
        local predictions "`y'_hat_1"
    }
    else if regexm("`model'","ologit") == 1 {
        levelsof `y' if `touse', local(y_list)
        foreach l in `y_list' {
            capture `outvar_drop_cmd' `y'_`l'
            capture confirm variable `y'_`l'
            if _rc == 0 {
                replace `y'_`l' = `y' == `l' if `touse' & `y' < .
            }
            else {
                gen `y'_`l' = `y' == `l' if `touse' & `y' < .
            }
            local predictions `predictions' `y'_hat_`l'
        }
        lasso linear `y'  `xalwaysmod' `wavemodopen'`x_maybe'`wavemodclose'  if `touse'
        local x_use = e(allvars_sel)
        estimates store model
    }
    else {
        quietly : levelsof `y' if `touse', local(y_list)
        foreach l in `y_list' {
            capture `outvar_drop_cmd' `y'_`l'
            capture confirm variable `y'_`l'
            if _rc == 0 {
                replace `y'_`l' = `y' == `l' if `touse' & `y' < .
            }
            else {
                gen `y'_`l' = `y' == `l' if `touse' & `y' < .
            }
            lasso linear `y'_`l'  `xalwaysmod' `wavemodopen'`x_maybe'`wavemodclose'  if `touse'
            local x_use_`l' = e(allvars_sel)
            local x_use : list x_use | x_use_`l'
            local predictions `predictions' `y'_hat_`l'
        }
    }
    
    capture drop `y'_hat_*
    if "`onlylasso'" == "" {
        `model' `y' `x_use' if `touse'
        estimates store model
    }
    estimates restore model
    predict `predictions' if `touse'
    foreach l in `y_list' {
        if "`onlypredictions'" != "" {
            replace `y'_`l' = `y'_hat_`l' if `touse' 
        }
        else {
            replace `y'_`l' = `y'_hat_`l' if `touse' & `y'_`l' >= .
        }
        
    }
    tokenize `y_list'
    local first_y `1'
    local out_vars ""
    foreach l in `y_list' {
        if `l' == `first_y' {
            if regexm("`model'","reg") == 1 {
                local out_vars `out_vars' `y'_`l'
                drop `y'_hat_`l'
            }
            else {
                drop `y'_`l' `y'_hat_`l'
            }
        }
        else {
            local out_vars `out_vars' `y'_`l'
            drop `y'_hat_`l'   
        }
    }
    return local out_vars = "`out_vars'"

end

