program define mi_impute_cmd_PanelPatch //imputation routine
    version 17

    di _newline "Transferring values for imputation set $MI_IMPUTE_user_m"
    local m = $MI_IMPUTE_user_m
    quietly : merge 1:1 $MI_IMPUTE_user_j $MI_IMPUTE_user_i $MI_IMPUTE_user_wave ///
        using "${MI_IMPUTE_user_ids_`m'}", nogen update replace
    
    quietly : merge m:1 _donorid _donorcluster $MI_IMPUTE_user_wave ///
        using "${MI_IMPUTE_user_donor_`m'}", nogen update
    
    

end


