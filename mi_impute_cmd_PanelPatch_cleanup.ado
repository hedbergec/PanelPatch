program define mi_impute_cmd_PanelPatch_cleanup
    drop _donorid _donorcluster 
    forvalues m = 1/$MI_IMPUTE_user_opt_add {
        erase `"${MI_IMPUTE_user_donor_`m'}"'
        erase `"${MI_IMPUTE_user_ids_`m'}"'
    }
end
