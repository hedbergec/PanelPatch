program define mi_impute_cmd_PanelPatch_init //initialization program for PanelPatch

    preserve 

    forvalues m = 1/$MI_IMPUTE_user_opt_add {
        
        use "$MI_IMPUTE_user_donordata" if _mj == `m' & _donormarker == 1 , clear
        
        drop _mj _donormarker _imputemarker
        local diagpath "$MI_IMPUTE_user_diagnosticdata"
        local savepath = "`diagpath'" + "PanelPatch_donors_`m'.dta"

        save "`savepath'", replace
        
        global MI_IMPUTE_user_donor_`m' "`savepath'"

        use "$MI_IMPUTE_user_donordata" if _mj == `m', clear

        keep _donorid _donorcluster _donormarker _imputemarker $MI_IMPUTE_user_j $MI_IMPUTE_user_i $MI_IMPUTE_user_wave
        local diagpath ""
        local savepath ""
        local diagpath "$MI_IMPUTE_user_diagnosticdata"
        local savepath = "`diagpath'" + "PanelPatch_ids_`m'.dta"
        keep if _donormarker == 1 | _imputemarker == 1
        save "`savepath'", replace

        global MI_IMPUTE_user_ids_`m' "`savepath'"

    }

    restore 

end








 

	
	
	









