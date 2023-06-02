program define PanelPatch_mk_compvars, rclass
    syntax varlist [if], [categorical] id(varlist) wave(varlist) [limit(real 0.1)]
    marksample touse, novarlist
    levelsof `wave' if `touse', local(wavelist)
    local waverange = wordcount("`wavelist'") - 1
    duplicates report `varlist' if `touse'
    return list 
    if r(unique_value) < 3 {
        local categorical ""
    }
    if "`categorical'" == ""  {
        capture confirm variable `varlist'_lag 
        if _rc == 0 {
            drop `varlist'_lag 
        }
        capture confirm variable `varlist'_lead 
        if _rc == 0 {
            drop `varlist'_lead
        }
        capture confirm variable `varlist'_intrp
        if _rc == 0 {
            drop `varlist'_intrp
        }
        tempvar lag lead mulag mulead mulag_use mulead_use firstlag firstlead

        
        gen `lag' = .
        gen `lead' = .

        forvalues l = 1/`waverange' {
            bysort `id' (`wave') : replace `lag' = `varlist'[_n-`l'] if `touse' & `lag' == .
            bysort `id' (`wave') : replace `lead' = `varlist'[_n+`l'] if `touse' & `lead' == .
            if `l' == 1 {
                gen `firstlag' = 1 if `lag' ! = .
                gen `firstlead' = 1 if `lead' != .
            }
        }
        bysort `wave' : egen `mulag' = mean(`lag') if `touse' & `firstlag' == 1
        bysort `wave' : egen `mulead' = mean(`lead') if `touse' & `firstlead' == 1
        bysort `wave' : egen `mulag_use' = mean(`mulag') if `touse' 
        bysort `wave' : egen `mulead_use' = mean(`mulead') if `touse' 
        replace `lag' = `mulag_use' if `touse' & `lag' == .
        replace `lead' = `mulead_use' if `touse' & `lead' == .

        gen `varlist'_lag = `lag'
        gen `varlist'_lead = `lead'
        egen `varlist'_intrp = rowmean(`varlist'_lag `varlist'_lead)

        drop `lag' `lead' `mulag' `mulead' `mulag_use' `mulead_use' `firstlag' `firstlead'
        return local lags = "`varlist'_lag"    
        return local leads = "`varlist'_lead"  
        return local intrps = "`varlist'_intrp"  
    }
    else {
        local return_varlist_lag  ""
        local return_varlist_lead  ""
        local return_varlist_intrp  ""
        local return_varlist_cat  ""
        levelsof `varlist' if `touse', local(catlist)
        tokenize "`catlist'"
        macro shift
        local catlist `*'
        local number_categories = wordcount("`catlist'")
        tokenize "`catlist'"
        foreach cat in `catlist' {
            capture confirm variable `varlist'_`cat'_lag 
            if _rc == 0 {
                drop `varlist'_`cat'_lag 
            }
            capture confirm variable `varlist'_`cat'_lead 
            if _rc == 0 {
                drop `varlist'_`cat'_lead
            }
            capture confirm variable `varlist'_`cat'_intrp 
            if _rc == 0 {
                drop `varlist'_`cat'_intrp
            }
            capture confirm variable `varlist'_cat_`cat' 
            if _rc == 0 {
                drop `varlist'_cat_`cat' 
            }
            tempvar lag lead mulag mulead mulag_use mulead_use firstlag firstlead

            gen `lag' = .
            gen `lead' = .
            gen `varlist'_cat_`cat' = `varlist' == `cat' if `touse' & `varlist' < .
            local return_varlist_cat = "`return_varlist_cat' `varlist'_cat_`cat'"
            forvalues l = 1/`waverange' {
                bysort `id' (`wave') : replace `lag' = `varlist'_cat_`cat'[_n-`l'] == 1 if `touse' & `lag' == . & `varlist'_cat_`cat'[_n-`l'] < .
                bysort `id' (`wave') : replace `lead' = `varlist'_cat_`cat'[_n+`l'] == 1 if `touse' & `lead' == . & `varlist'_cat_`cat'[_n+`l'] < .
                if `l' == 1 {
                    gen `firstlag' = 1 if `lag' ! = .
                    gen `firstlead' = 1 if `lead' != .
                }
            }
            bysort `wave' : egen `mulag' = mean(`lag') if `touse' & `firstlag' == 1
            bysort `wave' : egen `mulead' = mean(`lead') if `touse' & `firstlead' == 1
            bysort `wave' : egen `mulag_use' = mean(`mulag') if `touse' 
            bysort `wave' : egen `mulead_use' = mean(`mulead') if `touse' 
            replace `lag' = `mulag_use' if `touse' & `lag' == .
            replace `lead' = `mulead_use' if `touse' & `lead' == .

            gen `varlist'_`cat'_lag = `lag'
            gen `varlist'_`cat'_lead = `lead'

            drop `lag' `lead' `mulag' `mulead' `mulag_use' `mulead_use' `firstlag' `firstlead'
            local egen_list ""
            sum `varlist'_`cat'_lag, meanonly
            if r(mean) < `limit' | r(mean) > 1-`limit'  {
                drop `varlist'_`cat'_lag
            }
            else {
                local return_varlist_lag "`return_varlist_lag' `varlist'_`cat'_lag"
                local egen_list `egen_list' `varlist'_`cat'_lag
            }
            sum `varlist'_`cat'_lead, meanonly
            if r(mean) < `limit' | r(mean) > 1-`limit' {
                drop `varlist'_`cat'_lead
            }
            else {
                local return_varlist_lead "`return_varlist_lead' `varlist'_`cat'_lead"
                local egen_list `egen_list' `varlist'_`cat'_lead
            }
            if "`egen_list'" != "" {
                local return_varlist_intrp "`return_varlist_intrp' `varlist'_`cat'_intrp"
                egen `varlist'_`cat'_intrp = rowmean(`egen_list')
            }
            
        }
        return local lags = "`return_varlist_lag'"    
        return local leads = "`return_varlist_lead'" 
        return local intrps = "`return_varlist_intrp'" 
        return local cats = "`return_varlist_cat'"    
    }


      


end




