program define PanelPatch_running_response
    syntax varlist [if] , i(varlist) j(varlist) wave(varlist) GENerate(namelist)
    marksample touse 
    capture drop `generate'
    levelsof `wave' if `touse', local(wavelist)
    local numberwaves = wordcount("`wavelist'")
    tokenize `wavelist'
    gen `generate' = `varlist' if `wave' == `1' & `touse'
    forvalues w = 2/`numberwaves' {
        tempvar nextrow currentrow
        gen `currentrow' = `varlist' if `wave' < ``w'' & `touse'
        bysort `i' `j' : egen `nextrow' = min(`currentrow') if `touse'
        replace `generate' = `nextrow' if `wave' == ``w'' & `touse'
        replace `generate' = 0 if `generate' >= . & `wave' == ``w'' & `touse'
        drop `nextrow' `currentrow'
    }

end