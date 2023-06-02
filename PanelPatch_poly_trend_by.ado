program define PanelPatch_poly_trend_by, rclass //compute intercept, time and time poly coeffs by j and stores them in the three variables noted in gen
    syntax varlist [if] [in], /// varlist is outcome
        Time(varlist) GENerate(namelist) /// time is wave variable name and gen is stem for new variables
        j(varlist) polymax(numlist)
    marksample touse, novarlist
    quietly : levelsof `time' if /// gather unique values of time variable if 
        `touse' & `varlist' < . & `time' < ., local(waves) //useable and non-missing, store in waves
    local terms = min(wordcount("`waves'"),`polymax' + 1) //number poly terms is unique wave values, topped at polymax + 1 for constant
    local return_list = "" //start local for output variables
    tokenize `generate'
    forvalues k = 1/`terms' { //error checking 0 is intercept
        capture drop ``k''
        local out_`k' = "``k''"
    }
    *trend expression
    local trendx ""
    forvalues k = 2/`terms' { //make string of time variable with c. prefix 
        local trendx `trendx' c.`time'
    }
    local trendx = subinstr("`trendx'", " ", "##", .) //replace spaces with double pound (main effect and interaction)
    tempvar idvar 
    quietly : egen `idvar' = group(`j') if `touse'
    quietly : levelsof `idvar' if ///get list of by-var values if
        `touse' & `varlist' < . & `time' < ., /// useable and varlist is non-missing and time is nonmissing
        local(by_list) //store in macro by_list
    
    if "`by_list'" != "" { //if there are any by vals
        forvalues k = 1/`terms' { //make output vars for terms, 0 is intercept
            local return_list `return_list' `out_`k''
            quietly : gen `out_`k'' = .
        }
        foreach i in `by_list' { //loop over each by value 
            quietly : count if `touse' & `idvar' == `i' & `varlist' < . & `time' < .
            if r(N) >= `polymax' & r(N) < . { //check for number of non-missin obs to fit model
                quietly : reg `varlist' `trendx' if `touse' & `idvar' == `i' //run ols
                tempname b 
                matrix `b' = e(b)
                forvalues k = 1/`terms' { //make output vars for terms, last is intercept
                    quietly : replace `out_`k'' = `b'[1,`k'] if `touse' & `idvar' == `i' 
                }
            }
        }
        return local return_list `return_list'
    }
    *if there were no by values, nothing done and return list is empty
end

