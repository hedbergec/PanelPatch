
program define PanelPatch_gamma, rclass
    syntax varlist [if], [weightvar(varlist)] [j(varlist)] [raw]
    marksample touse, novarlist
    capture mi describe
    if _rc == 0 {
        local setcommand "mi "
    }
    if "`raw'" == "" {
        if "" != "`weightvar'" {
            local setcommand `setcommand' 
            `setcommand' svyset `j' [pw = `weightvar']
            local prefix "mi estimate : svy : "
            quietly : sum `weightvar'
            local totalweight = round(r(sum))
            local sigdigits = strlen("`totalweight'") + 2
        }
        else {
            local prefix "mi estimate : "
            local totalweight = _N
            local sigdigits = strlen("`totalweight'") + 2
        }
        local matname "b_mi"
    }
    else {
        local prefix ""
        local totalweight = _N
        local sigdigits = strlen("`totalweight'") + 2
        local matname "b"
    }
    
    tempvar cons
    quietly gen `cons' = 1 if `touse'
    quietly : `prefix' total `cons', over(`varlist')
    tempname totals
    matrix `totals' = e(`matname')
    tokenize `varlist' 
    local x = "`1'"
    local y = "`2'"
    quietly : levelsof `x' if `touse', local(xlist) 
    quietly : levelsof `y' if `touse', local(ylist) 
    quietly : sum `y'
    local maxr = r(max)
    quietly : sum `x'
    local maxc = r(max)
    local tabispecs = " "
    foreach r in `ylist' {
        foreach c in `xlist' {
            local current_total : display %`sigdigits'.0f `totals'[1,"c.`cons'@`c'.`x'#`r'.`y'"]
            local tabispecs =   " `tabispecs' `current_total' "  
            if `c' == `maxc' & `r' != `maxr' {
                local tabispecs " `tabispecs' \ "  
            }
        }
    }
    *di "`tabispecs'"
    tabi `tabispecs', gamma
    return scalar gamma = r(gamma)
end

