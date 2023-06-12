program define PanelPatch_corr, rclass
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
        }
        else {
            local prefix "mi estimate : "
        }
        local genprefix "mi xeq :"
        local matname "b_mi"
    }
    else {
        local prefix ""
        local genprefix ""
        local matname "b"
    }

    tokenize `varlist' 
    local x = "`1'"
    local y = "`2'"

    tempvar x2 y2 xy 

    quietly : `genprefix' gen `x2' = `x'^2
    quietly : `genprefix' gen `y2' = `y'^2
    quietly : `genprefix' gen `xy' = `x'*`y'

    quietly : `prefix' mean `x' `y' `x2' `y2' `xy' if `touse'
    tempname moments
    matrix `moments' = e(`matname')
    return scalar rho = (`moments'[1,5]-`moments'[1,1]*`moments'[1,2])/(sqrt(`moments'[1,3]-`moments'[1,1]^2)*sqrt(`moments'[1,4]-`moments'[1,2]^2))


end

