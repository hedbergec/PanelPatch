program define PanelPatch_corr, rclass
    syntax varlist [if], [weightvar(varlist)] 
    marksample touse, novarlist
    capture mi describe
    if _rc == 0 {
        local setcommand "mi "
        local prefix "mi estimate : "
    }
    if "" != "`weightvar'" {
        local setcommand `setcommand' 
        `setcommand' svyset [pw = `weightvar']
        local prefix "`prefix' svy : "
    }
    local genprefix "mi xeq :"
    local matname "b_mi"
    tokenize `varlist' 
    local x = "`1'"
    local y = "`2'"

    tempvar x2 y2 xy 

    quietly : `genprefix' gen double `x2' = `x'^2
    quietly : `genprefix' gen double `y2' = `y'^2
    quietly : `genprefix' gen double `xy' = `x'*`y'

    quietly : `prefix' mean `x' `y' `x2' `y2' `xy' if `touse'
    tempname moments
    matrix `moments' = e(`matname')
    return scalar rho = (`moments'[1,5]-`moments'[1,1]*`moments'[1,2])/(sqrt(`moments'[1,3]-`moments'[1,1]^2)*sqrt(`moments'[1,4]-`moments'[1,2]^2))

    

    

end

PanelPatch_corr ed1 ed2, weightvar(Baseweight_wgtadj) j(PreschoolID )