program define mi_impute_cmd_PanelPatch_parse //parsing program for abtlongdeck
    version 17 //written for stata 17
    syntax anything(equalok) /// will be a list of variables with an equal sign, left are things to be imputed, right are predictors
    [if] [, *] // if statment and options

    set type double //use double precision for all computations and new variable creation

    gettoken ivars xvars : anything, parse("=") //split the varaible list by equal sign
    gettoken eq xvars : xvars, parse("=") //get the predictors
    
    u_mi_impute_user_setup `if' , /// pass information to mi setup to post global macros
    ivars(`ivars') xvars(`xvars') `options' 
    set trace off
    /* 
    we need to parse our own options, the code below does this
    we declare which options are "ok", the use gettoken to grab each option(options) pattern
    the a quick check that each option grabbed is part of the OK list
    */

    local okoptions "i j wave donordata diagnosticdata"
    local the_options = "$MI_IMPUTE_user_options" //global macro passed by parse program
    quietly : gettoken left right : the_options, match(parns) bind //start grabbing the option(options) pattern
    while "`left'" != "" { //so long as there is something left (pun, see help gettoken for how it actually works) to parse
        if regexm("`left'", "^(.+)\(") local opt = regexs(1) //get what is before the (, which is the option name
        if regexm("`left'", "^.+\((.+)\)$") local opt_entries =  regexs(1) //get stuff in the ()
        local posof : list posof "`opt'" in okoptions //in the string okoptions, which one is it
        if `posof' != 0 { //if its on the list, its not pos=0, so post to global and local macros
            global MI_IMPUTE_user_`opt' = "`opt_entries'" //global macro
            local `opt' = "`opt_entries'" //local macro
            quietly : gettoken left right : right , match(parns) bind
        }
        else { //if its not on the list, error out
            di as error "option `opt' not allowed. options for abtdecklong are `okoptions'"
            error 198 
        }
    }

    

end




