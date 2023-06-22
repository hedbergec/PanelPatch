program define PanelPatchDiag
    syntax, i(varlist) j(varlist) wave(varlist) ///
        WAVEResponseflag(varlist) ///
        [weightvar(varlist)] /// 
        [SNumericvars(varlist)] ///
        [SUnorderedvars(varlist)] ///
        [SOrderedvars(varlist)]  ///
        [VNumericvars(varlist)] ///
        [VUnorderedvars(varlist)] ///
        [VOrderedvars(varlist)] 


    *Setting 
    mi svyset `j' [pw = `weightvar'_wgtadj]

    *WaveImpute Diagnostic Tables

    { //***** Table 1 Imputed due to wave nonresponse
    qui table `wave' if !missing(_donorid_master) & `waveresponseflag' == 0

    qui collect title "Table 1. Study Members with Imputed Waves due to Wave Nonresponse"
    collect layout
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) replace
    } //T1

    { //***** Table 2 Weighted Means for Repeating Numeric Variables Pre- and Post-Imputation
    tempvar wmvVNUM
    foreach v in `vnumericvars'{ //repeating numeric variables
        qui mi estimate: svy: mean `v', over(`wave') //imputed means
        mat matPOST`v' = r(table) //post matrix for variable
        
        tempname fr
        qui pwf
        frame copy `r(currentframe)' `fr'
        frames `fr'{ //changing frame to unset mi without deleting in working frame
        
        cap mi unset
        svyset `j' [pw = `weightvar']
        qui svy: mean `v', over(`wave') //non-imputed means using original sampling weight
        mat matPRE`v' = r(table) //pre matrix for variable
        } //end frame
        
        mat `wmvVNUM' = nullmat(`wmvVNUM')\matPRE`v'[1..2,1...]\matPOST`v'[1..2,1...] //matrix only keeping means and SDs for pre and post
    } //v; 

    *create macro for column names
    cap macro drop _cn
    quietly : glevelsof `wave', local(wavelist)
    foreach w in `wavelist'{
        local cn `"`cn' "Wave `w'""'
        local cn: list clean cn
    } //w; levels of wave

    mat colnames `wmvVNUM' = `cn' //label columns

    *create macro for row names
    cap macro drop _rn
    foreach var in `vnumericvars'{ //using var name as "equation" for labeling
        local rn = `"`rn' "`var'(pre):Weighted Mean" "`var'(pre):Clustered SE" "`var'(post):Weighted Mean" "`var'(post):Clustered SE""'
        local rn: list clean rn
    } // var; variable list

    mat rownames `wmvVNUM' = `rn' //label rows

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`wmvVNUM')", st_matrix("`wmvVNUM'"))
    mata : st_matrixrowstripe("r(`wmvVNUM')", st_matrixrowstripe("`wmvVNUM'"))
    mata : st_matrixcolstripe("r(`wmvVNUM')", st_matrixcolstripe("`wmvVNUM'"))

    qui collect clear
    qui collect get r(`wmvVNUM')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect style header roweq, title(hide) level(value)
    qui collect title "Table 2. Weighted Mean Values for Binary and Interval-Valued Repeating Variables Pre- and Post-Imputation"
    collect layout (roweq#rowname)(colname#cmdset)
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append
    } //T2

    { //***** Table 3 Weighted Frequency for Binary and Categorical Repeating Variables Pre- and Post-Imputation
    cap mat drop `wfVCAT'
    cap drop `cons'
    tempvar wfVCATmini wfVCATwide wfVCAT cons
    qui gen `cons' = 1

    foreach v in `snumericvars' `sorderedvars' `sunorderedvars' `vnumericvars' `vorderedvars' `vunorderedvars'{
        qui glevelsof `v'
        if `r(J)' < 10{
            local binCAT = "`binCAT' `v'"
            local binCAT: list clean binCAT
        } //if fewer than 10 levels
    } //v; all vars

    foreach v in `binCAT'{
        mi est: svy: total `cons', over(`wave' `v')
        mat matPOST`v' = r(table)
    } //v; categorical repeating vars

    tempname fr
    qui pwf
    frame copy `r(currentframe)' `fr'
    frames `fr'{ //changing frame to unset mi without deleting in working frame
        cap mi unset
        svyset `j' [pw = `weightvar']
        
        foreach v in `binCAT'{
            svy: total `cons', over(`wave' `v')
            mat matPRE`v' = r(table)
        }
    } //frames; using svyset in another frame

    foreach v in `binCAT'{
        local start = 1
        forv row = 1/`=colsof(matPRE`v')'{
            mat `wfVCATmini' = nullmat(`wfVCATmini')\(matPRE`v'[1,`row']\(matPRE`v'[2,`row']\matPOST`v'[1,`row'])\matPOST`v'[2,`row'])
        } //row; number of rows in matrix
        
        qui glevelsof `v'
        local lvlCNT = `r(J)'
        qui glevelsof `wave'
        local wvCNT = `r(J)'
        forv wv = 1/`wvCNT'{
            local end = `wv'*`lvlCNT'*4
            mat `wfVCATwide' = nullmat(`wfVCATwide'),(`wfVCATmini'[`start'..`end', 1])
            local start = `end' + 1
        } //wv; levels of wave
        
        mat `wfVCAT' = nullmat(`wfVCAT')\(`wfVCATwide')
        cap mat drop `wfVCATwide' `wfVCATmini'
    }  //v; repeating categorical variables

    *create macro for column names
    cap macro drop _cn
    quietly : glevelsof `wave', local(wavelist)
    foreach w in `wavelist'{
        local cn `"`cn' "Wave `w'""'
        local cn: list clean cn
    } //w; levels of wave

    mat colnames `wfVCAT' = `cn' //label columns

    *create macro for row names
    cap macro drop _rn
    foreach var in `binCAT'{ //using var name as "equation" for labeling
        qui glevelsof `var', local(lvls)
        foreach lvl in `lvls'{
            local rn = `"`rn' "`var'_`lvl' (pre):Weighted Frequency" "`var'_`lvl'(pre):Clustered SE" "`var'_`lvl'(post):Weighted Frequency" "`var'_`lvl'(post):Clustered SE""'
        } //lvl; level of cat var
        local rn: list clean rn
    } // var; variable list

    mat rownames `wfVCAT' = `rn' //label rows

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`wfVCAT')", st_matrix("`wfVCAT'"))
    mata : st_matrixrowstripe("r(`wfVCAT')", st_matrixrowstripe("`wfVCAT'"))
    mata : st_matrixcolstripe("r(`wfVCAT')", st_matrixcolstripe("`wfVCAT'"))

    qui collect clear
    qui collect get r(`wfVCAT')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect style header roweq, title(hide) level(value)
    qui collect title "Table 3. Weighted Frequency for Categorical Repeating Variables Pre- and Post-Imputation"
    collect layout (roweq#rowname)(colname#cmdset)
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append

    } //T3

    { //***** Table 4 Cross-Wave Linear and Quadratic Trends in Binary and Interval-Valued Variables Pre- and Post-Imputation
    tempvar cwlqNUM

    foreach v in `snumericvars' `vnumericvars'{
        qui mi est (_b[`wave']): svy: reg `v' `wave'
        mat postL`v' = r(table)[1,1]\r(table)[2,1]
        qui eststo postQ`i': mi est (_b[`wave'] + _b[c.`wave'#c.`wave']): svy: reg `v' c.`wave'##c.`wave'
        mat postQ`v' = r(table)[1,1]\r(table)[2,1]
        
        tempname fr`v'
        qui pwf
        frame copy `r(currentframe)' `fr`v''

        frames `fr`v''{ //changing frame to unset mi without deleting in working frame
            cap mi unset
            qui svyset `j' [pw = `weightvar']

            qui svy: reg `v' `wave'
            qui lincom _b[`wave']
            mat preL`v' = r(estimate)\r(se)
            
            qui svy: reg `v' c.`wave'##c.`wave'
            qui lincom _b[`wave'] + _b[c.`wave'#c.`wave']
            mat preQ`v' = r(estimate)\r(se)
        } //frames; using svyset in another frame

        mat `cwlqNUM' = nullmat(`cwlqNUM')\(preL`v',postL`v',preQ`v',postQ`v')
    } //v; numeric vars

    *create macro for column names
    mat colnames `cwlqNUM' = "Linear (pre)" "Linear (post)" "Quadratic (pre)" "Quadratic (post)" //label columns

    *create macro for row names
    cap macro drop _rn
    foreach var in `snumericvars' `vnumericvars'{ //using var name as "equation" for labeling
        local rn = `"`rn' "`var':Trend Estimate" "`var':Clustered SE""'
        local rn: list clean rn
    } // var; variable list

    mat rownames `cwlqNUM' = `rn' //label rows

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`cwlqNUM')", st_matrix("`cwlqNUM'"))
    mata : st_matrixrowstripe("r(`cwlqNUM')", st_matrixrowstripe("`cwlqNUM'"))
    mata : st_matrixcolstripe("r(`cwlqNUM')", st_matrixcolstripe("`cwlqNUM'"))

    qui collect clear
    qui collect get r(`cwlqNUM')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect style header roweq, title(hide) level(value)
    qui collect title "Table 4. Cross-Wave Linear and Quadratic Trends in Binary and Interval-Valued Variables Pre- and Post-Imputation"
    collect layout (roweq#rowname)(colname#cmdset)
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append

    } //T4

    { //***** Table 5 Gamma Statistic for Ordered Categorical Repeating Variables with Wave
    tempvar gocr
    cap mat drop `gocr'
    foreach v in `vorderedvars'{
        PanelPatch_gamma WaveID `v', weightvar(`weightvar'_wgtadj)
        mat `gocr' = nullmat(`gocr')\(`r(gamma)')
    } //v; unordered repeating vars

    *column names
    mat colnames `gocr' = "Goodman and Kruskal's Gamma" //label columns

    *create macro for row names
    cap macro drop _rn
    foreach var in `vorderedvars'{ //using var name as "equation" for labeling
        local rn = `"`rn' "`var'""'
        local rn: list clean rn
    } // var; variable list

    macro li _rn
    mat rownames `gocr' = `rn' //label rows

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`gocr')", st_matrix("`gocr'"))
    mata : st_matrixrowstripe("r(`gocr')", st_matrixrowstripe("`gocr'"))
    mata : st_matrixcolstripe("r(`gocr')", st_matrixcolstripe("`gocr'"))

    qui collect clear
    qui collect get r(`gocr')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect style header roweq, title(hide) level(value)
    qui collect title "Table 5. Gamma Statistic between Wave and Ordered Categorical Variables"
    collect layout (rowname)(colname#cmdset)
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append
    } //T5

    { //***** Table 6 Cross-Wave Linear Trends in Each Level of Unordered Categorical Variables Pre- and Post-Imputation
    tempvar cwlCAT

    foreach v in `sunorderedvars' `vunorderedvars'{
        cap drop `v'?
        qui tab `v', gen(`v')
        
        qui glevelsof `v', local(lvls)
        foreach lvl in `lvls'{
            qui mi register imputed `v'`lvl'
        }

        qui mi query
        local M = `r(M)'
        forv m = 1/`M'{
            cap drop _`m'_`v'?
            qui tab _`m'_`v', gen(_`m'_`v')
        }
        
        foreach lvl in `lvls'{
            qui mi est (_b[`wave']): svy: reg `v'`lvl' `wave'
            mat postL`v'`lvl' = r(table)[1,1]\r(table)[2,1]
        }
        
        tempname fr`v'
        qui pwf
        frame copy `r(currentframe)' `fr`v''

        frames `fr`v''{ //changing frame to unset mi without deleting in working frame
            cap mi unset
            qui svyset `j' [pw = `weightvar']
            
            qui glevelsof `v', local(lvls)
            foreach lvl in `lvls'{
                qui svy: reg `v'`lvl' `wave' 
                qui lincom _b[`wave']
                mat preL`v'`lvl' = r(estimate)\r(se)
            }
        } //frames; using svyset in another frame
        
        qui glevelsof `v', local(lvls)
        foreach lvl in `lvls'{
            mat `cwlCAT' = nullmat(`cwlCAT')\(preL`v'`lvl',postL`v'`lvl')
        }
    } //v; numeric vars

    *column names
    mat colnames `cwlCAT' = "Linear (pre)" "Linear (post)" //label columns

    *create macro for row names
    cap macro drop _rn
    foreach var in `sunorderedvars' `vunorderedvars'{ //using var name as "equation" for labeling
        qui glevelsof `var', local(lvls)
        foreach lvl in `lvls'{
            local rn = `"`rn' "`var'(`lvl'):Trend Estimate" "`var'(`lvl'):Clustered SE""'
            local rn: list clean rn
        }
    } // var; variable list

    mat rownames `cwlCAT' = `rn' //label rows

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`cwlCAT')", st_matrix("`cwlCAT'"))
    mata : st_matrixrowstripe("r(`cwlCAT')", st_matrixrowstripe("`cwlCAT'"))
    mata : st_matrixcolstripe("r(`cwlCAT')", st_matrixcolstripe("`cwlCAT'"))

    qui collect clear
    qui collect get r(`cwlCAT')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect style header roweq, title(hide) level(value)
    qui collect title "Table 6. Cross-Wave Linear Trends in Each Level of Unordered Categorical Variables Pre- and Post-Imputation"
    collect layout (roweq#rowname)(colname#cmdset)
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append
    } //T6

    { //***** Table 7 The correlation of each interval-valued repeating variable and Wave before and after imputation
    local fWAVE = `:word 1 of `lvls''
    qui glevelsof `wave', local(lvls)
    local lWAVE: list sizeof local(lvls)
    local lWAVE = `:word `lWAVE' of `lvls''
    tempvar ivrCORR ivrCORRpost ivrCORRpre
    cap mat drop `ivrCORR' `ivrCORRpost' `ivrCORRpre'
    foreach v in `vnumericvars'{
        qui glevelsof `v'
        if r(J)>10{
            tempvar v1 v5
            qui gen `v1' = `v' if `wave' == `fWAVE'
            qui bys `i': ereplace `v1' = max(`v1')
            qui gen `v5' = `v' if `wave' == `lWAVE'
            qui bys `i': ereplace `v5' = max(`v5')
            PanelPatch_corr `v1' `v5' if `wave' == `fWAVE', weightvar(`weightvar'_wgtadj)
            mat `ivrCORRpost' = nullmat(`ivrCORRpost')\(`r(rho)')
            tempname fr
            qui pwf
            frame copy `r(currentframe)' `fr'
            frames `fr'{ //changing frame to unset mi without deleting in working frame
                cap mi unset
                svyset `j' [pw = `weightvar']
                
                tempvar x y x2 y2 xy 
                    
                qui: gen double `x' = `v' if `wave' == `fWAVE'
                qui bys `i': ereplace `x' = max(`x')
                quietly : gen double `x2' = `x'^2
                qui: gen double `y' = `v' if `wave' == `lWAVE'
                qui bys `i': ereplace `y' = max(`y')
                quietly : gen double `y2' = `y'^2
                quietly : gen double `xy' = `x'*`y'

                quietly : svy: mean `x' `y' `x2' `y2' `xy' if `wave' == `fWAVE'
                tempname moments
                matrix `moments' = e(b)			
                mat `ivrCORRpre' = nullmat(`ivrCORRpre')\(`moments'[1,5]-`moments'[1,1]*`moments'[1,2])/(sqrt(`moments'[1,3]-`moments'[1,1]^2)*sqrt(`moments'[1,4]-`moments'[1,2]^2))
            } // end frame
        } //if more than 10 categories, treat as interval-valued
    } //v; numeric repeating vars

    mat `ivrCORR' = (`ivrCORRpre'), (`ivrCORRpost')


    *column names
    mat colnames `ivrCORR' = "Pre" "Post" //label columns

    *create macro for row names
    cap macro drop _rn
    foreach var in `vnumericvars'{
        qui glevelsof `var'
        if r(J)>10{//using var name as "equation" for labeling
            local rn = `"`rn' "`var'""'
            local rn: list clean rn
        } //if more than 10 categories, treat as interval-valued
    } //v; numeric repeating vars

    mat rownames `ivrCORR' = `rn' //label rows

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`ivrCORR')", st_matrix("`ivrCORR'"))
    mata : st_matrixrowstripe("r(`ivrCORR')", st_matrixrowstripe("`ivrCORR'"))
    mata : st_matrixcolstripe("r(`ivrCORR')", st_matrixcolstripe("`ivrCORR'"))

    qui collect clear
    qui collect get r(`ivrCORR')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect style header roweq, title(hide) level(value)
    qui collect title "Table 7. Corralation between Wave and Interval-Valued Repeating Variables Pre- and Post-Imputation"
    collect layout (rowname)(colname#cmdset)
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append
    } //T7

    { //***** Table 8 R-squared for Linear Regression of Binary, Ordered Categorical, and Interval-Valued Variables on the set of Stable Variables in the Final Wave Pre- and Post-Imputation
    tempvar R2all preR2 postR2
    cap mat drop `R2all' `preR2' `postR2'
    qui glevelsof `wave', local(lvls)
    local fWAVE: list sizeof local(lvls)
    local fWAVE = `:word `fWAVE' of `lvls''
    local regCONT: list clean snumericvars
    local regCAT: list sordervars | sunorderedvars

    qui mi query
    local M = r(M)
    scalar r2 = 0

    foreach v in `snumericvars' `sorderedvars' `sunorderedvars' `vnumericvars' `vorderedvars' `vunorderedvars'{
        local regCONT: list regCONT - v
        local regCAT: list regCAT - v
        qui mi xeq 1/`M': svy: reg `v' `regCONT' b1.(`regCAT') if `wave' == `fWAVE'; scalar r2 = r2 + atanh(sqrt(e(r2)))
        scalar r2 = tanh(r2/`M')^2 //R2 using Fisher's z over imputed data
        mat `postR2' = nullmat(`postR2')\(`=r2')
        
        cap frames drop `fr`v''
        tempname fr`v'
        qui pwf
        frame copy `r(currentframe)' `fr`v''

        frames `fr`v''{ //changing frame to unset mi without deleting in working frame
            cap mi unset
            qui svyset `j' [pw = `weightvar']

            qui svy: reg `v' `regCONT' b1.(`regCAT') if `wave' == `fWAVE'
            mat `preR2' = nullmat(`preR2')\e(r2)
        } //frames; using svyset in another frame
        
    } // v; all binary, ordered cat, and continuous vars

    mat `R2all' = (`preR2'), (`postR2')

    *column names
    mat colnames `R2all' = "Pre" "Post" //label columns

    *create macro for row names
    cap macro drop _rn
    foreach var in `snumericvars' `sorderedvars' `sunorderedvars' `vnumericvars' `vorderedvars' `vunorderedvars'{ //using var name as "equation" for labeling
        local rn = `"`rn' "`var'""'
        local rn: list clean rn
    } // var; variable list

    mat rownames `R2all' = `rn' //label rows

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`R2all')", st_matrix("`R2all'"))
    mata : st_matrixrowstripe("r(`R2all')", st_matrixrowstripe("`R2all'"))
    mata : st_matrixcolstripe("r(`R2all')", st_matrixcolstripe("`R2all'"))

    qui collect clear
    qui collect get r(`R2all')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect title "Table 8. R-squared for Linear Regression of Variables on the set of Stable Variables in the Final Wave Pre- and Post-Imputation"
    collect layout (rowname)(colname[Pre Post])
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append
    } //T8

    { //***** Table 9 Kish Design Effect Pre- and Post-Weight Adjustment
    qui glevelsof `wave', local(lvls)
    local fWAVE: list sizeof local(lvls)
    local fWAVE = `:word `fWAVE' of `lvls''

    qui sum `weightvar' if `wave' == `fWAVE'
    local varWGT = r(sd)*r(sd)
    local meanWGT = r(mean)
    local kishDEFF = 1 + (`varWGT'/(`meanWGT'*`meanWGT'))

    qui sum `weightvar'_wgtadj if `wave' == `fWAVE' & !missing(_donorid_master)
    local varWGTadj = r(sd)*r(sd)
    local meanWGTadj = r(mean)
    local kishDEFFadj = 1 + (`varWGTadj'/(`meanWGTadj'*`meanWGTadj'))

    tempvar kD
    mat `kD' = `kishDEFF', `kishDEFFadj'

    *column names
    mat colnames `kD' = "Original Weight" "Adjusted Weight" //label columns

    *row names
    mat rownames `kD' = "Kish Design Effect"

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`kD')", st_matrix("`kD'"))
    mata : st_matrixrowstripe("r(`kD')", st_matrixrowstripe("`kD'"))
    mata : st_matrixcolstripe("r(`kD')", st_matrixcolstripe("`kD'"))

    qui collect clear
    qui collect get r(`kD')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect style header roweq, title(hide) level(value)
    qui collect title "Table 9. Kish Design Effect Pre- and Post-Weight Adjustment"
    collect layout (roweq#rowname)(colname#cmdset)
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append
    } //T9

    { //***** Table 10 Weighted Means for Repeating Numeric Variables Pre- and Post-Weight Adjustment
    tempvar awmvVNUM
    foreach v in `vnumericvars'{ //repeating numeric variables
        qui mi estimate: svy: mean `v', over(`wave') //imputed means
        mat matPOST`v' = r(table) //post matrix for variable
        
        cap frame drop `fr'
        tempname fr
        qui pwf
        frame copy `r(currentframe)' `fr'
        frames `fr'{ //changing frame to unset mi without deleting in working frame
            mi svyset `j' [pw = `weightvar']
            
            qui mi estimate: svy: mean `v', over(`wave') //imputed means
            mat matPRE`v' = r(table) //pre matrix for variable
        } //frames; using svyset in another frame

        mat `awmvVNUM' = nullmat(`awmvVNUM')\matPRE`v'[1..2,1...]\matPOST`v'[1..2,1...] //matrix only keeping means and SDs for pre and post
    } //v; 

    *create macro for column names
    cap macro drop _cn
    quietly : glevelsof `wave', local(wavelist)
    foreach w in `wavelist'{
        local cn `"`cn' "Wave `w'""'
        local cn: list clean cn
    } //w; levels of wave

    mat colnames `awmvVNUM' = `cn' //label columns

    *create macro for row names
    cap macro drop _rn
    foreach var in `vnumericvars'{ //using var name as "equation" for labeling
        local rn = `"`rn' "`var'(pre):Weighted Mean" "`var'(pre):Clustered SE" "`var'(post):Weighted Mean" "`var'(post):Clustered SE""'
        local rn: list clean rn
    } // var; variable list

    mat rownames `awmvVNUM' = `rn' //label rows

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`awmvVNUM')", st_matrix("`awmvVNUM'"))
    mata : st_matrixrowstripe("r(`awmvVNUM')", st_matrixrowstripe("`awmvVNUM'"))
    mata : st_matrixcolstripe("r(`awmvVNUM')", st_matrixcolstripe("`awmvVNUM'"))

    qui collect clear
    qui collect get r(`awmvVNUM')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect style header roweq, title(hide) level(value)
    qui collect title "Table 10. Weighted Means for Repeating Numeric Variables Pre- and Post-Weight Adjustment"
    collect layout (roweq#rowname)(colname#cmdset)
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append
    } //T10

    { //***** Table 11 Weighted Frequency for Binary and Categorical Repeating Variables Pre- and Post-Weight Adjustment
    cap mat drop `awfVCAT'
    cap drop `cons'
    tempvar awfVCATmini awfVCATwide awfVCAT cons
    qui gen `cons' = 1

    foreach v in `snumericvars' `sorderedvars' `sunorderedvars' `vnumericvars' `vorderedvars' `vunorderedvars'{
        qui glevelsof `v'
        if `r(J)' < 10{
            local binCAT = "`binCAT' `v'"
            local binCAT: list clean binCAT
        } //if fewer than 10 levels
    } //v; all vars

    foreach v in `binCAT'{
        mi est: svy: total `cons', over(`wave' `v')
        mat matPOST`v' = r(table)
    } //v; categorical repeating vars

    cap frame drop `fr'
    tempname fr
    qui pwf
    frame copy `r(currentframe)' `fr'
    frames `fr'{ //changing frame to unset mi without deleting in working frame
        mi svyset `j' [pw = `weightvar']
        
        foreach v in `binCAT'{
            mi est: svy: total `cons', over(`wave' `v')
            mat matPRE`v' = r(table)
        }
    } //frames; using svyset in another frame

    foreach v in `binCAT'{
        local start = 1
        forv row = 1/`=colsof(matPRE`v')'{
            mat `awfVCATmini' = nullmat(`awfVCATmini')\(matPRE`v'[1,`row']\(matPRE`v'[2,`row']\matPOST`v'[1,`row'])\matPOST`v'[2,`row'])
        } //row; number of rows in matrix
        
        qui glevelsof `v'
        local lvlCNT = `r(J)'
        qui glevelsof `wave'
        local wvCNT = `r(J)'
        forv wv = 1/`wvCNT'{
            local end = `wv'*`lvlCNT'*4
            mat `awfVCATwide' = nullmat(`awfVCATwide'),(`awfVCATmini'[`start'..`end', 1])
            local start = `end' + 1
        } //wv; levels of wave
        
        mat `awfVCAT' = nullmat(`awfVCAT')\(`awfVCATwide')
        cap mat drop `awfVCATwide' `awfVCATmini'
    }  //v; repeating categorical variables

    *create macro for column names
    cap macro drop _cn
    quietly : glevelsof `wave', local(wavelist)
    foreach w in `wavelist'{
        local cn `"`cn' "Wave `w'""'
        local cn: list clean cn
    } //w; levels of wave

    mat colnames `awfVCAT' = `cn' //label columns

    *create macro for row names
    cap macro drop _rn
    foreach var in `binCAT'{ //using var name as "equation" for labeling
        qui glevelsof `var', local(lvls)
        foreach lvl in `lvls'{
            local rn = `"`rn' "`var'_`lvl' (pre):Weighted Frequency" "`var'_`lvl'(pre):Clustered SE" "`var'_`lvl'(post):Weighted Frequency" "`var'_`lvl'(post):Clustered SE""'
        } //lvl; level of cat var
        local rn: list clean rn
    } // var; variable list

    mat rownames `awfVCAT' = `rn' //label rows

    *Convert matrix to rclass for -collect-
    mata : st_matrix("r(`awfVCAT')", st_matrix("`awfVCAT'"))
    mata : st_matrixrowstripe("r(`awfVCAT')", st_matrixrowstripe("`awfVCAT'"))
    mata : st_matrixcolstripe("r(`awfVCAT')", st_matrixcolstripe("`awfVCAT'"))

    qui collect clear
    qui collect get r(`awfVCAT')
    qui collect style header cmdset, title(hide) level(hide)
    qui collect style header colname, title(hide) level(value)
    qui collect style header rowname, title(hide) level(value)
    qui collect style header roweq, title(hide) level(value)
    qui collect title "Table 11. Weighted Frequency for Binary and Categorical Repeating Variables Pre- and Post-Weight Adjustment"
    collect layout (roweq#rowname)(colname#cmdset)
    qui collect export PanelPatch_DiagnosticTables.txt, as(txt) append
    } //T11
end