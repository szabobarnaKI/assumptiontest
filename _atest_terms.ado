*! _atest_terms 0.1.0  classify the terms of the last estimation
*! Internal helper of assumptiontest. Reads the column names of e(b) (first
*! equation only). Returns
*!   r(cont) plain continuous variables (not categorical, not spline)
*!   r(cat)  variables entered as factors (i.var, 2.var, ...)
*!   r(spl)  variables belonging to a spline() group (as declared by the user)
*!   r(inter) 1 if interaction terms are present
*!   r(terms) design terms of the first equation (base/omitted levels removed)
*! Spline groups declared by the user are in globals ATEST_nsp, ATEST_sp1 ...
program define _atest_terms, rclass
    version 17

    tempname b
    matrix `b' = e(b)
    local cn : colnames `b'
    local ce : coleq `b'
    local eq1 : word 1 of `ce'
    local k = colsof(`b')

    local cont ""
    local cat  ""
    local terms ""
    local inter 0
    forvalues j = 1/`k' {
        local t : word `j' of `cn'
        local q : word `j' of `ce'
        if "`q'" != "`eq1'" continue
        if "`t'" == "_cons" continue
        * design columns for influence calculations: skip base/omitted levels
        if !ustrregexm("`t'", "[0-9]+[bo]\.") local terms `terms' `t'
        if strpos("`t'", "#") > 0 {
            local inter 1
            continue
        }
        local d = strpos("`t'", ".")
        if `d' > 0 {
            local pre = substr("`t'", 1, `d' - 1)
            local nm  = substr("`t'", `d' + 1, .)
            if "`pre'" == "c" {
                local cont `cont' `nm'
            }
            else {
                local cat `cat' `nm'
            }
        }
        else {
            local cont `cont' `t'
        }
    }
    local cont : list uniq cont
    local cat  : list uniq cat

    * ---- spline groups
    local spl ""
    forvalues g = 1/$ATEST_nsp {
        local spl `spl' ${ATEST_sp`g'}
    }
    local spl : list uniq spl
    local cont : list cont - spl

    return local cont  `cont'
    return local cat   `cat'
    return local spl   `spl'
    return local terms `terms'
    return scalar inter = `inter'
end
