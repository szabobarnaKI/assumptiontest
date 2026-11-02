*! _atest_dfbeta 0.3.1  influence screening from predict, dfbeta (Cox models)
*! Internal helper of assumptiontest. Flags observations whose removal would
*! shift any coefficient by more than outcut() standard errors.
*! Base/omitted factor levels (zero variance) are skipped. The dfbeta variables
*! are matched to e(b) columns whether predict creates one per column of e(b)
*! or one per estimated coefficient only.
*! Must be called while a Cox model (stcox) is the active estimation.
program define _atest_dfbeta
    version 17
    syntax , OUTcut(real) [NOte(string)]
    if "$ATEST_nooutliers" != "" exit

    local outtxt = strtrim(string(`outcut', "%9.3g"))
    if substr("`outtxt'", 1, 1) == "." local outtxt "0`outtxt'"

    tempname V
    matrix `V' = e(V)
    local k = colsof(`V')
    local pos ""
    forvalues j = 1/`k' {
        if `V'[`j', `j'] > 0 & `V'[`j', `j'] < . local pos `pos' `j'
    }
    local npos : word count `pos'

    local stub __atest_dfb
    capture drop `stub'*
    capture quietly predict double `stub'*, dfbeta
    if _rc {
        local rc = _rc
        _atest_post, name(outliers) type(notrun) ///
            desc("Influential observations (dfbeta)") ///
            note("predict, dfbeta not available (rc=`rc')")
        capture drop `stub'*
        exit
    }
    capture quietly ds `stub'*
    local ndfb : word count `r(varlist)'
    if `ndfb' != `k' & `ndfb' != `npos' {
        _atest_post, name(outliers) type(notrun) ///
            desc("Influential observations (dfbeta)") ///
            note("`ndfb' dfbeta variables for `k' coefficients (`npos' estimated): cannot be matched")
        capture drop `stub'*
        exit
    }

    local maxinf 0
    tempvar flagv
    quietly generate byte `flagv' = 0 if $ATEST_esamp
    forvalues i = 1/`npos' {
        local col : word `i' of `pos'
        local vix = cond(`ndfb' == `k', `col', `i')
        local se = sqrt(`V'[`col', `col'])
        quietly replace `flagv' = 1 if abs(`stub'`vix' / `se') > `outcut' & `stub'`vix' < . & $ATEST_esamp
        quietly summarize `stub'`vix' if $ATEST_esamp, meanonly
        local m = max(abs(r(min)), abs(r(max))) / `se'
        if `m' > `maxinf' local maxinf `m'
    }
    quietly count if `flagv' == 1
    local nout = r(N)
    capture drop `stub'*
    if `nout' > 0 {
        _atest_post, name(outliers) flag stat(`maxinf') statname(max|dfbeta|/SE) ///
            desc("Influential observations: `nout' observation(s) shift a coefficient by > `outtxt' SE") ///
            note(`"`note'"')
    }
    else {
        _atest_post, name(outliers) stat(`maxinf') statname(max|dfbeta|/SE) ///
            desc("Influential observations: none shifts a coefficient by > `outtxt' SE") ///
            note(`"`note'"')
    }
end
