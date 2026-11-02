*! _atest_ld 0.5.1  influence screening with likelihood displacement (Cox models)
*! Internal helper of assumptiontest. predict, ldisplace gives the (approximate)
*! likelihood displacement LD_i of each observation: the amount by which the
*! log likelihood at the full-data estimate falls (x2) when the estimate is
*! replaced by the estimate without observation i. An observation is flagged
*! when LD_i exceeds the median of the chi-squared distribution with as many df
*! as there are estimated coefficients (base/omitted factor levels excluded),
*! i.e. when deleting it moves the estimate outside the 50% confidence region
*! of the full-data estimate (the Cook and Weisberg 1982 calibration, which is
*! also behind the Cook's distance > 1 rule). Independent of alpha().
*! ldcut() overrides the cut-off.
*! Must be called while a Cox model (stcox) is the active estimation.
program define _atest_ld
    version 17
    syntax [, NOte(string)]
    if "$ATEST_nooutliers" != "" exit

    * number of estimated coefficients: columns with a positive variance
    tempname V
    matrix `V' = e(V)
    local k 0
    forvalues j = 1/`= colsof(`V')' {
        if `V'[`j', `j'] > 0 & `V'[`j', `j'] < . local ++k
    }
    if `k' == 0 local k = colsof(e(b))

    tempvar ld
    capture quietly predict double `ld' if $ATEST_esamp, ldisplace
    if _rc {
        local rc = _rc
        _atest_post, name(outliers) type(notrun) ///
            desc("Influential observations (likelihood displacement)") ///
            note("predict, ldisplace not available (rc=`rc')")
        exit
    }
    local cut = $ATEST_ldcut
    if `cut' <= 0 local cut = invchi2(`k', 0.5)
    quietly summarize `ld' if $ATEST_esamp
    local mx = r(max)
    quietly count if `ld' > `cut' & `ld' < . & $ATEST_esamp
    local nout = r(N)
    local cuttxt = strtrim(string(`cut', "%9.3g"))
    if substr("`cuttxt'", 1, 1) == "." local cuttxt "0`cuttxt'"
    local mxtxt  = string(`mx', "%6.3f")
    if `nout' > 0 {
        _atest_post, name(outliers) flag stat(`mx') statname(max LD) ///
            desc("Influential observations: `nout' observation(s) with likelihood displacement > `cuttxt'") ///
            note(`"`note'"')
    }
    else {
        _atest_post, name(outliers) stat(`mx') statname(max LD) ///
            desc("Influential observations: no likelihood displacement > `cuttxt'") ///
            note(`"`note'"')
    }
end
