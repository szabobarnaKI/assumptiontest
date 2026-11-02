*! _atest_linear 0.2.0  automatic linearity tests for plain continuous covariates
*! Internal helper of assumptiontest. Method (global ATEST_lin):
*!   rcs       (default) restricted cubic spline with 3 knots (Stata mkspline, cubic
*!             nknots(3), knots at the 10th/50th/90th percentiles). The model already
*!             contains the linear term, so only the single non-linear spline term is
*!             added; its Wald test (1 df) is the linearity test.
*!   quadratic a centred quadratic term is added; its Wald test is the linearity test.
*! All terms are added in ONE refit of the model (robust/cluster/strata/frailty
*! options are kept, so the Wald tests respect the fitted structure). If a 3-knot
*! spline cannot be built for a variable (too few distinct values), the quadratic
*! term is used for that variable and the output says so.
program define _atest_linear
    version 17
    syntax [, CONTlist(string) SCale(string)]

    local method = cond("$ATEST_lin" == "", "rcs", "$ATEST_lin")
    local addlist ""
    local vlist ""
    local made ""
    local i 0
    foreach v of local contlist {
        capture confirm numeric variable `v'
        if _rc continue
        * dichotomous 0/1 variables cannot be non-linear
        quietly count if $ATEST_esamp & `v' != 0 & `v' != 1 & !missing(`v')
        if r(N) == 0 continue
        quietly summarize `v' if $ATEST_esamp, meanonly
        if r(min) == r(max) continue
        local mean = r(mean)
        local ++i
        local how_`i' quad
        if "`method'" == "rcs" {
            capture quietly mkspline __atrcs`i'_ = `v' if $ATEST_esamp, cubic nknots(3)
            if _rc == 0 {
                capture confirm numeric variable __atrcs`i'_2
                if _rc == 0 {
                    local how_`i' rcs
                    local made `made' __atrcs`i'_1 __atrcs`i'_2
                    local add_`i' __atrcs`i'_2
                }
            }
        }
        if "`how_`i''" == "quad" {
            tempvar q`i'
            quietly generate double `q`i'' = (`v' - `mean')^2 if $ATEST_esamp
            local add_`i' `q`i''
        }
        local addlist `addlist' `add_`i''
        local vlist `vlist' `v'
        local v_`i' `v'
    }

    if "`vlist'" == "" {
        display as text "  (no non-binary plain continuous covariates: linearity not tested)"
        exit
    }

    _atest_refit, extra(`addlist')
    local rc = r(rc)
    local n = `i'
    forvalues j = 1/`n' {
        local v `v_`j''
        local lab = cond("`how_`j''" == "rcs", ///
            "restricted cubic spline, 3 knots", "quadratic term")
        if "`method'" == "rcs" & "`how_`j''" == "quad" {
            local lab "quadratic term; 3-knot spline not possible"
        }
        if `rc' {
            _atest_post, name(lin_`v') type(notrun) ///
                desc("Linearity of the log-`scale' in `v'") ///
                note("refit with added terms failed (rc=`rc')")
            continue
        }
        capture quietly test `add_`j''
        if _rc {
            _atest_post, name(lin_`v') type(notrun) ///
                desc("Linearity of the log-`scale' in `v'") ///
                note("added term omitted by the estimator")
        }
        else {
            local c2 = r(chi2)
            local pp = r(p)
            _atest_post, name(lin_`v') critical ///
                desc("Linearity of the log-`scale' in `v' (`lab'; Wald test of the non-linear term)") ///
                stat(`c2') df(1) p(`pp') statname(chi2(1))
        }
    }
    capture drop `made'
    _atest_restore
end
