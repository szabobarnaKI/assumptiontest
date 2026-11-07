*! _atest_linear 0.3.1  automatic linearity tests for plain continuous covariates
*! Internal helper of assumptiontest. Method (global ATEST_lin):
*!   rcs       (default) restricted cubic spline with 3 knots (Stata mkspline, cubic
*!             nknots(3), knots at the 10th/50th/90th percentiles). The model already
*!             contains the linear term, so only the single non-linear spline term is
*!             added; its test (1 df) is the linearity test.
*!   quadratic a centred quadratic term is added; its test is the linearity test.
*! Default: all terms are added in ONE refit of the model (robust/cluster/strata/
*! frailty options are kept, so the Wald tests respect the fitted structure).
*! Option lr: per variable, a profile (penalized) likelihood-ratio test: full model
*! vs the same design with the added coefficient constrained to 0 (firthlogit,
*! where Wald tests are unreliable under separation); Wald if constraints fail. After regress the Wald test is an F test.
*! scale(): hazard, time, mean, odds (shown as log-...) or identity (shown as mean).
*! If a 3-knot spline cannot be built for a variable (too few distinct values),
*! the quadratic term is used for that variable and the output says so.
program define _atest_linear
    version 17
    syntax [, CONTlist(string) SCale(string) LR]

    local method = cond("$ATEST_lin" == "", "rcs", "$ATEST_lin")
    local sclab = cond("`scale'" == "identity", "mean", "log-`scale'")
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
    local n = `i'

    * ------------------------------------------------------------------
    * likelihood-ratio version: one refit per variable
    * ------------------------------------------------------------------
    if "`lr'" != "" {
        * profile (penalized) likelihood-ratio test (Heinze and Schemper 2002):
        * full model with the added term vs the SAME design with its coefficient
        * constrained to 0. Comparing with the fitted model instead would be wrong
        * for penalized likelihoods, because the Firth penalty depends on the design.
        * If the estimator does not accept constraints(), the Wald test is used.
        forvalues j = 1/`n' {
            local v `v_`j''
            local lab = cond("`how_`j''" == "rcs", "restricted cubic spline, 3 knots", "quadratic term")
            if "`method'" == "rcs" & "`how_`j''" == "quad" local lab "quadratic term; 3-knot spline not possible"
            _atest_refit, extra(`add_`j'')
            local rc = r(rc)
            local ll1 = .
            local se1 = 0
            local wc = .
            local wp = .
            if `rc' == 0 {
                local ll1 = e(ll)
                capture local se1 = _se[`add_`j'']
                if _rc local se1 = 0
                capture quietly test `add_`j''
                if _rc == 0 {
                    local wc = r(chi2)
                    local wp = r(p)
                }
            }
            _atest_restore
            if `rc' | `se1' == 0 | missing(`se1') {
                local why = cond(`rc', "refit with added term failed (rc=`rc')", "added term omitted")
                _atest_post, name(lin_`v') type(notrun) ///
                    desc("Linearity of the `sclab' in `v'") note("`why'")
                continue
            }
            local ll0 = .
            local rc0 111
            quietly constraint free
            local cn = r(free)
            quietly constraint define `cn' `add_`j'' = 0
            _atest_refit, extra(`add_`j'') addopts(constraints(`cn'))
            local rc0 = r(rc)
            if `rc0' == 0 local ll0 = e(ll)
            _atest_restore
            constraint drop `cn'
            if `rc0' == 0 & !missing(`ll0') & !missing(`ll1') {
                local c2 = max(0, 2 * (`ll1' - `ll0'))
                local pp = chi2tail(1, `c2')
                _atest_post, name(lin_`v') critical ///
                    desc("Linearity of the `sclab' in `v' (`lab'; profile likelihood-ratio test of the non-linear term)") ///
                    stat(`c2') df(1) p(`pp') statname(LR chi2(1))
            }
            else if !missing(`wp') {
                _atest_post, name(lin_`v') critical ///
                    desc("Linearity of the `sclab' in `v' (`lab'; Wald test of the non-linear term)") ///
                    stat(`wc') df(1) p(`wp') statname(chi2(1)) ///
                    note("constrained refit failed (rc=`rc0'): Wald test used")
            }
            else {
                _atest_post, name(lin_`v') type(notrun) ///
                    desc("Linearity of the `sclab' in `v'") note("no likelihood-ratio or Wald test available")
            }
        }
        capture drop `made'
        exit
    }

    * ------------------------------------------------------------------
    * Wald version: all terms in one refit
    * ------------------------------------------------------------------
    _atest_refit, extra(`addlist')
    local rc = r(rc)
    forvalues j = 1/`n' {
        local v `v_`j''
        local lab = cond("`how_`j''" == "rcs", ///
            "restricted cubic spline, 3 knots", "quadratic term")
        if "`method'" == "rcs" & "`how_`j''" == "quad" {
            local lab "quadratic term; 3-knot spline not possible"
        }
        if `rc' {
            _atest_post, name(lin_`v') type(notrun) ///
                desc("Linearity of the `sclab' in `v'") ///
                note("refit with added terms failed (rc=`rc')")
            continue
        }
        capture quietly test `add_`j''
        if _rc {
            _atest_post, name(lin_`v') type(notrun) ///
                desc("Linearity of the `sclab' in `v'") ///
                note("added term omitted by the estimator")
        }
        else {
            if r(chi2) < . {
                local c2 = r(chi2)
                local sn "chi2(1)"
            }
            else {
                local c2 = r(F)
                local sn "F(1,`r(df_r)')"
            }
            local pp = r(p)
            _atest_post, name(lin_`v') critical ///
                desc("Linearity of the `sclab' in `v' (`lab'; Wald test of the non-linear term)") ///
                stat(`c2') df(1) p(`pp') statname(`sn')
        }
    }
    capture drop `made'
    _atest_restore
end
