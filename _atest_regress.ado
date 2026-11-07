*! _atest_regress 0.1.2  assumption tests after regress (linear regression, OLS)
*! Module of assumptiontest, called as  _atest_`e(cmd)'  by the dispatcher.
*! Handles robust/cluster/hc2/hc3 variance (diagnostics on a refit with the
*! default OLS variance: same coefficients and residuals) and weights.
program define _atest_regress
    version 17

    local vce `e(vce)'
    if "`vce'" == "" local vce ols
    local has_robust = ("`vce'" != "ols")
    local wtype `e(wtype)'
    local depvar `e(depvar)'

    display as text "  Model structure: linear regression (OLS), variance: `vce'" _c
    if "`wtype'" != "" display as text ", `wtype's" _c
    display
    display

    * ------------------------------------------------------------------
    * 1. Spline groups (information only) and linearity of the mean
    * ------------------------------------------------------------------
    _atest_splinfo
    _atest_linear, contlist($ATEST_cont) scale(identity)

    * ------------------------------------------------------------------
    * 2. Residual diagnostics. estat hettest/ovtest/vif are not allowed after
    *    robust/cluster variance; the coefficients and residuals do not depend
    *    on the variance option, so they are computed on a default-variance refit.
    * ------------------------------------------------------------------
    local note ""
    local rc 0
    if `has_robust' {
        _atest_refit, remove(vce robust cluster hc2 hc3)
        local rc = r(rc)
        local note "diagnostic refit with default (OLS) variance"
    }
    if `rc' {
        foreach t in homosked reset normal vif outliers {
            _atest_post, name(`t') type(notrun) desc("Residual diagnostic `t'") ///
                note("default-variance refit failed (rc=`rc')")
        }
    }
    else {
        * ---- constant error variance: Breusch-Pagan/Cook-Weisberg, iid form (Koenker)
        capture quietly estat hettest, iid
        if _rc == 0 {
            local c2 = r(chi2)
            local dd = r(df)
            local pp = r(p)
            if `has_robust' {
                _atest_post, name(homosked) type(info) ///
                    desc("Constant error variance (Breusch-Pagan, iid form); not assumed by the `vce' variance") ///
                    stat(`c2') df(`dd') p(`pp') statname(chi2(`dd')) note("`note'")
            }
            else {
                _atest_post, name(homosked) critical ///
                    desc("Constant error variance (Breusch-Pagan/Cook-Weisberg test against fitted values, iid form)") ///
                    stat(`c2') df(`dd') p(`pp') statname(chi2(`dd'))
            }
        }
        else {
            local rch = _rc
            _atest_post, name(homosked) type(notrun) desc("Constant error variance (Breusch-Pagan)") ///
                note("estat hettest failed (rc=`rch'; e.g. not allowed with pweights)")
        }

        * ---- functional form: Ramsey RESET with powers of the fitted values
        capture quietly estat ovtest
        if _rc == 0 {
            local ff = r(F)
            local d1 = r(df)
            local d2 = r(df_r)
            local pp = r(p)
            _atest_post, name(reset) ///
                desc("Functional form / omitted terms (Ramsey RESET, powers of the fitted values)") ///
                stat(`ff') df(`d1') p(`pp') statname(F(`d1',`d2')) note("`note'")
        }
        else {
            local rch = _rc
            _atest_post, name(reset) type(notrun) desc("Functional form (Ramsey RESET)") ///
                note("estat ovtest failed (rc=`rch')")
        }

        * ---- normality of the residuals (skewness/kurtosis test)
        *      (own note: the unadjusted-test remark must not carry over to later rows)
        local snote "`note'"
        tempvar res
        quietly predict double `res' if $ATEST_esamp, residuals
        capture quietly sktest `res'
        local rck = _rc
        local sc = .
        local sp = .
        if `rck' == 0 {
            capture local sc = r(chi2)
            capture local sp = r(P_chi2)
            if missing(`sp') {
                tempname U
                capture matrix `U' = r(Utest)
                if _rc == 0 {
                    local sc = `U'[1, colsof(`U') - 1]
                    local sp = `U'[1, colsof(`U')]
                }
            }
        }
        * the adjusted statistic is missing for extreme kurtosis: unadjusted test instead
        local sn "adj chi2(2)"
        if `rck' == 0 & missing(`sp') {
            capture quietly sktest `res', noadjust
            if _rc == 0 {
                capture local sc = r(chi2)
                capture local sp = r(P_chi2)
                if missing(`sp') {
                    tempname U2
                    capture matrix `U2' = r(Utest)
                    if _rc == 0 {
                        local sc = `U2'[1, colsof(`U2') - 1]
                        local sp = `U2'[1, colsof(`U2')]
                    }
                }
                local sn "chi2(2)"
                local snote = strtrim("`snote' unadjusted test (adjusted statistic not computable)")
            }
        }
        if `rck' == 0 & !missing(`sp') {
            _atest_post, name(normal) ///
                desc("Normally distributed errors (skewness/kurtosis test of the residuals; matters mainly in small samples)") ///
                stat(`sc') df(2) p(`sp') statname(`sn') note("`snote'")
        }
        else {
            _atest_post, name(normal) type(notrun) desc("Normally distributed errors") ///
                note("sktest of the residuals failed or gave no p-value (rc=`rck')")
        }

        * ---- multicollinearity: largest variance inflation factor
        *      (interaction terms and spline() variables are collinear by construction: skipped)
        capture quietly estat vif
        if _rc == 0 {
            local spl $ATEST_spl
            local mx 0
            local mxn ""
            local i 1
            while `"`r(name_`i')'"' != "" {
                local nm `"`r(name_`i')'"'
                local vv = r(vif_`i')
                local skip = (strpos(`"`nm'"', "#") > 0)
                if `: list posof `"`nm'"' in spl' > 0 local skip 1
                if !`skip' & `vv' > `mx' & `vv' < . {
                    local mx = `vv'
                    local mxn `"`nm'"'
                }
                local ++i
            }
            local vtxt = strtrim(string($ATEST_vifcut, "%9.3g"))
            if `mx' > $ATEST_vifcut {
                _atest_post, name(vif) flag stat(`mx') statname(max VIF) ///
                    desc("Multicollinearity: largest variance inflation factor > `vtxt' (`mxn')")
            }
            else {
                _atest_post, name(vif) stat(`mx') statname(max VIF) ///
                    desc("Multicollinearity: no variance inflation factor > `vtxt'")
            }
        }
        else {
            local rch = _rc
            _atest_post, name(vif) type(notrun) desc("Multicollinearity (variance inflation factors)") ///
                note("estat vif failed (rc=`rch')")
        }

        * ---- influential observations (Cook's distance) and outliers (studentized residuals)
        if "$ATEST_nooutliers" == "" {
            tempvar cd rs
            capture quietly predict double `cd' if $ATEST_esamp, cooksd
            local rc1 = _rc
            capture quietly predict double `rs' if $ATEST_esamp, rstudent
            local rc2 = _rc
            local ctxt = strtrim(string($ATEST_cookcut, "%6.3g"))
            if substr("`ctxt'", 1, 1) == "." local ctxt "0`ctxt'"
            if `rc1' == 0 {
                quietly summarize `cd' if $ATEST_esamp
                local mx = r(max)
                quietly count if `cd' > $ATEST_cookcut & `cd' < . & $ATEST_esamp
                local nfl = r(N)
                if `nfl' > 0 {
                    _atest_post, name(outliers) flag stat(`mx') statname(max Cook's D) ///
                        desc("Influential observations: `nfl' observation(s) with Cook's distance > `ctxt'") note("`note'")
                }
                else {
                    _atest_post, name(outliers) stat(`mx') statname(max Cook's D) ///
                        desc("Influential observations: no Cook's distance > `ctxt'") note("`note'")
                }
            }
            else {
                _atest_post, name(outliers) type(notrun) desc("Influential observations (Cook's distance)") ///
                    note("predict, cooksd failed (rc=`rc1')")
            }
            if `rc2' == 0 {
                quietly count if abs(`rs') > $ATEST_rescut & `rs' < . & $ATEST_esamp
                local nres = r(N)
                _atest_post, name(res_out) type(info) ///
                    desc("`nres' observation(s) with |studentized residual| > $ATEST_rescut (information only)")
            }
        }
    }
    _atest_restore

    * ------------------------------------------------------------------
    * 3. Assumptions that cannot be tested on the data
    * ------------------------------------------------------------------
    local indtxt "Independence of the errors (no unmodelled clustering or serial correlation)"
    if inlist("`vce'", "cluster") local indtxt "Errors independent apart from the clustering handled by the cluster variance"
    _atest_post, name(indep) type(untestable) desc("`indtxt'")
    _atest_post, name(cond_mean) type(untestable) ///
        desc("Conditional mean correctly specified: no omitted confounders, covariates exogenous and measured without error")
end
