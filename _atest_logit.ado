*! _atest_logit 0.1.2  assumption tests after logit, logistic and firthlogit
*! Module of assumptiontest, called as  _atest_`e(cmd)'  by the dispatcher
*! (_atest_logistic and _atest_firthlogit call this module).
*! firthlogit (SSC, Coveney): penalized likelihood (Firth 1993). Separation is
*! not a problem there, Wald tests are unreliable under separation, so the
*! linearity tests are penalized likelihood-ratio tests; influence is computed
*! from the penalized fit (observation-level Pregibon delta-beta).
program define _atest_logit
    version 17

    local cmd `e(cmd)'
    local firth = ("`cmd'" == "firthlogit")
    local vce `e(vce)'
    local has_robust = inlist("`vce'", "robust", "cluster")
    local depvar `e(depvar)'

    * ---- outcome, events, number of estimated coefficients
    tempvar y01
    quietly generate byte `y01' = (`depvar' != 0) if $ATEST_esamp & !missing(`depvar')
    quietly count if $ATEST_esamp
    local N = r(N)
    quietly count if `y01' == 1 & $ATEST_esamp
    local nev = r(N)
    local nnon = `N' - `nev'
    tempname V
    matrix `V' = e(V)
    local cn : colnames `V'
    local kx 0
    forvalues j = 1/`= colsof(`V')' {
        local t : word `j' of `cn'
        if "`t'" != "_cons" & `V'[`j', `j'] > 0 & `V'[`j', `j'] < . local ++kx
    }

    local mtxt "logistic regression (`cmd')"
    if `firth' local mtxt "logistic regression, Firth penalized likelihood (firthlogit)"
    display as text "  Model structure: `mtxt', events " as result `nev' ///
        as text " / non-events " as result `nnon' _c
    if `has_robust' display as text ", robust/cluster variance" _c
    if "`e(wtype)'" != "" display as text ", `e(wtype)'s" _c
    display
    display

    * ---- predicted probabilities (firthlogit predicts the linear predictor only)
    tempvar xb p
    local havep 0
    capture quietly predict double `xb' if $ATEST_esamp, xb
    if _rc capture quietly predict double `xb' if $ATEST_esamp
    if _rc {
        local rcp = _rc
        _atest_post, name(predict) type(notrun) desc("Linear predictor") ///
            note("predict failed (rc=`rcp'); calibration and influence not run")
    }
    else {
        quietly generate double `p' = invlogit(`xb') if $ATEST_esamp
        local havep 1
    }

    * ------------------------------------------------------------------
    * 1. Spline groups (information only)
    * ------------------------------------------------------------------
    _atest_splinfo

    * ------------------------------------------------------------------
    * 2. Separation and sparse data
    * ------------------------------------------------------------------
    if `firth' {
        _atest_post, name(separation) type(info) ///
            desc("Separation: finite estimates under (quasi-)complete separation are provided by the Firth penalty")
    }
    else {
        local ncds = cond(e(N_cds) < ., e(N_cds), 0)
        local ncdf = cond(e(N_cdf) < ., e(N_cdf), 0)
        * variables dropped as perfect predictors (rules recorded by logit)
        local rulevars ""
        tempname RU
        capture matrix `RU' = e(rules)
        if _rc == 0 {
            local rn : rownames `RU'
            forvalues i = 1/`= rowsof(`RU')' {
                local any 0
                forvalues j = 1/`= colsof(`RU')' {
                    if `RU'[`i', `j'] != 0 & `RU'[`i', `j'] < . local any 1
                }
                if `any' local rulevars `rulevars' `: word `i' of `rn''
            }
        }
        if `ncds' + `ncdf' > 0 | "`rulevars'" != "" {
            local stxt ""
            if `ncds' + `ncdf' > 0 local stxt "`ncds' success(es) and `ncdf' failure(s) completely determined"
            if "`rulevars'" != "" {
                if "`stxt'" != "" local stxt "`stxt'; "
                local stxt "`stxt'perfect predictor(s) dropped: `rulevars'"
            }
            _atest_post, name(separation) critical flag ///
                desc("No (quasi-)complete separation: `stxt'") ///
                note("maximum-likelihood estimates may not exist; consider firthlogit")
        }
        else {
            _atest_post, name(separation) critical ///
                desc("No (quasi-)complete separation (no completely determined outcomes, no perfect predictors)")
        }
    }
    if `kx' > 0 {
        local epv = min(`nev', `nnon') / `kx'
        local etxt = strtrim(string($ATEST_epvcut, "%9.3g"))
        local epvtxt = strtrim(string(`epv', "%9.2f"))
        if `epv' < $ATEST_epvcut {
            _atest_post, name(epv) flag stat(`epv') statname(EPV) ///
                desc("Sparse data: `epvtxt' events per estimated coefficient (min(events, non-events)/`kx') < `etxt'")
        }
        else {
            _atest_post, name(epv) stat(`epv') statname(EPV) ///
                desc("Events per estimated coefficient (min(events, non-events)/`kx') >= `etxt'")
        }
    }

    * ------------------------------------------------------------------
    * 3. Linearity of the log-odds in continuous covariates
    * ------------------------------------------------------------------
    if `firth' {
        _atest_linear, contlist($ATEST_cont) scale(odds) lr
    }
    else {
        _atest_linear, contlist($ATEST_cont) scale(odds)
    }

    * ------------------------------------------------------------------
    * 4. Calibration (Hosmer-Lemeshow) and link (linktest)
    * ------------------------------------------------------------------
    if `havep' _atest_hl, pr(`p') y(`y01')
    * link test (Pregibon 1980), computed directly: -linktest- drops its _hat
    * and _hatsq variables when it finishes, so they cannot be tested afterwards
    local rcl 111
    if `havep' {
        tempvar hsq
        quietly generate double `hsq' = `xb'^2 if $ATEST_esamp
        local wgt ""
        if "`e(wtype)'" != "" local wgt "[`e(wtype)' `e(wexp)']"
        capture quietly `cmd' `depvar' `xb' `hsq' `wgt' if $ATEST_esamp
        local rcl = _rc
        if `rcl' == 0 {
            capture quietly test `hsq'
            local rcl = _rc
            if `rcl' == 0 {
                local c2 = r(chi2)
                local pp = r(p)
                _atest_post, name(linktest) ///
                    desc("Link / specification (link test: Wald test of the squared linear predictor)") ///
                    stat(`c2') df(1) p(`pp') statname(chi2(1))
            }
        }
    }
    if `rcl' {
        _atest_post, name(linktest) type(notrun) desc("Link / specification (link test)") ///
            note("refit on the linear predictor failed after `cmd' (rc=`rcl')")
    }
    _atest_restore

    * ------------------------------------------------------------------
    * 5. Influential observations (Pregibon delta-beta) and large residuals
    * ------------------------------------------------------------------
    if "$ATEST_nooutliers" == "" & `havep' {
        local ctxt = strtrim(string($ATEST_cookcut, "%6.3g"))
        if substr("`ctxt'", 1, 1) == "." local ctxt "0`ctxt'"
        tempvar db rs
        local note ""
        local rc 0
        if `firth' {
            * observation-level delta-beta from the penalized fit:
            * h = p(1-p) x'Vx (Vx from stdp), r = Pearson residual, db = r^2 h/(1-h)^2
            tempvar sd h r
            capture quietly predict double `sd' if $ATEST_esamp, stdp
            local rc = _rc
            if `rc' == 0 {
                quietly generate double `h' = `p' * (1 - `p') * `sd'^2 if $ATEST_esamp
                quietly generate double `r' = (`y01' - `p') / sqrt(`p' * (1 - `p')) if $ATEST_esamp
                quietly generate double `db' = `r'^2 * `h' / (1 - `h')^2 if $ATEST_esamp & `h' < 1
                quietly generate double `rs' = `r' / sqrt(1 - `h') if $ATEST_esamp & `h' < 1
                local note "observation-level, penalized fit"
            }
            else {
                local note "predict, stdp failed after firthlogit (rc=`rc')"
            }
        }
        else {
            if `has_robust' {
                _atest_refit, remove(vce robust cluster)
                local rc = r(rc)
                local note "diagnostic refit without robust/cluster variance"
            }
            if `rc' == 0 {
                capture quietly predict double `db' if $ATEST_esamp, dbeta
                local rc = _rc
                if `rc' == 0 capture quietly predict double `rs' if $ATEST_esamp, rstandard
            }
            else {
                local note "default-variance refit failed (rc=`rc')"
            }
            _atest_restore
        }
        capture confirm variable `db'
        if _rc == 0 {
            quietly summarize `db' if $ATEST_esamp
            local mx = r(max)
            quietly count if `db' > $ATEST_cookcut & `db' < . & $ATEST_esamp
            local nfl = r(N)
            if `nfl' > 0 {
                _atest_post, name(outliers) flag stat(`mx') statname(max delta-beta) ///
                    desc("Influential observations: `nfl' observation(s) with Pregibon delta-beta > `ctxt'") note("`note'")
            }
            else {
                _atest_post, name(outliers) stat(`mx') statname(max delta-beta) ///
                    desc("Influential observations: no Pregibon delta-beta > `ctxt'") note("`note'")
            }
        }
        else {
            _atest_post, name(outliers) type(notrun) desc("Influential observations (Pregibon delta-beta)") ///
                note("`note'")
        }
        capture confirm variable `rs'
        if _rc == 0 {
            quietly count if abs(`rs') > $ATEST_rescut & `rs' < . & $ATEST_esamp
            local nres = r(N)
            _atest_post, name(res_out) type(info) ///
                desc("`nres' observation(s) with |standardized Pearson residual| > $ATEST_rescut (information only)")
        }
    }

    * ------------------------------------------------------------------
    * 6. Assumptions that cannot be tested on the data
    * ------------------------------------------------------------------
    local indtxt "Independence of observations"
    if `has_robust' local indtxt "Independence of observations (clustering handled by the robust/cluster variance)"
    _atest_post, name(indep) type(untestable) desc("`indtxt'")
    _atest_post, name(cond_mean) type(untestable) ///
        desc("Log-odds correctly specified: no important omitted covariates or interactions beyond those tested")
end
