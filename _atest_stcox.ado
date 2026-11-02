*! _atest_stcox 0.2.0  assumption tests after stcox
*! Module of assumptiontest, called as  _atest_`e(cmd)'  by the dispatcher.
*! Covers: stratified models (strata()), robust/cluster variance, shared
*! frailty (shared()), time-varying effects via tvc()/texp(), and multiple-record
*! data produced by stsplit.
program define _atest_stcox
    version 17

    local cl = lower(`"`e(cmdline)'"')
    local has_tvc    = (strpos(`"`cl'"', "tvc(") > 0)
    local has_robust = (strpos(`"`cl'"', "robust") > 0 | strpos(`"`cl'"', "cluster") > 0)
    local has_shared = (strpos(`"`cl'"', "shared(") > 0)
    local has_strata = (strpos(`"`cl'"', "strata(") > 0)
    local stid : char _dta[st_id]
    local multirec ""
    if "`stid'" != "" {
        tempvar tg
        quietly egen byte `tg' = tag(`stid') if $ATEST_esamp
        quietly count if `tg' == 1
        local nsubj = r(N)
        quietly count if $ATEST_esamp
        if r(N) > `nsubj' local multirec "`stid'"
    }
    local ties = cond("`e(ties)'" == "", "efron", "`e(ties)'")

    display as text "  Model structure: ties=`ties'" _c
    if `has_strata' display as text ", stratified" _c
    if `has_shared' display as text ", shared frailty" _c
    if `has_tvc'    display as text ", tvc()/texp()" _c
    if `has_robust' display as text ", robust/cluster variance" _c
    if "`multirec'" != "" display as text ", multiple-record data (id: `multirec')" _c
    display
    display

    * 1. Spline groups: information only (linearity is not assumed there)
    _atest_splinfo

    * ------------------------------------------------------------------
    * 2. Proportional hazards
    * ------------------------------------------------------------------
    if `has_tvc' {
        capture quietly test [tvc]
        if _rc == 0 {
            _atest_post, name(ph_tvc) critical ///
                desc("Proportional hazards for the tvc() variables (Wald test of the time-varying terms)") ///
                stat(`r(chi2)') df(`r(df)') p(`r(p)') statname(chi2(`r(df)'))
        }
        _atest_post, name(ph_other) type(untestable) ///
            desc("Proportional hazards for covariates NOT listed in tvc() (not tested when tvc() is used)")
    }
    else {
        local note ""
        local rc 0
        if `has_robust' | `has_shared' {
            _atest_refit, remove(vce robust cluster shared frailty)
            local rc = r(rc)
            local note "diagnostic refit without robust/cluster/shared()"
        }
        if `rc' {
            _atest_post, name(ph_global) type(notrun) ///
                desc("Proportional hazards (Schoenfeld residuals)") ///
                note("diagnostic refit failed (rc=`rc')")
        }
        else {
            _atest_phtest, note("`note'")
        }
        _atest_restore
    }

    * ------------------------------------------------------------------
    * 3. Linearity of continuous covariates in the log hazard
    * ------------------------------------------------------------------
    _atest_linear, contlist($ATEST_cont) scale(hazard)

    * ------------------------------------------------------------------
    * 4. Calibration (martingale residuals) and influential observations
    * ------------------------------------------------------------------
    if `has_tvc' {
        * martingale residuals are not available after tvc(); a refit without
        * the time-varying terms would describe a different model, so not run
        _atest_post, name(gb_fit) type(notrun) ///
            desc("Calibration across risk groups (Gronnesby-Borgan type)") ///
            note("predict, mgale not available after tvc(); split the data with stsplit and use time-by-covariate interactions instead")
    }
    else {
        tempvar xb mg
        capture quietly predict double `xb' if $ATEST_esamp, xb
        local rc1 = _rc
        capture quietly predict double `mg' if $ATEST_esamp, mgale
        local rc2 = _rc
        if `rc1' | `rc2' {
            _atest_post, name(gb_fit) type(notrun) ///
                desc("Calibration across risk groups (Gronnesby-Borgan type)") ///
                note("martingale residuals not available for this model")
        }
        else {
            _atest_gb, xb(`xb') obs(_d) mg(`mg')
        }
    }
    if `has_tvc' | `has_robust' | `has_shared' {
        * predict ldisplace/dfbeta are not available after tvc(), robust/cluster
        * or shared(): influence is assessed on a diagnostic refit without them
        local what ""
        if `has_tvc'    local what "`what'tvc()/texp() "
        if `has_robust' local what "`what'robust/cluster "
        if `has_shared' local what "`what'shared() "
        local what = strtrim("`what'")
        _atest_refit, remove(tvc texp vce robust cluster shared frailty)
        local rc = r(rc)
        if `rc' == 0 {
            _atest_infl, note("diagnostic refit without `what'")
        }
        else {
            _atest_post, name(outliers) type(notrun) ///
                desc("Influential observations") note("diagnostic refit without `what' failed (rc=`rc')")
        }
        _atest_restore
    }
    else {
        _atest_infl
    }

    * ------------------------------------------------------------------
    * 5. Assumptions that cannot be tested on the data
    * ------------------------------------------------------------------
    _atest_post, name(cens) type(untestable) ///
        desc("Censoring is non-informative (independent of the event process)")
    if `has_robust' {
        _atest_post, name(indep) type(untestable) ///
            desc("Observations independent apart from the clustering handled by the robust/cluster variance")
    }
    else {
        _atest_post, name(indep) type(untestable) ///
            desc("Independence of observations (no unmodelled clustering)")
    }
    if "`multirec'" != "" {
        _atest_post, name(multirec) type(untestable) ///
            desc("Multiple records per subject (id `multirec') carry no extra information beyond the covariate history")
    }
    if `has_shared' {
        _atest_post, name(frailty) type(untestable) ///
            desc("Frailty distribution is correctly specified and independent of the covariates")
    }
    if `has_strata' {
        _atest_post, name(strata) type(untestable) ///
            desc("Strata capture all baseline-hazard differences; proportionality is assumed within strata only")
    }
    _atest_post, name(omitted) type(untestable) ///
        desc("No important omitted covariates; correct time scale and handling of ties (`ties')")
end
