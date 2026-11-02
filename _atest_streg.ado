*! _atest_streg 0.3.0  assumption tests after streg (exponential and Weibull)
*! Module of assumptiontest, called as  _atest_`e(cmd)'  by the dispatcher.
*! Other distributions are flagged and not tested. Handles robust/cluster
*! variance, shared frailty, ancillary()/strata(), PH or AFT metric (time).
program define _atest_streg
    version 17

    local dist `e(dist)'
    if "`dist'" == "" {
        local dist = subinstr(subinstr("`e(cmd)'", "_fr", "", .), "_sw", "", .)
    }
    * with shared()/frailty() e(cmd) carries a "het" suffix (weibullhet, ereghet)
    if substr("`dist'", -3, 3) == "het" local dist = substr("`dist'", 1, strlen("`dist'") - 3)
    * e(cmd) holds Stata's internal names; map them to the dist() names
    if "`dist'" == "ereg"      local dist exponential
    if "`dist'" == "lnormal"   local dist lognormal
    if "`dist'" == "llogistic" local dist loglogistic
    if "`dist'" == "gamma"     local dist ggamma
    local frm  `e(frm2)'
    local cl = lower(`"`e(cmdline)'"')
    local has_robust = (strpos(`"`cl'"', "robust") > 0 | strpos(`"`cl'"', "cluster") > 0)
    local has_shared = (strpos(`"`cl'"', "shared(") > 0 | strpos(`"`cl'"', "frailty(") > 0)
    local has_anc    = (strpos(`"`cl'"', "ancillary(") > 0 | strpos(`"`cl'"', "strata(") > 0)
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

    display as text "  Model structure: dist=`dist', metric=`frm'" _c
    if `has_anc'    display as text ", ancillary()/strata()" _c
    if `has_shared' display as text ", shared frailty" _c
    if `has_robust' display as text ", robust/cluster variance" _c
    if "`multirec'" != "" display as text ", multiple-record data (id: `multirec')" _c
    display
    display

    * ---- only exponential and Weibull are covered
    if !inlist("`dist'", "exponential", "weibull") {
        _atest_post, name(dist_unsupported) flag ///
            desc("Distribution `dist' is not covered by assumptiontest (only exponential and Weibull); no tests were run")
        exit
    }

    local scale = cond("`frm'" == "time", "time", "hazard")

    * ------------------------------------------------------------------
    * 1. Spline groups (information only)
    * ------------------------------------------------------------------
    _atest_splinfo

    * ------------------------------------------------------------------
    * 2. Shape of the baseline hazard
    * ------------------------------------------------------------------
    if "`dist'" == "weibull" {
        _atest_streg_lnp, anc(`has_anc')
        if r(rc) == 0 {
            local dtxt "Weibull shape p = 1 (a constant hazard, i.e. exponential, would suffice if not rejected)"
            if `has_anc' local dtxt "Weibull shape p = 1 in every ancillary()/strata() group (joint Wald test of the ln_p equation)"
            _atest_post, name(shape1) type(info) desc("`dtxt'") ///
                stat(`r(chi2)') df(`r(df)') p(`r(p)') statname(chi2(`r(df)'))
        }
        * generalized gamma nests the Weibull (kappa = 1); AFT metric only
        if `has_anc' {
            _atest_post, name(weib_fit) type(notrun) ///
                desc("Weibull vs generalized gamma") ///
                note("skipped: ancillary()/strata() used")
        }
        else if `has_shared' {
            _atest_post, name(weib_fit) type(notrun) ///
                desc("Weibull vs generalized gamma") ///
                note("skipped: streg does not allow shared()/frailty() with dist(ggamma)")
        }
        else {
            local addo ""
            if "`frm'" == "hazard" local addo "time"
            _atest_refit, remove(dist distribution nohr) addopts(dist(ggamma) `addo')
            local rc = r(rc)
            local rcmd `"`r(cmdline)'"'
            _atest_wald, par(kappa) null(1)
            local rc2 = r(rc)
            if `rc' | `rc2' {
                _atest_post, name(weib_fit) type(notrun) ///
                    desc("Weibull adequacy (nested in the generalized gamma, kappa = 1)") ///
                    note("generalized gamma refit failed (rc=`rc') or kappa test unavailable (rc=`rc2'); refit: `rcmd'")
            }
            else {
                _atest_post, name(weib_fit) critical ///
                    desc("Weibull adequacy (Wald test of kappa = 1 in the generalized gamma)") ///
                    stat(`r(chi2)') df(1) p(`r(p)') statname(chi2(1))
            }
            _atest_restore
        }
    }
    else {
        * exponential: constant hazard tested against the Weibull
        _atest_refit, remove(dist distribution nohr) addopts(dist(weibull))
        local rc = r(rc)
        local rcmd `"`r(cmdline)'"'
        _atest_streg_lnp, anc(`has_anc')
        local rc2 = r(rc)
        if `rc' | `rc2' {
            _atest_post, name(const_haz) type(notrun) ///
                desc("Constant hazard (exponential) vs Weibull") ///
                note("Weibull refit failed (rc=`rc') or shape test unavailable (rc=`rc2'); refit: `rcmd'")
        }
        else {
            local dtxt "Constant hazard over time (Wald test of Weibull shape p = 1)"
            if `has_anc' local dtxt "Constant hazard over time in every strata()/ancillary() group (joint Wald test of Weibull ln_p = 0)"
            _atest_post, name(const_haz) critical desc("`dtxt'") ///
                stat(`r(chi2)') df(`r(df)') p(`r(p)') statname(chi2(`r(df)'))
        }
        _atest_restore
    }

    * ------------------------------------------------------------------
    * 3. Proportional hazards and influence, on the equivalent Cox model
    *    (both distributions are PH models; the Cox-based Schoenfeld test is
    *    distribution free). Weights and all options are dropped in the refit.
    * ------------------------------------------------------------------
    * strata() of streg (dummies in all equations) becomes strata() of stcox
    local cxo ""
    local cxt "on equivalent stcox model"
    if ustrregexm(`"`e(cmdline)'"', "strata\(([^)]*)\)") {
        local stv = strtrim(ustrregexs(1))
        local cxo "strata(`stv')"
        local cxt "on equivalent stcox model, strata(`stv')"
    }
    _atest_refit, newcmd(stcox) remove(all) nowt addopts(`cxo')
    local rc = r(rc)
    if `rc' {
        _atest_post, name(ph_global) type(notrun) ///
            desc("Proportional hazards (Schoenfeld residuals of the equivalent Cox model)") ///
            note("stcox refit failed (rc=`rc')")
        _atest_post, name(outliers) type(notrun) ///
            desc("Influential observations (equivalent Cox model)") ///
            note("stcox refit failed (rc=`rc')")
    }
    else {
        _atest_phtest, note("`cxt'")
        _atest_infl, note("`cxt'")
    }
    _atest_restore

    * ------------------------------------------------------------------
    * 4. Linearity, calibration
    * ------------------------------------------------------------------
    _atest_linear, contlist($ATEST_cont) scale(`scale')

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

    * ------------------------------------------------------------------
    * 5. Assumptions that cannot be tested on the data
    * ------------------------------------------------------------------
    _atest_post, name(cens) type(untestable) ///
        desc("Censoring is non-informative (independent of the event process)")
    local indtxt "Independence of observations"
    if `has_robust' local indtxt "Independence of observations (clustering handled by the robust/cluster variance)"
    _atest_post, name(indep) type(untestable) desc("`indtxt'")
    if `has_shared' {
        _atest_post, name(frailty) type(untestable) ///
            desc("Frailty distribution is correctly specified and independent of the covariates")
    }
    _atest_post, name(omitted) type(untestable) ///
        desc("No important omitted covariates; correct time scale")
end

* ---- Wald test of Weibull shape p = 1 (ln_p = 0) on the active estimates.
*      anc(1): ln_p has covariates (ancillary()/strata()); then all ln_p
*      coefficients INCLUDING the constant are tested jointly (test [ln_p]
*      alone leaves out the constant). Returns r(rc) r(chi2) r(df) r(p).
program define _atest_streg_lnp, rclass
    version 17
    syntax , ANC(integer)
    return scalar rc = 111
    if `anc' {
        capture quietly test [ln_p]
        if _rc == 0 {
            capture quietly test [ln_p]_cons = 0, accum
        }
        else {
            capture quietly test [ln_p]_cons = 0
        }
        if _rc == 0 {
            return scalar rc   = 0
            return scalar chi2 = r(chi2)
            return scalar df   = r(df)
            return scalar p    = r(p)
            exit
        }
    }
    _atest_wald, par(ln_p)
    if r(rc) == 0 {
        local c2 = r(chi2)
        local pp = r(p)
        return scalar rc   = 0
        return scalar chi2 = `c2'
        return scalar df   = 1
        return scalar p    = `pp'
    }
end
