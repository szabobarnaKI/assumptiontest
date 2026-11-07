*! _atest_poisson 0.1.1  assumption tests after poisson
*! Module of assumptiontest, called as  _atest_`e(cmd)'  by the dispatcher.
*! With robust/cluster variance (Poisson quasi-likelihood) the equidispersion
*! assumption is not needed for valid inference and is reported as information.
program define _atest_poisson
    version 17

    local vce `e(vce)'
    local has_robust = inlist("`vce'", "robust", "cluster")
    local depvar `e(depvar)'

    display as text "  Model structure: Poisson regression" _c
    if `has_robust' display as text ", robust/cluster variance" _c
    if "`e(offset)'" != "" display as text ", offset/exposure: `e(offset)'" _c
    if "`e(wtype)'" != "" display as text ", `e(wtype)'s" _c
    display
    display

    * quantities used by several tests
    tempvar mu
    quietly predict double `mu' if $ATEST_esamp, n
    quietly count if $ATEST_esamp
    local N = r(N)
    local kpar = e(k)

    * ------------------------------------------------------------------
    * 1. Spline groups (information only)
    * ------------------------------------------------------------------
    _atest_splinfo

    * ------------------------------------------------------------------
    * 2. Equidispersion Var(y) = mu: Cameron-Trivedi auxiliary regression
    *    ((y-mu)^2 - y)/mu on mu without constant (alternative Var = mu + a*mu^2);
    *    robust SE, F test (two-sided: over- and under-dispersion)
    * ------------------------------------------------------------------
    tempvar u
    quietly generate double `u' = ((`depvar' - `mu')^2 - `depvar') / `mu' if $ATEST_esamp
    capture quietly regress `u' `mu' if $ATEST_esamp, noconstant vce(robust)
    local rc = _rc
    if `rc' == 0 {
        local a = _b[`mu']
        local dfr = e(df_r)
        capture quietly test `mu'
        local rc = _rc
    }
    if `rc' == 0 {
        local ff = r(F)
        local pp = r(p)
        local dir = cond(`a' > 0, "over-dispersion", "under-dispersion")
        local atxt = strtrim(string(`a', "%9.3g"))
        if substr("`atxt'", 1, 1) == "." local atxt "0`atxt'"
        if substr("`atxt'", 1, 2) == "-." local atxt = "-0" + substr("`atxt'", 2, .)
        if `has_robust' {
            _atest_post, name(overdisp) type(info) ///
                desc("Equidispersion Var(y) = mu (Cameron-Trivedi; estimated `dir', a = `atxt'); not assumed by the robust variance") ///
                stat(`ff') df(1) p(`pp') statname(F(1,`dfr'))
        }
        else {
            _atest_post, name(overdisp) critical ///
                desc("Equidispersion Var(y) = mu (Cameron-Trivedi test; estimated `dir', a = `atxt')") ///
                stat(`ff') df(1) p(`pp') statname(F(1,`dfr')) ///
                note("if rejected: vce(robust) or nbreg")
        }
    }
    else {
        _atest_post, name(overdisp) type(notrun) desc("Equidispersion Var(y) = mu") ///
            note("auxiliary regression failed (rc=`rc')")
    }
    _atest_restore

    * ------------------------------------------------------------------
    * 3. Linearity of the log mean in continuous covariates
    * ------------------------------------------------------------------
    _atest_linear, contlist($ATEST_cont) scale(mean)

    * ------------------------------------------------------------------
    * 4. Goodness of fit (Pearson chi-squared) and excess zeros
    * ------------------------------------------------------------------
    tempvar pr
    quietly generate double `pr' = (`depvar' - `mu')^2 / `mu' if $ATEST_esamp
    quietly summarize `pr' if $ATEST_esamp, meanonly
    local X2 = r(sum)
    local dfp = `N' - `kpar'
    if `dfp' > 0 {
        local pp = chi2tail(`dfp', `X2')
        _atest_post, name(pearson) ///
            desc("Goodness of fit: Pearson chi-squared (approximate, unreliable for very small means)") ///
            stat(`X2') df(`dfp') p(`pp') statname(chi2(`dfp'))
    }

    tempvar p0 z0 vv
    quietly generate double `p0' = exp(-`mu') if $ATEST_esamp
    quietly generate byte `z0' = (`depvar' == 0) if $ATEST_esamp
    quietly summarize `z0' if $ATEST_esamp, meanonly
    local obs0 = r(sum)
    quietly summarize `p0' if $ATEST_esamp, meanonly
    local exp0 = r(sum)
    quietly generate double `vv' = `p0' * (1 - `p0') if $ATEST_esamp
    quietly summarize `vv' if $ATEST_esamp, meanonly
    local var0 = r(sum)
    if `var0' > 0 {
        local zz = (`obs0' - `exp0') / sqrt(`var0')
        local pz = 1 - normal(`zz')
        local exptxt = string(`exp0', "%9.1f")
        _atest_post, name(zeros) ///
            desc("No excess zeros (observed `obs0' vs expected `exptxt' zero counts; one-sided)") ///
            stat(`zz') p(`pz') statname(z)
    }
    else {
        _atest_post, name(zeros) type(notrun) desc("Excess zeros") ///
            note("expected zero probabilities are all 0 or 1")
    }

    * ------------------------------------------------------------------
    * 5. Outliers: Cook's distance and standardized deviance residuals from
    *    the equivalent GLM (Poisson, log link), refitted without robust/cluster
    * ------------------------------------------------------------------
    if "$ATEST_nooutliers" == "" {
        _atest_refit, newcmd(glm) remove(irr vce robust cluster) ///
            addopts(family(poisson) link(log))
        local rc = r(rc)
        if `rc' == 0 {
            tempvar cd rd
            capture quietly predict double `cd' if $ATEST_esamp, cooksd
            local rc1 = _rc
            capture quietly predict double `rd' if $ATEST_esamp, deviance standardized
            local rc2 = _rc
            if `rc1' == 0 {
                quietly summarize `cd' if $ATEST_esamp
                local mx = r(max)
                local ctxt  = strtrim(string($ATEST_cookcut, "%6.3g"))
                if substr("`ctxt'", 1, 1) == "." local ctxt "0`ctxt'"
                quietly count if `cd' > $ATEST_cookcut & `cd' < . & $ATEST_esamp
                local nfl = r(N)
                if `nfl' > 0 {
                    _atest_post, name(outliers) flag stat(`mx') statname(max Cook's D) ///
                        desc("Influential observations: `nfl' observation(s) with Cook's distance > `ctxt' (GLM refit)")
                }
                else {
                    _atest_post, name(outliers) stat(`mx') statname(max Cook's D) ///
                        desc("Influential observations: no Cook's distance > `ctxt' (GLM refit)")
                }
            }
            else {
                _atest_post, name(outliers) type(notrun) ///
                    desc("Influential observations (Cook's distance)") note("predict, cooksd failed after glm")
            }
            if `rc2' == 0 {
                quietly count if abs(`rd') > $ATEST_rescut & `rd' < . & $ATEST_esamp
                local nres = r(N)
                _atest_post, name(dev_out) type(info) ///
                    desc("`nres' observation(s) with |standardized deviance residual| > $ATEST_rescut (information only)")
            }
        }
        else {
            _atest_post, name(outliers) type(notrun) ///
                desc("Influential observations (Cook's distance)") note("equivalent glm refit failed (rc=`rc')")
        }
        _atest_restore
    }

    * ------------------------------------------------------------------
    * 6. Assumptions that cannot be tested on the data
    * ------------------------------------------------------------------
    _atest_post, name(cond_mean) type(untestable) ///
        desc("Conditional mean correctly specified: no omitted covariates or interactions beyond those tested")
    local indtxt "Independence of observations"
    if `has_robust' local indtxt "Independence of observations (clustering handled by the robust/cluster variance)"
    _atest_post, name(indep) type(untestable) desc("`indtxt'")
    if "`e(offset)'" != "" {
        _atest_post, name(exposure) type(untestable) ///
            desc("Counts are proportional to the exposure/offset (coefficient fixed at 1)")
    }
end
