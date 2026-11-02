*! _atest_nbreg 0.3.2  assumption tests after nbreg (default NB2 form)
*! Module of assumptiontest, called as  _atest_`e(cmd)'  by the dispatcher.
*! nbreg with dispersion(constant) (NB1) is flagged and not tested.
program define _atest_nbreg
    version 17

    local cl = lower(`"`e(cmdline)'"')
    local has_robust = (strpos(`"`cl'"', "robust") > 0 | strpos(`"`cl'"', "cluster") > 0)
    local depvar `e(depvar)'
    local alpha  = e(alpha)

    display as text "  Model structure: negative binomial (NB2), alpha = " ///
        as result %6.4f `alpha' _c
    if `has_robust' display as text ", robust/cluster variance" _c
    if "`e(offset)'" != "" display as text ", offset/exposure: `e(offset)'" _c
    display
    display

    if "`e(dispers)'" != "mean" {
        _atest_post, name(dist_unsupported) flag ///
            desc("nbreg with dispersion(`e(dispers)') is not covered (only the default NB2 form); no tests were run")
        exit
    }

    * quantities used by several tests
    tempvar mu V
    quietly predict double `mu' if $ATEST_esamp, n
    quietly generate double `V' = `mu' + `alpha' * `mu'^2 if $ATEST_esamp
    quietly count if $ATEST_esamp
    local N = r(N)
    local kpar = e(k)

    * ------------------------------------------------------------------
    * 1. Spline groups (information only)
    * ------------------------------------------------------------------
    _atest_splinfo

    * ------------------------------------------------------------------
    * 2. Over-dispersion: negative binomial vs Poisson (LR test of alpha = 0)
    * ------------------------------------------------------------------
    local lr = .
    local odnote ""
    if e(chi2_c) < . {
        local lr = e(chi2_c)
    }
    else if e(ll_c) < . {
        local lr = 2 * (e(ll) - e(ll_c))
    }
    else if `has_robust' {
        * nbreg does not store the comparison test with robust/cluster variance;
        * the point estimates (and alpha) do not depend on the variance option,
        * so the LR test is taken from a refit with the default variance
        _atest_refit, remove(vce robust cluster nolrtest)
        if r(rc) == 0 {
            if e(chi2_c) < . {
                local lr = e(chi2_c)
            }
            else if e(ll_c) < . {
                local lr = 2 * (e(ll) - e(ll_c))
            }
            local odnote "diagnostic refit without robust/cluster variance"
        }
        _atest_restore
    }
    if `lr' < . {
        local p = chi2tail(1, max(`lr', 0)) / 2
        if `p' < $ATEST_alpha {
            _atest_post, name(overdisp) type(info) ///
                desc("Over-dispersion present: Poisson is rejected in favour of negative binomial (LR test of alpha = 0)") ///
                stat(`lr') df(1) p(`p') statname(chibar2(01)) note("`odnote'")
        }
        else {
            _atest_post, name(overdisp) flag ///
                desc("Over-dispersion not demonstrated: Poisson regression might suffice (LR test of alpha = 0)") ///
                stat(`lr') df(1) p(`p') statname(chibar2(01)) note("`odnote'")
        }
    }
    else {
        _atest_post, name(overdisp) type(notrun) ///
            desc("Over-dispersion (LR test of alpha = 0)") ///
            note("comparison likelihood not available (pweights, or the refit failed)")
    }

    * ------------------------------------------------------------------
    * 3. Linearity of the log mean in continuous covariates
    * ------------------------------------------------------------------
    _atest_linear, contlist($ATEST_cont) scale(mean)

    * ------------------------------------------------------------------
    * 4. Variance function Var = mu + alpha*mu^2 (score-type test)
    *    regress ((y-mu)^2 - V)/V on alpha*mu^2*ln(mu)/V, robust SE
    * ------------------------------------------------------------------
    tempvar u z
    quietly generate double `u' = ((`depvar' - `mu')^2 - `V') / `V' if $ATEST_esamp
    quietly generate double `z' = `alpha' * `mu'^2 * ln(`mu') / `V' if $ATEST_esamp
    capture quietly regress `u' `z' if $ATEST_esamp, noconstant vce(robust)
    local rc = _rc
    if `rc' == 0 {
        capture quietly test `z'
        if _rc == 0 {
            local c2 = r(F)
            local pv = r(p)
            _atest_post, name(varfun) ///
                desc("Variance function Var = mu + alpha*mu^2 (approximate score-type test)") ///
                stat(`c2') df(1) p(`pv') statname(F(1))
        }
        else {
            local rc 1
        }
    }
    if `rc' {
        _atest_post, name(varfun) type(notrun) ///
            desc("Variance function Var = mu + alpha*mu^2") note("auxiliary regression failed")
    }
    _atest_restore

    * ------------------------------------------------------------------
    * 5. Goodness of fit (Pearson chi-squared) and excess zeros
    * ------------------------------------------------------------------
    tempvar pr
    quietly generate double `pr' = (`depvar' - `mu')^2 / `V' if $ATEST_esamp
    quietly summarize `pr' if $ATEST_esamp, meanonly
    local X2 = r(sum)
    local dfp = `N' - `kpar'
    if `dfp' > 0 {
        local pp = chi2tail(`dfp', `X2')
        _atest_post, name(pearson) ///
            desc("Goodness of fit: Pearson chi-squared (approximate, unreliable for very small means)") ///
            stat(`X2') df(`dfp') p(`pp') statname(chi2(`dfp'))
    }

    tempvar p0 z0
    capture quietly predict double `p0' if $ATEST_esamp, pr(0)
    if _rc == 0 {
        quietly generate byte `z0' = (`depvar' == 0) if $ATEST_esamp
        quietly summarize `z0' if $ATEST_esamp, meanonly
        local obs0 = r(sum)
        quietly summarize `p0' if $ATEST_esamp, meanonly
        local exp0 = r(sum)
        tempvar vv
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
    }
    else {
        _atest_post, name(zeros) type(notrun) ///
            desc("Excess zeros") note("predicted zero probabilities not available")
    }

    * ------------------------------------------------------------------
    * 6. Outliers: influential observations, and large Pearson residuals
    * ------------------------------------------------------------------
    if "$ATEST_nooutliers" == "" {
        * equivalent GLM (negative binomial, log link, alpha by ML): Cook's distance
        * and standardized deviance residuals (Pregibon 1981; McCullagh and Nelder 1989)
        _atest_refit, newcmd(glm) remove(dispersion irr nolrtest vce robust cluster) ///
            addopts(family(nbinomial ml) link(log))
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
                local mxtxt = string(`mx', "%6.3f")
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
                        desc("Influential observations: no Cook's distance > `ctxt' (largest `mxtxt'; GLM refit)")
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
    * 7. Assumptions that cannot be tested on the data
    * ------------------------------------------------------------------
    _atest_post, name(cond_mean) type(untestable) ///
        desc("Conditional mean correctly specified: no omitted covariates or interactions beyond those tested")
    _atest_post, name(nb2) type(untestable) ///
        desc("Unobserved heterogeneity is gamma-distributed with mean 1 (NB2 mixture)")
    local indtxt "Independence of observations"
    if `has_robust' local indtxt "Independence of observations (clustering handled by the robust/cluster variance)"
    _atest_post, name(indep) type(untestable) desc("`indtxt'")
    if "`e(offset)'" != "" | "`e(exposure)'" != "" {
        _atest_post, name(exposure) type(untestable) ///
            desc("Counts are proportional to the exposure/offset (coefficient fixed at 1)")
    }
end
