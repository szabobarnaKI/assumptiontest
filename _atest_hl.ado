*! _atest_hl 0.1.0  Hosmer-Lemeshow goodness-of-fit test from predicted probabilities
*! Internal helper of assumptiontest. Observations are grouped by quantiles of
*! the predicted probability (as estat gof, group(10)); statistic
*! sum (O_g - E_g)^2 / (E_g (1 - E_g/n_g)), referred to chi2(G-2).
*! Works after any binary model from its predicted probabilities (also after
*! firthlogit, where estat gof is not available). Weights are not used.
program define _atest_hl
    version 17
    syntax , PR(varname) Y(varname) [NQ(integer 10) NOte(string)]

    tempvar grp
    quietly xtile `grp' = `pr' if $ATEST_esamp, nq(`nq')
    quietly levelsof `grp' if $ATEST_esamp, local(lv)
    local G 0
    local stat 0
    foreach k of local lv {
        quietly count if `grp' == `k'
        local ng = r(N)
        quietly summarize `y' if `grp' == `k', meanonly
        local O = r(sum)
        quietly summarize `pr' if `grp' == `k', meanonly
        local E = r(sum)
        local den = `E' * (1 - `E' / `ng')
        if `den' > 0 & `den' < . {
            local stat = `stat' + (`O' - `E')^2 / `den'
            local ++G
        }
    }
    if `G' < 3 {
        _atest_post, name(hl_gof) type(notrun) ///
            desc("Calibration: Hosmer-Lemeshow goodness of fit") ///
            note("fewer than 3 groups of predicted risk")
        exit
    }
    local df = `G' - 2
    local p  = chi2tail(`df', `stat')
    _atest_post, name(hl_gof) ///
        desc("Calibration: Hosmer-Lemeshow goodness of fit (`G' groups of predicted risk)") ///
        stat(`stat') df(`df') p(`p') statname(chi2(`df')) note(`"`note'"')
end
