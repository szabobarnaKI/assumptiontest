*! _atest_gb 0.1.0  Gronnesby-Borgan-type calibration test from martingale residuals
*! Internal helper of assumptiontest. Observations are grouped by quantiles of
*! the linear predictor; observed events O_k are compared with expected events
*! E_k = O_k - sum(martingale residual). Statistic: sum (O_k-E_k)^2 / O_k,
*! referred to chi2(G-1). Approximate, intended as an automatic screening test.
program define _atest_gb
    version 17
    syntax , XB(varname) OBS(varname) MG(varname) [NQ(integer 10)]

    tempvar grp
    quietly xtile `grp' = `xb' if $ATEST_esamp, nq(`nq')
    quietly levelsof `grp', local(lv)
    local G 0
    local stat 0
    foreach k of local lv {
        quietly summarize `obs' if `grp' == `k', meanonly
        local O = r(sum)
        quietly summarize `mg' if `grp' == `k', meanonly
        local E = `O' - r(sum)
        if `O' >= 1 {
            local stat = `stat' + (`O' - `E')^2 / `O'
            local ++G
        }
    }
    if `G' < 3 {
        _atest_post, name(gb_fit) type(notrun) ///
            desc("Calibration across risk groups (Gronnesby-Borgan type)") ///
            note("fewer than 3 risk groups with events")
        exit
    }
    local df = `G' - 1
    local p  = chi2tail(`df', `stat')
    _atest_post, name(gb_fit) ///
        desc("Calibration of observed vs expected events across `G' risk groups (Gronnesby-Borgan type)") ///
        stat(`stat') df(`df') p(`p') statname(chi2(`df'))
end
