*! _atest_infl 0.3.0  influence screening for Cox models: dispatch on influence()
*! Internal helper of assumptiontest. influence(ld) (default) uses likelihood
*! displacement, influence(dfbeta) the standardised dfbeta. note() is passed on
*! and printed with the result (e.g. when a diagnostic refit is used).
program define _atest_infl
    version 17
    syntax [, NOte(string)]
    if "$ATEST_nooutliers" != "" exit
    if "$ATEST_infl" == "dfbeta" {
        _atest_dfbeta, outcut($ATEST_outcut) note(`"`note'"')
    }
    else {
        _atest_ld, note(`"`note'"')
    }
end
