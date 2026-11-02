*! _atest_splinfo 0.1.0  report spline groups declared with spline()
*! Internal helper of assumptiontest. Linearity is NOT assumed for variables in
*! a spline group, so the joint Wald test of the non-linear terms is informational.
*! Must be called while the original estimation is active.
program define _atest_splinfo
    version 17
    forvalues g = 1/$ATEST_nsp {
        local gv ${ATEST_sp`g'}
        local gbase : word 1 of `gv'
        local ghi : list gv - gbase
        capture quietly test `ghi'
        if _rc == 0 {
            local df = r(df)
            _atest_post, name(spl_`gbase') type(info) ///
                desc("Spline `gbase': joint Wald test of the non-linear terms (linearity is not assumed here)") ///
                stat(`r(chi2)') df(`df') p(`r(p)') statname(chi2(`df'))
        }
        _atest_post, name(splfit_`gbase') type(untestable) ///
            desc("Spline `gbase' adequately captures the functional form (knot number/placement)")
    }
end
