*! _atest_phtest 0.1.0  run estat phtest on the ACTIVE stcox results and post all results
*! Internal helper of assumptiontest. Posts one test per covariate plus the
*! global test (all critical). If graphpath() was given, a Schoenfeld-residual
*! plot is written to disk for covariates failing the test, and dropped from memory.
program define _atest_phtest
    version 17
    syntax [, NOte(string) ]

    capture quietly estat phtest, detail
    if _rc {
        local rc = _rc
        _atest_post, name(ph_global) type(notrun) ///
            desc("Proportional hazards (Schoenfeld residuals)") ///
            note("estat phtest failed (rc=`rc')")
        exit
    }
    local gchi = r(chi2)
    local gdf  = r(df)
    local gp   = r(p)
    tempname P
    capture matrix `P' = r(phtest)
    local havemat = (_rc == 0)
    local failvars ""
    if `havemat' {
        local rn : rownames `P'
        local nr = rowsof(`P')
        forvalues i = 1/`nr' {
            local rname : word `i' of `rn'
            local pp = `P'[`i', 4]
            local cc = `P'[`i', 2]
            local dd = `P'[`i', 3]
        * skip base-level / omitted rows (missing statistics)
        if missing(`pp') | missing(`cc') continue
            local rid = subinstr("`rname'", ".", "_", .)
        _atest_post, name(ph_`rid') critical ///
                desc("Proportional hazards for `rname' (Schoenfeld residuals)") ///
                stat(`cc') df(`dd') p(`pp') statname(chi2(`dd'))
            if `pp' < $ATEST_alpha local failvars `failvars' `rname'
        }
    }
    _atest_post, name(ph_global) critical ///
        desc("Proportional hazards, global test (Schoenfeld residuals)") ///
        stat(`gchi') df(`gdf') p(`gp') statname(chi2(`gdf')) note("`note'")

    * graphs only for failing covariates, saved to disk and dropped at once
    if `"$ATEST_graphpath"' != "" & "`failvars'" != "" {
        foreach v of local failvars {
            * candidate spellings of the covariate: as listed, then without the factor prefix
            local dot = strpos("`v'", ".")
            local bare = cond(`dot' > 0, substr("`v'", `dot' + 1, .), "`v'")
            local done 0
            local why ""
            foreach cand in `v' `bare' {
                if `done' continue
                * the graph must be drawn (not -nodraw-) for graph export to find its window
                capture quietly estat phtest, plot(`cand') name(_atest_g, replace)
                local rcA = _rc
                local rcB = _rc
                if `rcB' == 0 {
                    local fn = subinstr("`v'", ".", "_", .)
                    capture quietly graph export `"$ATEST_graphpath/ph_`fn'.png"', name(_atest_g) replace width(1200)
                    local rcC = _rc
                    if `rcC' == 0 {
                        display as text "      graph saved: $ATEST_graphpath/ph_`fn'.png"
                        display as text "      (a flat smooth = proportional hazards; a clear trend over time = violation)"
                        local done 1
                    }
                    else {
                        local why "`why' [plot(`cand') ok, graph export rc=`rcC']"
                    }
                    capture graph drop _atest_g
                }
                else {
                    local why "`why' [plot(`cand') rc=`rcA'/`rcB']"
                }
            }
            if !`done' {
                display as text "      (graph for `v' not created:`why')"
            }
        }
    }
end
