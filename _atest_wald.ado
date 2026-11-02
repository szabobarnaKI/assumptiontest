*! _atest_wald 0.1.0  Wald test of an ancillary parameter, trying several spellings
*! Internal helper of assumptiontest. Usage: _atest_wald, par(ln_p) [null(0)]
*! Tries [par]_cons, _b[/par], [/par] and _b[par:_cons] against the active estimates.
*! Returns r(rc) (0 = a test was obtained), r(chi2), r(p), r(spec).
program define _atest_wald, rclass
    version 17
    syntax , PAR(string) [NULL(real 0)]

    local cand1 "[`par']_cons = `null'"
    local cand2 "_b[/`par'] = `null'"
    local cand3 "[/`par'] = `null'"
    local cand4 "_b[`par':_cons] = `null'"
    return scalar rc = 111
    forvalues j = 1/4 {
        capture quietly test `cand`j''
        if _rc == 0 {
            return scalar rc   = 0
            return scalar chi2 = r(chi2)
            return scalar p    = r(p)
            return local  spec `"`cand`j''"'
            exit
        }
    }
end
