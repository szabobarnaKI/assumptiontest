*! assumptiontest 0.3.0  model-assumption tests after stcox, streg and nbreg
*! Checks the model type of the last estimation, runs the relevant automatic
*! assumption tests, lists untestable assumptions, returns all results in r(),
*! preserves the estimation results (including e(sample)), and stops with an
*! error if a critical assumption is violated (unless -continue- is given).
*!
*! r() handling: the r() results present before the call (for example r(table)
*! from the estimation output) are kept. By default the assumptiontest results
*! are added to them (assumptiontest's value wins on a name clash); with
*! -noreturn- r() is left exactly as it was. r() is set before the error stop.
*!
*! Extending: for another estimation command  xyz  provide  _atest_xyz.ado
*! (see the help file, section "Extending assumptiontest").
program define assumptiontest
    version 17
    syntax [, NORETurn * ]

    * set the caller's r() aside; it is put back on every exit path below
    capture macro drop ATESTRC
    tempname rold rnew
    _return hold `rold'
    capture noisily _atest_main, `options'
    local rc = _rc
    local stoprc 0
    if "$ATESTRC" != "" local stoprc = $ATESTRC
    capture macro drop ATESTRC

    if `rc' {
        * error before or outside the tests (syntax, no estimation, ...)
        _return restore `rold'
        exit `rc'
    }
    if "`noreturn'" != "" {
        _return restore `rold'
    }
    else {
        _return hold `rnew'
        _return restore `rold'
        _atest_merge `rnew'
    }
    if `stoprc' exit `stoprc'
end

* ---- r() = caller's r() plus the assumptiontest results (held as `rnew').
*      Called with the caller's r() active; on a name clash the assumptiontest
*      result replaces the earlier one.
program define _atest_merge, rclass
    version 17
    args rnew
    return add
    _return restore `rnew'
    local sc : r(scalars)
    local mc : r(macros)
    local mt : r(matrices)
    foreach s of local sc {
        return scalar `s' = r(`s')
    }
    foreach m of local mc {
        return local `m' `"`r(`m')'"'
    }
    foreach m of local mt {
        tempname M
        matrix `M' = r(`m')
        return matrix `m' = `M'
    }
end

* ---- the tests themselves; returns its own results in r()
program define _atest_main, rclass
    version 17

    syntax [, CONTinue ALpha(real 0.05) GRAPHPath(string) NOOUTliers ///
        OUTcut(real 1) RESCut(real 3) INFLuence(string) LDcut(real -1) COOKcut(real 1) ///
        LINearity(string) CRITical(string) NONCRITical(string) * ]

    * ---- an estimation must be active
    if "`e(cmd)'" == "" {
        display as error "no estimation results found; fit a model first"
        exit 301
    }
    local cmd `e(cmd)'
    * stcox stores e(cmd) = cox (stcox_fr with shared frailty)
    if inlist("`cmd'", "cox", "stcox_fr") local cmd stcox
    * streg stores the distribution (ereg, weibull, ...) in e(cmd); e(cmd2) holds the command
    if inlist("`e(cmd2)'", "stcox", "streg") local cmd `e(cmd2)'
    if "`e(prefix)'" == "svy" {
        display as error "models fit with the svy: prefix are not covered by assumptiontest"
        exit 199
    }
    capture quietly findfile _atest_`cmd'.ado
    if _rc {
        display as error "assumptiontest does not cover estimation command {bf:`cmd'}"
        display as error "covered: stcox, streg, nbreg (see {help assumptiontest##extend:Extending assumptiontest})"
        exit 199
    }
    if "`influence'" == "" local influence ld
    if !inlist("`influence'", "ld", "dfbeta") {
        display as error "influence() must be ld or dfbeta"
        exit 198
    }
    if "`linearity'" == "" local linearity rcs
    if !inlist("`linearity'", "rcs", "quadratic") {
        display as error "linearity() must be rcs or quadratic"
        exit 198
    }
    if `alpha' <= 0 | `alpha' >= 1 {
        display as error "alpha() must be between 0 and 1"
        exit 198
    }

    * ---- spline() groups (repeatable): base variable first, then higher-order terms
    macro drop ATEST_*
    global ATEST_nsp 0
    local cn : colnames e(b)
    local rest = strtrim(`"`options'"')
    local allspl ""
    while `"`rest'"' != "" {
        local p = strpos(lower(`"`rest'"'), "spline(")
        if `p' != 1 {
            display as error `"option {bf:`rest'} not allowed"'
            exit 198
        }
        local after = substr(`"`rest'"', 8, .)
        local q = strpos(`"`after'"', ")")
        if `q' == 0 {
            display as error "spline(): missing closing parenthesis"
            exit 198
        }
        local grp = strtrim(substr(`"`after'"', 1, `q' - 1))
        local rest = strtrim(substr(`"`after'"', `q' + 1, .))
        if `: word count `grp'' < 2 {
            display as error "spline() needs the base variable followed by at least one higher-order term"
            exit 198
        }
        foreach v of local grp {
            capture confirm numeric variable `v'
            if _rc {
                display as error "spline(): {bf:`v'} is not a numeric variable in the dataset"
                exit 198
            }
            if `: list posof "`v'" in cn' == 0 {
                display as error "spline(): {bf:`v'} is not a term of the fitted model"
                exit 198
            }
            if `: list posof "`v'" in allspl' > 0 {
                display as error "spline(): {bf:`v'} appears in more than one group"
                exit 198
            }
        }
        local allspl `allspl' `grp'
        global ATEST_nsp = $ATEST_nsp + 1
        global ATEST_sp${ATEST_nsp} `grp'
    }

    * ---- graph folder
    if `"`graphpath'"' != "" {
        local last = substr(`"`graphpath'"', -1, 1)
        if inlist(`"`last'"', "/", "\") local graphpath = substr(`"`graphpath'"', 1, strlen(`"`graphpath'"') - 1)
        capture mkdir `"`graphpath'"'
        capture confirm new file `"`graphpath'/_atest_check.tmp"'
        if _rc {
            display as error `"graphpath(): cannot write to {bf:`graphpath'}"'
            exit 603
        }
    }

    * ---- global state used by the modules (removed again at the end)
    global ATEST_alpha     `alpha'
    global ATEST_outcut    `outcut'
    global ATEST_rescut    `rescut'
    global ATEST_infl      `influence'
    global ATEST_ldcut     `ldcut'
    global ATEST_cookcut   `cookcut'
    global ATEST_lin       `linearity'
    global ATEST_graphpath `"`graphpath'"'
    global ATEST_crit      `critical'
    global ATEST_noncrit   `noncritical'
    global ATEST_nooutliers `nooutliers'
    global ATEST_n 0
    global ATEST_names ""
    capture matrix drop _atest_res

    * ---- preserve the estimation results and e(sample)
    capture quietly estimates drop __atest_b
    quietly estimates store __atest_b
    tempvar esamp
    quietly generate byte `esamp' = e(sample)
    global ATEST_base  __atest_b
    global ATEST_esamp `esamp'
    quietly count if `esamp'
    local Nsample = r(N)

    * ---- classify the model terms
    _atest_terms
    global ATEST_cont `r(cont)'
    global ATEST_cat  `r(cat)'
    global ATEST_spl  `r(spl)'

    * ---- header
    local title `"`e(title)'"'
    if `"`title'"' == "" local title "`cmd'"
    local alphatxt = strtrim(string(`alpha', "%9.0g"))
    if substr("`alphatxt'", 1, 1) == "." local alphatxt "0`alphatxt'"
    display
    display as text "{hline 78}"
    display as text "assumptiontest: " as result `"`title'"' as text "  (`cmd')"
    display as text "  Estimation sample: " as result `Nsample' as text " observations;" ///
        as text "  alpha = " as result "`alphatxt'" as text ";  influence: " as result "`influence'" ///
        as text ";  linearity: " as result "`linearity'"
    local contdisp "$ATEST_cont"
    local catdisp  "$ATEST_cat"
    if "`contdisp'" == "" local contdisp "(none)"
    if "`catdisp'"  == "" local catdisp  "(none)"
    display as text "  Plain continuous: " as result "`contdisp'" ///
        as text "  Categorical: " as result "`catdisp'"
    forvalues g = 1/$ATEST_nsp {
        display as text "  Spline group `g': " as result "${ATEST_sp`g'}" as text " (base first)"
    }
    display as text "{hline 78}"
    display as text "  [VIOLATED] = critical assumption rejected (stops execution unless -continue-)"
    display as text "  [FLAG] = non-critical test rejected;  [NOT TESTABLE] = cannot be tested on the data"
    display as text "{hline 78}"

    * ---- run the module (errors are caught so that state can be cleaned up)
    capture noisily _atest_`cmd'
    local modrc = _rc

    * ---- restore estimation results and e(sample)
    capture _atest_restore
    quietly count if e(sample) & `esamp'
    local nback = r(N)
    quietly count if e(sample)
    local nback2 = r(N)
    local restored = (`nback' == `Nsample' & `nback2' == `Nsample')

    * ---- collect results
    local ntests = $ATEST_n
    local names  $ATEST_names
    tempname R
    local nviol 0
    local nflag 0
    local nuntest 0
    local nnotrun 0
    local violated ""
    local flagged  ""
    local untest   ""
    local notrun   ""
    if `ntests' > 0 {
        matrix `R' = _atest_res
        matrix rownames `R' = `names'
        matrix colnames `R' = critical code stat df p
        forvalues i = 1/`ntests' {
            local nm : word `i' of `names'
            local crit = `R'[`i', 1]
            local code = `R'[`i', 2]
            if `code' == 1 & `crit' == 1 {
                local ++nviol
                local violated `violated' `nm'
            }
            else if `code' == 1 {
                local ++nflag
                local flagged `flagged' `nm'
            }
            else if `code' == 2 {
                local ++nuntest
                local untest `untest' `nm'
            }
            else if `code' == 3 {
                local ++nnotrun
                local notrun `notrun' `nm'
            }
        }
        capture matrix drop atest_results
        matrix atest_results = `R'
    }

    * ---- summary
    display as text "{hline 78}"
    display as text "Summary:  " as result `nviol' as text " critical violation(s), " ///
        as result `nflag' as text " flag(s), " as result `nuntest' as text " not testable, " ///
        as result `nnotrun' as text " not run"
    if `nviol' > 0 display as error "  Critical violations: `violated'"
    if `nflag' > 0 display as result "  Flags: `flagged'"
    if `nnotrun' > 0 display as result "  Not run: `notrun'"
    if `restored' {
        display as text "  Estimation results and e(sample) restored; postestimation commands can be used."
    }
    else {
        display as error "  Warning: e(sample) could not be fully restored; re-run the estimation before postestimation."
    }
    display as text "{hline 78}"

    * ---- clean up (globals, stored estimates, helper variables)
    local alpha_ret = $ATEST_alpha
    capture estimates drop __atest_b
    capture drop __atest_dfb*
    capture matrix drop _atest_res
    macro drop ATEST_*

    * ---- return values
    return clear
    return local cmd        `"`cmd'"'
    return local model      `"`title'"'
    return local violated   `"`violated'"'
    return local flagged    `"`flagged'"'
    return local untestable `"`untest'"'
    return local notrun     `"`notrun'"'
    return local tests      `"`names'"'
    return scalar alpha        = `alpha_ret'
    return scalar N            = `Nsample'
    return scalar n_tests      = `ntests'
    return scalar n_violated   = `nviol'
    return scalar n_flagged    = `nflag'
    return scalar n_untestable = `nuntest'
    return scalar n_notrun     = `nnotrun'
    return scalar passed       = (`nviol' == 0)
    if `ntests' > 0 {
        forvalues i = 1/`ntests' {
            local nm : word `i' of `names'
            local pv = `R'[`i', 5]
            if `pv' < . return scalar p_`nm' = `pv'
        }
        return matrix results = `R'
    }

    * ---- stop rules: decided here, carried out by the wrapper after r() is set
    global ATESTRC 0
    if `modrc' {
        display as error "assumptiontest: a test module stopped with return code `modrc'"
        global ATESTRC `modrc'
    }
    else if `nviol' > 0 {
        if "`continue'" == "" {
            display as error "assumptiontest: `nviol' critical assumption(s) violated (`violated')"
            display as error "  results are kept in r() and matrix atest_results; use option {bf:continue} to proceed anyway"
            global ATESTRC 459
        }
        else {
            display as text "  (option continue specified: execution continues despite violated assumptions)"
        }
    }
end
