*! _atest_post 0.1.0  records and displays one assumption-test result
*! Internal helper of assumptiontest. State is kept in globals ATEST_* and
*! the global matrix _atest_res (columns: critical, code, stat, df, p).
*! code: 0 = passed, 1 = violated/flagged, 2 = not testable, 3 = not run, 4 = info
program define _atest_post
    version 17
    syntax , NAME(string) DESC(string) [ P(string) STAT(string) DF(string) ///
        STATName(string) CRITical FLAG TYPE(string) NOte(string) ]

    * ---- clean row name (valid matrix row name, unique, <=32 chars)
    local nm = subinstr("`name'", ".", "_", .)
    local nm = subinstr("`nm'", "#", "X", .)
    local nm = subinstr("`nm'", " ", "_", .)
    local nm = substr("`nm'", 1, 30)
    local base `nm'
    local k 1
    local usednames $ATEST_names
    while `: list posof "`nm'" in usednames' > 0 {
        local ++k
        local nm `base'`k'
    }

    * ---- critical / non-critical, with user overrides
    local crit = ("`critical'" != "")
    local critlist $ATEST_crit
    local noncritlist $ATEST_noncrit
    if `: list posof "`name'" in critlist' > 0    local crit 1
    if `: list posof "`name'" in noncritlist' > 0 local crit 0

    * ---- numeric parts
    foreach v in p stat df {
        if "``v''" == "" {
            local `v' .
        }
        else {
            local `v' = real("``v''")
        }
    }

    * ---- code
    local code 0
    if "`type'" == "untestable" {
        local code 2
    }
    else if "`type'" == "notrun" {
        local code 3
    }
    else if "`type'" == "info" {
        local code 4
    }
    else if "`flag'" != "" {
        local code 1
    }
    else if `p' < . {
        if `p' < $ATEST_alpha {
            local code 1
        }
    }

    * ---- store
    tempname row
    matrix `row' = (`crit', `code', `stat', `df', `p')
    if $ATEST_n == 0 {
        matrix _atest_res = `row'
    }
    else {
        matrix _atest_res = _atest_res \ `row'
    }
    global ATEST_n = $ATEST_n + 1
    global ATEST_names `"$ATEST_names `nm'"'
    global ATEST_d$ATEST_n `"`desc'"'

    * ---- display
    local tag  "OK"
    local col  "text"
    if `code' == 1 & `crit' {
        local tag "VIOLATED"
        local col "error"
    }
    else if `code' == 1 {
        local tag "FLAG"
        local col "result"
    }
    else if `code' == 2 {
        local tag "NOT TESTABLE"
    }
    else if `code' == 3 {
        local tag "NOT RUN"
        local col "result"
    }
    else if `code' == 4 {
        local tag "INFO"
    }
    local res ""
    if `stat' < . {
        local sn = cond("`statname'" == "", "stat", "`statname'")
        local res `"`sn' = `=string(`stat', "%9.3f")'"'
    }
    if `p' < . {
        local ptxt = cond(`p' < 0.001, "<0.001", string(`p', "%5.3f"))
        local res `"`res'  p = `ptxt'"'
    }
    local ctag = cond(`crit', "critical", "non-critical")
    display as `col' "  [`tag']" as text " `desc'"
    if `"`res'"' != "" | `"`note'"' != "" {
        display as text "      `res'" _c
        if `code' < 2 display as text "  (`ctag')" _c
        if `"`note'"' != "" display as text "  `note'" _c
        display
    }
end
