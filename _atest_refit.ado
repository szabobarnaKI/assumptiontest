*! _atest_refit 0.3.0  re-run the last estimation command in a modified form
*! Internal helper of assumptiontest. Splits e(cmdline) into its parts
*! (varlist, weight, if/in, options), optionally removes options, adds
*! extra covariates / options, switches the command, and runs it quietly.
*! Options with arguments are removed with balanced parentheses, so nested
*! arguments such as texp(ln(_t)) or vce(cluster id) are removed completely.
*! Commands that do not store e(cmdline) (some community-contributed ones):
*! the command line is rebuilt from e(cmd), e(depvar), the terms of the first
*! equation of e(b), the weights, and  if <estimation sample>.
*! Returns r(rc) (return code of the refit) and r(cmdline) (command executed).
program define _atest_refit, rclass
    version 17
    syntax [, NEWcmd(string) EXTRA(string) REMove(string) ADDopts(string asis) ///
        NOWT NOIFIN ]

    local oldcmd `e(cmd)'
    if inlist("`oldcmd'", "cox", "stcox_fr") local oldcmd stcox
    if inlist("`e(cmd2)'", "stcox", "streg") local oldcmd `e(cmd2)'
    local newcmd = cond("`newcmd'" == "", "`oldcmd'", "`newcmd'")

    local s = strtrim(`"`e(cmdline)'"')
    if `"`s'"' == "" {
        * ---- rebuild the command line from the stored results
        tempname bb
        matrix `bb' = e(b)
        local cn : colnames `bb'
        local ce : coleq `bb'
        local eq1 : word 1 of `ce'
        local xs ""
        local hascons 0
        forvalues j = 1/`= colsof(`bb')' {
            local t : word `j' of `cn'
            local q : word `j' of `ce'
            if "`q'" != "`eq1'" continue
            if "`t'" == "_cons" {
                local hascons 1
                continue
            }
            local xs `xs' `t'
        }
        local s `"`oldcmd' `e(depvar)' `xs'"'
        if "`e(wtype)'" != "" local s `"`s' [`e(wtype)' `e(wexp)']"'
        local s `"`s' if $ATEST_esamp"'
        if !`hascons' local s `"`s', noconstant"'
    }

    * ---- split the command line at top level (parentheses depth 0)
    local sp = strpos(`"`s'"', " ")
    local rest = strtrim(substr(`"`s'"', `sp' + 1, .))
    local n = strlen(`"`rest'"')
    local depth 0
    local p_ifin 0
    local p_wt   0
    local p_op   0
    forvalues i = 1/`n' {
        local c = substr(`"`rest'"', `i', 1)
        if "`c'" == "(" {
            local ++depth
        }
        else if "`c'" == ")" {
            local --depth
        }
        else if `depth' == 0 {
            if "`c'" == "," {
                local p_op `i'
                continue, break
            }
            if "`c'" == "[" & `p_wt' == 0 local p_wt `i'
            if `p_ifin' == 0 {
                local w4 = substr(`"`rest'"', `i', 4)
                if inlist("`w4'", " if ", " if(", " in ") local p_ifin `i'
            }
        }
    }
    local firststart = `n' + 1
    foreach p in p_ifin p_wt p_op {
        if ``p'' > 0 & ``p'' < `firststart' local firststart ``p''
    }
    local vars = strtrim(substr(`"`rest'"', 1, `firststart' - 1))
    local ifin ""
    local wt   ""
    local opts ""
    if `p_ifin' > 0 {
        local e = `n' + 1
        foreach p in p_wt p_op {
            if ``p'' > `p_ifin' & ``p'' < `e' local e ``p''
        }
        local ifin = strtrim(substr(`"`rest'"', `p_ifin', `e' - `p_ifin'))
    }
    if `p_wt' > 0 {
        local e = `n' + 1
        foreach p in p_ifin p_op {
            if ``p'' > `p_wt' & ``p'' < `e' local e ``p''
        }
        local wt = strtrim(substr(`"`rest'"', `p_wt', `e' - `p_wt'))
    }
    if `p_op' > 0 {
        local opts = strtrim(substr(`"`rest'"', `p_op' + 1, .))
    }

    * ---- remove options
    if "`remove'" == "all" {
        local opts ""
    }
    else {
        foreach o of local remove {
            if inlist("`o'", "robust", "time", "nohr", "noshow", "irr", "nolrtest", "hc2", "hc3", "or") {
                local opts = ustrregexra(`"`opts'"', "(^|\s)`o'(\s|$)", " ")
            }
            else {
                * option with an argument: delete from its name to the matching ")"
                local guard 0
                while ustrregexm(`"`opts'"', "(^|\s)`o'\s*\(") & `guard' < 20 {
                    local ++guard
                    local m = ustrregexs(0)
                    local st = strpos(`"`opts'"', `"`m'"')
                    local op = `st' + strlen(`"`m'"') - 1
                    local L = strlen(`"`opts'"')
                    local d 0
                    local en = `L'
                    forvalues j = `op'/`L' {
                        local c = substr(`"`opts'"', `j', 1)
                        if "`c'" == "(" {
                            local ++d
                        }
                        else if "`c'" == ")" {
                            local --d
                            if `d' == 0 {
                                local en `j'
                                continue, break
                            }
                        }
                    }
                    local opts = substr(`"`opts'"', 1, `st' - 1) + " " + substr(`"`opts'"', `en' + 1, .)
                }
            }
        }
    }
    local opts = strtrim(itrim(`"`opts'"'))
    if `"`addopts'"' != "" local opts `"`opts' `addopts'"'

    if "`nowt'"   != "" local wt ""
    if "`noifin'" != "" local ifin ""

    local full `"`newcmd' `vars' `extra' `wt' `ifin'"'
    if `"`opts'"' != "" local full `"`full', `opts'"'

    capture quietly `full'
    return scalar rc = _rc
    return local cmdline `"`full'"'
end
