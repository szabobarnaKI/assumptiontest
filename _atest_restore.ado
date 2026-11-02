*! _atest_restore 0.1.0  put the original estimation results (and e(sample)) back
*! Internal helper of assumptiontest. Uses globals ATEST_base (name of the
*! stored estimates) and ATEST_esamp (0/1 variable holding the estimation sample).
program define _atest_restore
    version 17
    quietly estimates restore $ATEST_base
    capture quietly ereturn repost, esample($ATEST_esamp)
end
