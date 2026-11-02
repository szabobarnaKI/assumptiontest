# assumptiontest
A Stata command that runs and evaluates model assumption tests after stcox, streg (exponential and Weibull) and nbreg: proportional hazards, linearity, distributional form, calibration, over-dispersion and influential observations, with automatic decisions, optional graphs saved to disk only, and an error stop if a critical assumption is violated. The estimation results and e(sample) are preserved for postestimation commands.

To install the command, run the following command in Stata:
##
    net from "https://szabobarnaki.github.io/assumptiontest/"
