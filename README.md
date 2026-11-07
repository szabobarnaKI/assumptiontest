# assumptiontest
A Stata command that runs and evaluates model assumption tests after regress, logit, logistic, firthlogit (SSC), poisson, nbreg, stcox, streg (exponential and Weibull): proportional hazards, linearity, distributional form, equidispersion, constant error variance, separation, calibration, over-dispersion and influential observations, with automatic decisions, optional graphs saved to disk only, and an error stop if a critical assumption is violated. The estimation results and e(sample) are preserved for postestimation commands.

To install the command, run the following command in Stata:
##
    net from "https://szabobarnaki.github.io/assumptiontest/"
