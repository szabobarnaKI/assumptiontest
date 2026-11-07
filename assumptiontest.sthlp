{smcl}
{* *! assumptiontest 0.4.1}{...}
{title:Title}

{phang}
{bf:assumptiontest} {hline 2} Test the assumptions of the most recently fitted stcox, streg, nbreg, poisson,
regress, logit, logistic or firthlogit model


{title:Syntax}

{p 8 17 2}
{cmd:assumptiontest}
[{cmd:,} {it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt cont:inue}}do not stop with an error when a critical assumption is violated{p_end}
{synopt:{opt spline(varlist)}}declare one group of spline variables; repeatable{p_end}
{synopt:{opt al:pha(#)}}significance level for all tests; default {cmd:alpha(0.05)}{p_end}
{synopt:{opt graphp:ath(path)}}write diagnostic graphs to this folder only; none are kept in memory{p_end}
{synopt:{opt lin:earity(method)}}linearity test: {cmd:rcs} (default) or {cmd:quadratic}{p_end}
{synopt:{opt infl:uence(measure)}}influence measure (stcox, streg): {cmd:ld} (default) or {cmd:dfbeta}{p_end}
{synopt:{opt ldcut(#)}}likelihood-displacement cut-off; default is the chi-squared median (50% confidence region){p_end}
{synopt:{opt outcut(#)}}dfbeta threshold in standard errors ({cmd:influence(dfbeta)}); default {cmd:outcut(1)}{p_end}
{synopt:{opt cook:cut(#)}}Cook's distance / delta-beta cut-off (nbreg, poisson, regress, logit models); default {cmd:cookcut(1)}{p_end}
{synopt:{opt nooutliers}}skip the influence/outlier screening{p_end}
{synopt:{opt rescut(#)}}residual threshold for the informational outlier count; default {cmd:rescut(3)}{p_end}
{synopt:{opt vif:cut(#)}}largest acceptable variance inflation factor (regress); default {cmd:vifcut(10)}{p_end}
{synopt:{opt epv:cut(#)}}smallest acceptable number of events per coefficient (logit models); default {cmd:epvcut(10)}{p_end}
{synopt:{opt critical(namelist)}}treat these tests as critical{p_end}
{synopt:{opt noncritical(namelist)}}treat these tests as non-critical{p_end}
{synopt:{opt noret:urn}}leave {cmd:r()} exactly as it was before the call; do not add the results{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:assumptiontest} identifies the last estimation command from {cmd:e(cmd)} (and {cmd:e(cmd2)}
for {cmd:streg}), runs the automatic assumption tests relevant for it, and prints every assumption
of the model, including those that cannot be tested on the data, with the result of
each test. The {cmd:r()} results that existed before the call (for example {cmd:r(table)}
from the estimation output) are kept, and the test results are added to them; with
{cmd:noreturn}, {cmd:r()} is left unchanged. See {help assumptiontest##results:Stored results}.

{pstd}
Covered: {helpb regress}, {helpb logit}, {helpb logistic}, the
community-contributed {cmd:firthlogit} (Coveney; {cmd:ssc install firthlogit}), {helpb poisson}, {helpb nbreg} (default form; {cmd:dispersion(constant)}
is flagged and not tested), {helpb stcox}, {helpb streg} (exponential and Weibull; any other
distribution is flagged and not tested). Special structures handled: {cmd:strata()}, robust/cluster
variance, shared frailty, {cmd:tvc()}/{cmd:texp()}, {cmd:ancillary()}, {cmd:offset()}/{cmd:exposure()},
and multiple-record data from {helpb stsplit}.

{pstd}
Some diagnostics are not available for the fitted model itself (for example {cmd:predict, ldisplace}
after robust/cluster variance or {cmd:tvc()}, or {cmd:estat hettest} after {cmd:regress, vce(robust)}). They are then computed on a {it:diagnostic refit} of the
same model without the feature concerned, and the output says so in the note of that row.

{pstd}
The estimation results {bf:and} {cmd:e(sample)} are restored before the command
ends, so postestimation commands can be used afterwards. The dataset is never
modified (only temporary variables are used).

{pstd}
By default the command stops with error 459 after printing all results if a
{it:critical} test is rejected. Non-critical rejections are printed as
{bf:[FLAG]} and never stop execution. Option {cmd:continue} suppresses the stop.
No multiple-testing correction is applied.


{title:Options}

{phang}
{opt spline(varlist)} declares variables that belong to one spline (for example
the output of {helpb mkspline}). The {bf:base variable comes first}, higher-order terms
follow. The option may be repeated, once per spline. Variables in a spline group are
not tested for linearity; a joint Wald test of the higher-order terms is shown for
information. Whether the knots are adequate cannot be tested and is listed as such.

{phang}
{opt graphpath(path)} writes diagnostic graphs (currently the Schoenfeld-residual
plots for covariates that fail the proportional-hazards test after {cmd:stcox} or {cmd:streg};
no graphs are made after the other commands) to {it:path} as PNG
files. Each graph is dropped from memory immediately after export, so repeated calls
never reach Stata's limit on graphs in memory. Each graph must be drawn before it can be exported,
so a Graph window may flash briefly while the command runs. Without this option, no graphs are
drawn. {it:path} is created if it does not exist.

{phang}
{opt critical(namelist)} and {opt noncritical(namelist)} override the default
classification below; use the test names shown in the first column of {cmd:r(results)}
(for example {cmd:noncritical(ph_age lin_bmi)} or {cmd:critical(normal)}).

{phang}
{opt noreturn} leaves {cmd:r()} exactly as it was before {cmd:assumptiontest} was
called: nothing is added, removed or changed, so {cmd:r(table)} and any other results
of the preceding estimation remain available. The command then only prints its output.
The results table is still kept in matrix {cmd:atest_results}. Without {cmd:noreturn}
(the default), the earlier {cmd:r()} results are kept as well and the results listed
under {help assumptiontest##results:Stored results} are added.

{phang}
{opt linearity(method)} selects how linearity of plain continuous covariates is tested.
{cmd:rcs} (default) adds a restricted cubic spline with 3 knots (Stata {helpb mkspline},
{cmd:cubic nknots(3)}, knots at the 10th, 50th and 90th percentiles). The model already contains
the linear term, so only the single non-linear spline term is added and its test (1 df)
is the linearity test (a Wald chi-squared test; an F test after {cmd:regress}; a profile penalized likelihood-ratio
test after {cmd:firthlogit}). {cmd:quadratic} adds a centred quadratic term instead. If a 3-knot spline
cannot be built for a variable (too few distinct values), the quadratic term is used for that variable
and the output says so.

{phang}
{opt influence(measure)} selects the influence measure for outlier screening. {cmd:ld} (default) uses
{it:likelihood displacement}, LD_i = 2[l(b) - l(b_(i))], approximated for Cox models by
{cmd:predict, ldisplace}. An observation is flagged
when LD_i exceeds the median of the chi-squared distribution with as many degrees of freedom as there are
estimated coefficients (base levels of factor variables are not counted), that is, when deleting it moves the estimate outside the 50% confidence region
of the full-data estimate (Cook and Weisberg 1982; the same calibration underlies the Cook's distance > 1 rule used for
the other commands). The cut-off does not depend on {opt alpha()}. In large samples single observations rarely move the estimate this far, so the
screening flags only observations that on their own change the estimates substantially. {opt ldcut(#)} replaces this cut-off. These
options apply to {cmd:stcox} and {cmd:streg}; for the other commands see {opt cookcut(#)} below. {cmd:dfbeta} flags observations whose
deletion would shift any single coefficient by more than {opt outcut(#)} standard errors (default 1).

{phang}
{opt cookcut(#)} sets the cut-off for Cook's distance (default 1, clearly influential;
4/N is a common stricter screening value). After {cmd:nbreg} and {cmd:poisson}, Cook's distance comes from the equivalent
generalized linear model, {cmd:glm, family(nbinomial ml)} or {cmd:family(poisson)}, {cmd:link(log)}, refitted without robust/cluster
variance; after {cmd:regress} from {cmd:predict, cooksd}; after {cmd:logit} and {cmd:logistic} the same cut-off is applied to
Pregibon's delta-beta ({cmd:predict, dbeta}, by covariate pattern), the logistic analogue of Cook's distance; after {cmd:firthlogit}
to an observation-level delta-beta computed from the penalized fit.

{phang}
{opt rescut(#)} sets the threshold for the informational count of large residuals: |standardized deviance residual|
(nbreg, poisson), |studentized residual| (regress), |standardized Pearson residual| (logit models).

{phang}
{opt vifcut(#)} sets the variance-inflation-factor threshold for the multicollinearity flag after {cmd:regress} (default 10).
Interaction terms and variables in {cmd:spline()} groups are skipped because they are collinear by construction.

{phang}
{opt epvcut(#)} sets the events-per-coefficient threshold for the sparse-data flag after the logit models
(default 10, Peduzzi et al. 1996). Events per coefficient = min(number of events, number of non-events) / number of
estimated coefficients (constant and base levels not counted).


{title:Tests and default classification}

{pstd}
{bf:Critical} tests address the defining assumptions of the model. Everything else is
a non-critical flag.

{p2colset 5 30 32 2}{...}
{syntab:{it:regress}}
{p2col:{cmd:lin_}{it:var}}critical: linearity of the mean (F test of the non-linear spline or quadratic term){p_end}
{p2col:{cmd:homosked}}critical: constant error variance, Breusch-Pagan/Cook-Weisberg test against the fitted values ({cmd:estat hettest, iid}, Koenker's form, which does not assume normal errors). With robust/cluster/hc2/hc3 variance: information, computed on a refit with the default variance{p_end}
{p2col:{cmd:reset}}flag: Ramsey RESET ({cmd:estat ovtest}){p_end}
{p2col:{cmd:normal}}flag: normality of the residuals (skewness/kurtosis test, {cmd:sktest}); in large samples inference does not depend on it{p_end}
{p2col:{cmd:vif}}flag: largest variance inflation factor > {opt vifcut()}{p_end}
{p2col:{cmd:outliers}}flag: Cook's distance > {opt cookcut()}; {cmd:res_out} counts |studentized residual| > {opt rescut()} (information){p_end}

{syntab:{it:logit, logistic}}
{p2col:{cmd:separation}}critical: no complete or quasi-complete separation (no completely determined outcomes, {cmd:e(N_cds)}, {cmd:e(N_cdf)}, and no perfect predictors dropped, {cmd:e(rules)}); if violated, consider {cmd:firthlogit}{p_end}
{p2col:{cmd:epv}}flag: events per estimated coefficient < {opt epvcut()}{p_end}
{p2col:{cmd:lin_}{it:var}}critical: linearity of the log-odds (Wald test){p_end}
{p2col:{cmd:hl_gof}}flag: Hosmer-Lemeshow goodness of fit, 10 groups of predicted risk, chi2(G-2) (weights not used){p_end}
{p2col:{cmd:linktest}}flag: link test as in {helpb linktest} (refit on the linear predictor and its square; Wald test of the square), computed directly{p_end}
{p2col:{cmd:outliers}}flag: Pregibon's delta-beta > {opt cookcut()} (on a refit without robust/cluster variance if needed); {cmd:res_out} counts |standardized Pearson residual| > {opt rescut()} (information){p_end}

{syntab:{it:firthlogit}}
{p2col:{cmd:separation}}information: the Firth penalty gives finite estimates under separation{p_end}
{p2col:{cmd:lin_}{it:var}}critical: linearity of the log-odds; profile penalized likelihood-ratio test: the model with the added term is compared with the same model with that coefficient constrained to 0, so both share the same Firth penalty (Heinze and Schemper 2002); two refits per variable. Wald tests are unreliable under separation. If {cmd:firthlogit} does not accept {cmd:constraints()}, the Wald test is used and the note says so{p_end}
{p2col:{cmd:epv}, {cmd:hl_gof}, {cmd:linktest}}flags: as for logit{p_end}
{p2col:{cmd:outliers}}flag: observation-level delta-beta from the penalized fit, r^2 h/(1-h)^2 with Pearson residual r and leverage h = p(1-p) x'Vx{p_end}

{syntab:{it:poisson}}
{p2col:{cmd:overdisp}}critical: equidispersion, Var(y) = mu; Cameron-Trivedi auxiliary regression of ((y-mu)^2-y)/mu on mu without constant, robust SE, two-sided F test (detects over- and under-dispersion; the direction is shown). With robust/cluster variance: information (Poisson quasi-likelihood inference does not need it){p_end}
{p2col:{cmd:lin_}{it:var}}critical: linearity of the log mean{p_end}
{p2col:{cmd:pearson}, {cmd:zeros}}flags: as for nbreg (expected zeros exp(-mu)){p_end}
{p2col:{cmd:outliers}}flag: Cook's distance from the equivalent {cmd:glm} (Poisson, log link); {cmd:dev_out} as for nbreg{p_end}

{syntab:{it:nbreg}}
{p2col:{cmd:lin_}{it:var}}critical: linearity of the log mean (Wald test of the non-linear spline or quadratic term){p_end}
{p2col:{cmd:overdisp}}LR test of alpha=0 (one-sided, mixture chi-squared). Flag if Poisson cannot be rejected (a Poisson model might suffice); information if over-dispersion is present. With robust/cluster variance the test is taken from a refit with the default variance (same estimates){p_end}
{p2col:{cmd:varfun}}flag: approximate score-type test of the NB2 variance function, Var = mu + alpha*mu^2 (auxiliary regression with robust SE; F test){p_end}
{p2col:{cmd:pearson}}flag: Pearson chi-squared goodness of fit (approximate; poor for very small means){p_end}
{p2col:{cmd:zeros}}flag: one-sided z test for excess zeros{p_end}
{p2col:{cmd:outliers}}flag: Cook's distance from the equivalent {cmd:glm} (negative binomial, log link); {cmd:dev_out} counts large standardized deviance residuals (information){p_end}

{syntab:{it:stcox}}
{p2col:{cmd:ph_}{it:var}, {cmd:ph_global}}critical: proportional hazards, Schoenfeld residuals ({cmd:estat phtest}). With robust/cluster/shared(), a diagnostic refit without them is used. With {cmd:tvc()} the Wald test of the time-varying terms
({cmd:ph_tvc}) is used and other covariates are not testable.{p_end}
{p2col:{cmd:lin_}{it:var}}critical: linearity of the log hazard; Wald test of the non-linear term of a 3-knot restricted cubic spline (or of a quadratic term, {cmd:linearity(quadratic)}), all added in one refit (the fitted
robust/strata/frailty structure is kept){p_end}
{p2col:{cmd:gb_fit}}flag: observed vs expected events across risk groups (Gronnesby-Borgan type); not run after {cmd:tvc()}, because martingale residuals are not available and a model without the time-varying terms would be a different
model{p_end}
{p2col:{cmd:outliers}}flag: influence screening by likelihood displacement (default) or dfbeta; after {cmd:tvc()}, robust/cluster variance or shared(), computed on a diagnostic refit without them {p_end}
{p2col:{cmd:spl_}{it:var}}information: nonlinearity of a spline{p_end}

{syntab:{it:streg}}
{p2col:{cmd:ph_}{it:var}, {cmd:ph_global}}critical: Schoenfeld test on the equivalent {cmd:stcox} model (both distributions are PH models; options and weights are dropped in this refit, except {cmd:strata()}, which becomes
{cmd:stcox, strata()}){p_end}
{p2col:{cmd:const_haz}}critical, exponential: Wald test of Weibull shape p=1 (constant hazard); with {cmd:strata()}, joint test in all strata{p_end}
{p2col:{cmd:weib_fit}}critical, Weibull: Wald test of kappa=1 in the generalized gamma (AFT metric); not run if ancillary()/strata() or shared()/frailty() are used, or if the refit does not converge{p_end}
{p2col:{cmd:shape1}}information, Weibull: p=1 (exponential would suffice); with ancillary()/strata(), joint Wald test of all ln_p coefficients including the constant{p_end}
{p2col:{cmd:lin_}{it:var}}critical: linearity (as for stcox; log-hazard, or log-time with the {cmd:time} option){p_end}
{p2col:{cmd:gb_fit}, {cmd:outliers}}flags: as for stcox (influence is computed on the equivalent Cox model){p_end}

{pstd}
Linearity is tested only for plain continuous covariates (one refit per call; two refits per variable after {cmd:firthlogit}): variables entered as
factors, variables in a {cmd:spline()} group, 0/1 variables and interaction terms are
not tested. Every other variable in the model is treated as plain continuous, so a categorical
variable entered without the {cmd:i.} prefix is tested as if it were continuous.

{pstd}
{bf:Not testable} assumptions (always listed where they apply): non-informative censoring (survival models),
independence of observations/clustering, omitted covariates and correct specification of the
conditional mean (for regress also exogeneity and covariates measured without error), spline knot adequacy, frailty distribution (shared frailty), strata capturing all
baseline-hazard differences (stratified Cox models), multiple records per subject carrying no extra
information (multiple-record data), gamma heterogeneity (NB2), exposure/offset coefficient fixed at 1.


{title:Visual evaluation}

{pstd}
All tests are automatic. The only graphs are optional Schoenfeld-residual plots
({cmd:graphpath()}) for covariates whose PH test failed. Interpretation: the smooth in
the plot should be roughly flat around the line of the estimated coefficient;
a clear upward or downward trend over time means the hazard ratio changes with time.
The test decision does not depend on the graph.


{marker methods}{...}
{title:Suggested text for a methods section}

{pstd}
The text below describes what {cmd:assumptiontest} does. Delete the parts that do not apply to
the models you report, and fill in the bracketed items.

{pstd}
{it:Linear regression.} "Linearity of the mean in continuous covariates was assessed with a three-knot restricted cubic spline as above
(F test). Constant error variance was tested with the Breusch-Pagan/Cook-Weisberg test against the fitted values in Koenker's
form (Breusch and Pagan 1979; Cook and Weisberg 1983; Koenker 1981) [; with robust variance estimates this was reported for information only],
the functional form with Ramsey's RESET test (Ramsey 1969), normality of the residuals with the skewness-kurtosis test
(D'Agostino, Belanger and D'Agostino 1990), multicollinearity with variance inflation factors (threshold 10), and influential
observations with Cook's distance (Cook 1977)."

{pstd}
{it:Logistic regression.} "Complete and quasi-complete separation was checked from the observations whose outcome was completely
determined and from perfect predictors (Albert and Anderson 1984), and the number of events per estimated coefficient was compared with 10
(Peduzzi et al. 1996). Linearity of the log-odds was assessed with a three-knot restricted cubic spline as above. Calibration was assessed with the
Hosmer-Lemeshow test with ten groups of predicted risk (Hosmer, Lemeshow and Sturdivant 2013), the link function with the link test
(Pregibon 1980), and influential observations with Pregibon's delta-beta (Pregibon 1981). [For Firth's penalized-likelihood logistic
regression (Firth 1993; Heinze and Schemper 2002), linearity was tested with profile penalized likelihood-ratio tests and delta-beta was computed
from the penalized fit.]"

{pstd}
{it:Poisson regression.} "Equidispersion was tested with the auxiliary-regression test of Cameron and Trivedi (1990)
[; with robust variance estimates this assumption is not required and the test is reported for information].
Linearity of the log mean, overall fit, excess zeros and influential observations were assessed as for negative binomial regression,
with the Poisson generalized linear model."

{pstd}
{it:Negative binomial regression.} "Over-dispersion relative to the Poisson model was assessed with
the likelihood-ratio test of alpha = 0, using the mixture chi-squared reference distribution (Self and Liang 1987). Linearity
of the log mean was assessed with a three-knot restricted cubic spline as above. The negative binomial
variance function (variance = mean + alpha x mean squared) was checked with a score-type test, overall fit with
the Pearson chi-squared statistic, excess zeros by comparing observed and model-predicted numbers of zero counts
(Cameron and Trivedi 2013), and influential observations by Cook's distance and standardized deviance residuals from the equivalent
generalized linear model (Pregibon 1981; McCullagh and Nelder 1989). The independence of observations and the correct specification of the conditional mean
cannot be fully tested with the data."

{pstd}
{it:Cox regression.} "The proportional-hazards assumption was assessed for each covariate and
globally with tests based on scaled Schoenfeld residuals (Grambsch and Therneau 1994). For models with
robust or clustered variance estimates or shared frailty, the test was carried out on the corresponding
model without these features. [For models with time-varying coefficients, the Wald test of the
time-by-covariate interaction terms was used.] The linearity of the log hazard in continuous covariates
was assessed by adding a restricted cubic spline with three knots (at the 10th, 50th and 90th percentiles)
and testing its non-linear term with a Wald test (Harrell 2015); [alternatively: by adding a quadratic term]. Overall
calibration was assessed by comparing observed and expected numbers of events across ten risk-score
groups (Gronnesby and Borgan 1996; May and Hosmer 1998). Influential observations were identified by the
likelihood displacement (Cook 1986; Collett 2015), flagging observations whose deletion would move the
coefficient vector outside the 50% confidence region of the full-data estimate (Cook and Weisberg 1982) [; for models with robust variance,
shared frailty or time-varying coefficients, on the corresponding model without these features]. Non-informative censoring, independence
of observations and the absence of important omitted covariates cannot be tested with the data and were
judged on subject-matter grounds."

{pstd}
{it:Parametric survival models (exponential, Weibull).} "Proportional hazards was assessed on the
equivalent Cox model as above. The constant-hazard assumption of the exponential model was tested with a Wald test of
the Weibull shape parameter (p = 1); the adequacy of the Weibull model was tested with a Wald test of
kappa = 1 in the generalized gamma model, which contains the Weibull as a special case. Linearity, calibration
and influential observations were assessed as for Cox regression (influence on the equivalent Cox model)."

{pstd}
{it:All models.} "All tests were two-sided [except where noted] at the 5% level without adjustment for multiple
testing. Tests of the defining model assumptions (proportional hazards, linearity, distributional form, equidispersion,
constant error variance, absence of separation) were treated as critical;
other tests were treated as diagnostic flags. Analyses were performed in Stata [version] with the user-written
command assumptiontest [version, source]."

{pstd}
{bf:Limitations to mention if relevant.} The calibration tests, the variance-function test and the Pearson test
are approximate; Cook's distance, delta-beta, VIF and events-per-variable cut-offs are conventions. Because no multiple-testing correction is applied, occasional
rejections are expected by chance when many covariates are tested.

{pstd}
{bf:References for the methods text}

{phang}Albert A, Anderson JA. On the existence of maximum likelihood estimates in logistic regression models. {it:Biometrika} 1984; 71: 1-10.{p_end}
{phang}Breusch TS, Pagan AR. A simple test for heteroscedasticity and random coefficient variation. {it:Econometrica} 1979; 47: 1287-1294.{p_end}
{phang}Cameron AC, Trivedi PK. Regression-based tests for overdispersion in the Poisson model. {it:Journal of Econometrics} 1990; 46: 347-364.{p_end}
{phang}Cameron AC, Trivedi PK. {it:Regression Analysis of Count Data}, 2nd ed. Cambridge University Press, 2013.{p_end}
{phang}Collett D. {it:Modelling Survival Data in Medical Research}, 3rd ed. CRC Press, 2015.{p_end}
{phang}Cook RD. Detection of influential observation in linear regression. {it:Technometrics} 1977; 19: 15-18.{p_end}
{phang}Cook RD, Weisberg S. {it:Residuals and Influence in Regression}. Chapman and Hall, 1982.{p_end}
{phang}Cook RD, Weisberg S. Diagnostics for heteroscedasticity in regression. {it:Biometrika} 1983; 70: 1-10.{p_end}
{phang}Cook RD. Assessment of local influence. {it:Journal of the Royal Statistical Society B} 1986; 48: 133-169.{p_end}
{phang}D'Agostino RB, Belanger A, D'Agostino RB Jr. A suggestion for using powerful and informative tests of normality. {it:The American Statistician} 1990; 44: 316-321.{p_end}
{phang}Firth D. Bias reduction of maximum likelihood estimates. {it:Biometrika} 1993; 80: 27-38.{p_end}
{phang}Grambsch PM, Therneau TM. Proportional hazards tests and diagnostics based on weighted residuals. {it:Biometrika} 1994; 81: 515-526.{p_end}
{phang}Gronnesby JK, Borgan O. A method for checking regression models in survival analysis based on the risk score. {it:Lifetime Data Analysis} 1996; 2: 315-328.{p_end}
{phang}Harrell FE. {it:Regression Modeling Strategies}, 2nd ed. Springer, 2015.{p_end}
{phang}Heinze G, Schemper M. A solution to the problem of separation in logistic regression. {it:Statistics in Medicine} 2002; 21: 2409-2419.{p_end}
{phang}Hosmer DW, Lemeshow S, Sturdivant RX. {it:Applied Logistic Regression}, 3rd ed. Wiley, 2013.{p_end}
{phang}Koenker R. A note on studentizing a test for heteroscedasticity. {it:Journal of Econometrics} 1981; 17: 107-112.{p_end}
{phang}McCullagh P, Nelder JA. {it:Generalized Linear Models}, 2nd ed. Chapman and Hall, 1989.{p_end}
{phang}May S, Hosmer DW. A simplified method of calculating an overall goodness-of-fit test for the Cox proportional hazards
model. {it:Lifetime Data Analysis} 1998; 4: 109-120.{p_end}
{phang}Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR. A simulation study of the number of events per variable in logistic regression analysis.
{it:Journal of Clinical Epidemiology} 1996; 49: 1373-1379.{p_end}
{phang}Pregibon D. Goodness of link tests for generalized linear models. {it:Applied Statistics} 1980; 29: 15-24.{p_end}
{phang}Pregibon D. Logistic regression diagnostics. {it:Annals of Statistics} 1981; 9: 705-724.{p_end}
{phang}Ramsey JB. Tests for specification errors in classical linear least-squares regression analysis. {it:Journal of the Royal Statistical Society B} 1969; 31: 350-371.{p_end}
{phang}Self SG, Liang KY. Asymptotic properties of maximum likelihood estimators and likelihood ratio tests under nonstandard
conditions. {it:Journal of the American Statistical Association} 1987; 82: 605-610.{p_end}


{marker extend}{...}
{title:Extending assumptiontest}

{pstd}
For an estimation command {it:xyz}, create {cmd:_atest_xyz.ado} on the adopath;
{cmd:assumptiontest} calls it automatically when {cmd:e(cmd)} is {it:xyz}. A command that fits the same model as an
existing module can use a one-line module (as {cmd:_atest_logistic.ado} and {cmd:_atest_firthlogit.ado} call {cmd:_atest_logit}). If the command
stores a different name in {cmd:e(cmd)} (as {cmd:stcox} stores {cmd:cox} and {cmd:streg} stores
the distribution), add a mapping line near the top of {cmd:assumptiontest.ado} and the same line in
{cmd:_atest_refit.ado}. The module can use these helpers:

{phang2}{cmd:_atest_post, name(}{it:id}{cmd:) desc(}{it:text}{cmd:)} [{cmd:p(}{it:#}{cmd:) stat(}{it:#}{cmd:) df(}{it:#}{cmd:) statname(}{it:text}{cmd:) critical flag type(untestable|notrun|info) note(}{it:text}{cmd:)}]{break}
records and prints one result; {cmd:critical} marks it as critical, {cmd:type()} records assumptions that are not
testable, not run, or informational.{p_end}
{phang2}{cmd:_atest_refit, [newcmd() extra() remove() addopts() nowt noifin]} re-runs {cmd:e(cmdline)} in modified form
(options listed in {cmd:remove()} are deleted including nested parentheses; {cmd:remove(all)} deletes all options;
{cmd:nowt} and {cmd:noifin} drop the weight and the if/in qualifiers) and returns {cmd:r(rc)} and {cmd:r(cmdline)};
call {cmd:_atest_restore} afterwards. If the command does not store {cmd:e(cmdline)}, the command line is rebuilt from
{cmd:e(cmd)}, {cmd:e(depvar)}, the terms of {cmd:e(b)}, the weights and the estimation sample.{p_end}
{phang2}{cmd:_atest_linear, contlist() scale() [lr]}, {cmd:_atest_splinfo}, {cmd:_atest_gb}, {cmd:_atest_hl, pr() y()}, {cmd:_atest_infl [, note()]},
{cmd:_atest_phtest [, note()]} and {cmd:_atest_wald, par() [null()]}
are reusable test building blocks.{p_end}

{pstd}
Globals available to modules: {cmd:$ATEST_esamp} (0/1 estimation-sample variable),
{cmd:$ATEST_cont}, {cmd:$ATEST_cat}, {cmd:$ATEST_spl} (term lists), {cmd:$ATEST_nsp} and
{cmd:$ATEST_sp1}, {cmd:$ATEST_sp2}, ... (spline groups), {cmd:$ATEST_alpha},
{cmd:$ATEST_lin}, {cmd:$ATEST_infl}, {cmd:$ATEST_outcut}, {cmd:$ATEST_ldcut}, {cmd:$ATEST_cookcut},
{cmd:$ATEST_rescut}, {cmd:$ATEST_vifcut}, {cmd:$ATEST_epvcut}, {cmd:$ATEST_nooutliers}, {cmd:$ATEST_graphpath}.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:assumptiontest} does not clear {cmd:r()}. The {cmd:r()} results that existed before
the call, such as {cmd:r(table)} and {cmd:r(level)} from the estimation output, are kept
unchanged. By default the results below are added to them; if an earlier result has the
same name as one of these, the {cmd:assumptiontest} result replaces it. With
{cmd:noreturn}, nothing is added and {cmd:r()} is exactly as before the call.
{cmd:r()} is set in the same way when the command stops with error 459.

{pstd}
{cmd:assumptiontest} adds to {cmd:r()} (unless {cmd:noreturn} is specified):

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}estimation sample size{p_end}
{synopt:{cmd:r(alpha)}}significance level used{p_end}
{synopt:{cmd:r(n_tests)}}number of rows posted{p_end}
{synopt:{cmd:r(n_violated)}}critical tests rejected{p_end}
{synopt:{cmd:r(n_flagged)}}non-critical tests rejected{p_end}
{synopt:{cmd:r(n_untestable)}}assumptions not testable on the data{p_end}
{synopt:{cmd:r(n_notrun)}}tests that could not be run{p_end}
{synopt:{cmd:r(passed)}}1 if no critical test was rejected{p_end}
{synopt:{cmd:r(p_}{it:name}{cmd:)}}p-value of each test that has one{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}, {cmd:r(model)}}estimation command and model title{p_end}
{synopt:{cmd:r(tests)}}names of all rows{p_end}
{synopt:{cmd:r(violated)}, {cmd:r(flagged)}}names of rejected critical / non-critical tests{p_end}
{synopt:{cmd:r(untestable)}, {cmd:r(notrun)}}names of untestable assumptions / tests not run{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:r(results)}}one row per test: critical, code, stat, df, p (code 0 passed, 1 rejected, 2 not testable, 3 not run, 4 information){p_end}

{pstd}
The same table is also kept in matrix {cmd:atest_results} (also with {cmd:noreturn}
and after the error stop); it is replaced by the next call.

{pstd}
Errors that occur before the tests run (for example error 301, no estimation results,
or error 198, invalid option) leave {cmd:r()} unchanged.


{title:Examples}

{phang2}{cmd:. stcox age i.sex i.grp}{p_end}
{phang2}{cmd:. assumptiontest, graphpath("C:/results/graphs")}{p_end}

{phang2}{cmd:. mkspline agesp = age, cubic nknots(4)}{p_end}
{phang2}{cmd:. stcox agesp1 agesp2 agesp3 i.sex}{p_end}
{phang2}{cmd:. assumptiontest, spline(agesp1 agesp2 agesp3) continue}{p_end}

{phang2}{cmd:. nbreg visits age i.sex, exposure(years)}{p_end}
{phang2}{cmd:. assumptiontest, noncritical(lin_age) cookcut(0.01)}{p_end}

{phang2}{cmd:. poisson visits age i.sex, exposure(years) vce(robust)}{p_end}
{phang2}{cmd:. assumptiontest}{p_end}

{phang2}{cmd:. regress sbp age bmi i.sex}{p_end}
{phang2}{cmd:. assumptiontest, critical(normal) vifcut(5)}{p_end}

{phang2}{cmd:. logit death age i.nyha}{p_end}
{phang2}{cmd:. assumptiontest, epvcut(20)}{p_end}
{phang2}{cmd:. firthlogit death age i.nyha}{p_end}
{phang2}{cmd:. assumptiontest}{p_end}

{phang2}{cmd:. assumptiontest, linearity(quadratic) influence(dfbeta) outcut(0.5)}{p_end}

{pstd}
Keep {cmd:r()} from the estimation untouched (for example to export {cmd:r(table)} afterwards):

{phang2}{cmd:. stcox age i.sex}{p_end}
{phang2}{cmd:. assumptiontest, noreturn continue}{p_end}
{phang2}{cmd:. matrix list r(table)}{p_end}


{title:Remarks}

{pstd}
Several tests are approximations (Gronnesby-Borgan type calibration, Hosmer-Lemeshow, the NB2
variance-function score test, Pearson chi-squared) and are therefore non-critical by default. With many covariates and no
multiple-testing correction, occasional rejections are expected by chance; use
{cmd:critical()}/{cmd:noncritical()} to adapt the classification.

{pstd}
Models fit with an {cmd:svy:} prefix are not covered.
