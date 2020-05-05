# Dominance Analysis
## A Stata Implementaion

Dominance analysis (DA) determines the relative importance of independent variables in an estimation model based on contribution to an overall model fit statistic (see Grömping, 2007 for a discussion).  DA is an ensemble method in which importance determinations about independent variables are made by aggregating results across multiple models, though the method usually requires the ensemble contain each possible combination of the independent variables in the full model.  

The all possible combinations ensemble with p independent variables in the full model results in 2^p-1 models estimated.  That is, each combiation of p variables alterating between included versus excluded (i.e., the 2 base to the exponent) where the constant[s]-only model is omitted (i.e., the -1 representing the distinct combination where no independent variables are included; see Budescu, 1993).

`domin` is implemented as a flexible wrapper command that can be used with most Stata estimation commands that follow the standard `depvar indepvars` format and return a scalar-valued fit metric; commands that do not can be accommodated with a sub-wrapper command (an example of this is included below).

Some examples of the command as applied to Stata estimation commands are shown below.

# Initial Examples

## Simple Linear Regression-based DA

The default analysis for `domin` is `regress` with `fitstat(e(r2))` and these do not need to be typed (though if they are not, `domin` will throw a warning).  The results of the analysis are shown as below including general, conditional, and complete dominance results as well as the strongest dominance designations.

```
.    webuse auto
(1978 Automobile Data)

.     domin price mpg rep78 headroom
Regression type not entered in reg(). 
reg(regress) assumed.

Fitstat type not entered in fitstat(). 
fitstat(e(r2)) assumed.


Total of 7 regressions

General dominance statistics: Linear regression
Number of obs             =                      69
Overall Fit Statistic     =                  0.2575

            |      Dominance      Standardized      Ranking
 price      |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 mpg        |         0.2262      0.8787            1 
 rep78      |         0.0218      0.0847            2 
 headroom   |         0.0094      0.0366            3 
-------------------------------------------------------------------------------------
-more-
```

### General Dominance Statistics

General dominance statistics are the most commonly reported and easiest to interpret. General dominance statistics are derived as a weighted average marginal/incremental contribution to the overall fit statistic an independent variable makes across all models in which the independent variable is included.  For example, _rep78_ has a larger general dominance statistic than, and thus "generally dominates", independent variable _headroom_.  If general dominance statistics are equal for two independent variables, no general dominance designation can be made between those independent variables.

General dominance statistics distill the entire ensemble of models into a single value for each independent variable, which is why they are easiest to interpret. In addition, a useful property of the general dominance statistics is that they are an additive decomposition of the fit statistic associated with the full model (i.e., the general dominance statistics can be summed to obtain the value of the full model's fit statistic).  Thus, general dominance statistics are equivalent to Shapley values (see `findit shapley`).  General dominance statistics are the arithmetic average of all conditional dominance statistics discussed next.


```
Conditional dominance statistics
-------------------------------------------------------------------------------------

           #indepvars:  #indepvars:  #indepvars:
                    1            2            3
     mpg       0.2079       0.2262       0.2445
   rep78       0.0000       0.0218       0.0436
headroom       0.0124       0.0094       0.0065
-------------------------------------------------------------------------------------
-more-
```

### Conditional Dominance Statistics

Conditional dominance statistics are also derived from the all possible combinations ensemble.  Conditional dominance statistics are computed as the average incremental contributions to the overall model fit statistic within a single "order" for models in which the independent variable is included - where "order" refers to a distinct number of independent variables in the estimation model.  One order is thus all models that include one independent variable.  Another order is all models that include two independent variables, and so on to p - or the order including only the model with all p independent variables.  Each independent variable will then have p different conditional dominance statistics.  In the example above, there are three conditional dominance statistics for each independent variable because there are three independent variables

The evidence conditional dominance statistics provide with respect to relative importance is stronger than that provided by general dominance statistics.  Because general dominance statistics are the arithmetic average of all p conditional dominance statistics, conditional dominance statistics, considered as a set, provide more information about each independent variable or, alternatively, are less "averaged" than general dominance statistics. Conditional dominance statistics also provide information about independent variable redundancy, collinearity, and suppression effects as the user can see how the inclusion of any independent variable is, on average, affected by the inclusion of other independent variables in the estimation model in terms of their effect on model fit.  In the above conditional dominance matrix, observe the difference between the patterns of results for _rep78_ and _headroom_.  _rep78_ shows stronger suppression-like effects in that it grows in predictive usefulness with more independent variables.  By contrast, _headroom_ shows the opposite pattern shrinking in importance.

For example, _mpg_ has larger conditional dominance statistics than independent variable _rep78_ across all three orders and thus "conditionally dominates".  To be more specific, for an independent variable conditionally dominate another, its conditional dominancer statistic must be larger than the other across all p orders.  If, at any order, the conditional dominance statistics for two independent variables are equal or there is a change rank order no conditional dominance designation can be made between those independent variables.  Conditional dominance imples general dominance as well, but the reverse is not true.  An independent variable can generally dominate another, but not conditionally dominate it.  For instance, _rep78_ generally dominance but does not conditionally dominate _headroom_.

```
Complete dominance designation
-------------------------------------------------------------------------------------

                      dominated?:  dominated?:  dominated?:
                             mpg        rep78     headroom
     dominates?:mpg            0            1            1
   dominates?:rep78           -1            0            0
dominates?:headroom           -1            0            0
-------------------------------------------------------------------------------------
-more-
```

### Complete Dominance Designations

Complete dominance designations are the final designation derived from the all possible combinations ensemble.  Complete dominance designations are made by comparing all possible incremental contributions to model fit for two independent variables. The evidence the complete dominance designation provides with respect to relative importance is the strongest possible, and supercedes that of general and conditional dominance.  Complete dominance is the strongest evidence as it is completely un-averaged and pits each independent variable against one another in every possible comparison.  Thus, it is not possible for some good incremental contributions to compensate for some poorer incremental contributions as can occur when such data are averaged.  Complete dominance then provides information on a property of the entire ensemble of models, as it relates to a comparison between two independent variables.

For example, _mpg_ has a larger incremental contribution to model fit than _headroom_ across all possible comparisons and "completely dominates" it. As with conditional dominance designations, for an independent variable to completely dominate another, the incremental contribution to fit associated with that independent variable for each of the possible 2^(p-2) comparisons with another must all be larger than another.  If, for any comparison, the incremental contribution to fit for two independent variables are equal or there is a change in rank order, no complete dominance designation can be made between those independent variables.  Complete dominance imples both general and conditional dominance, but, again, the reverse is not true. By comparison to general and conditional dominance designations, the complete dominance designation has no natural statistic.  That said, domin returns a complete dominance matrix which reads from the left to right.  Thus, a value of 1 means that the indepdendent variable in the row completely dominates the independent variable in the column.  Conversely, a value of -1 means the opposite, that the independent variable in the row is completely dominated by the independent variable in the column.  A 0 value means no complete dominance designation could be made as the comparison independent variables' incremental contributions differ in relative magnitude from model to model.

```
Strongest dominance designations

mpg completely dominates rep78
mpg completely dominates headroom
rep78 generally dominates headroom
```

### Strongest Dominance Designtions

Finally, if all three dominance statistics are reported (i.e., `noconditional` and `nocomplete` options are not used), a "strongest dominance designations" list is reported.  The strongest dominance designations list reports the strongest dominance designation between all pairwise, independent variable comparisons.


## Ordered Logistic Regression-based DA

A model like `ologit` is easy to accommodate in `domin` like below.

```
. domin rep78 trunk weight length, reg(ologit) fitstat(e(r2_p)) all(turn)

Total of 7 regressions

General dominance statistics: Ordered logistic regression
Number of obs             =                      69
Overall Fit Statistic     =                  0.1209
All Subsets Fit Stat.     =                  0.1003

            |      Dominance      Standardized      Ranking
 rep78      |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 trunk      |         0.0082      0.0680            2 
 weight     |         0.0021      0.0174            3 
 length     |         0.0102      0.0845            1 
-------------------------------------------------------------------------------------
Conditional dominance statistics
-------------------------------------------------------------------------------------

         #indepvars:  #indepvars:  #indepvars:
                  1            2            3
 trunk       0.0131       0.0077       0.0039
weight       0.0027       0.0016       0.0020
length       0.0139       0.0097       0.0070
-------------------------------------------------------------------------------------
Complete dominance designation
-------------------------------------------------------------------------------------

                    dominated?:  dominated?:  dominated?:
                         trunk       weight       length
 dominates?:trunk            0            1           -1
dominates?:weight           -1            0           -1
dominates?:length            1            1            0
-------------------------------------------------------------------------------------

Strongest dominance designations

length completely dominates trunk
trunk completely dominates weight
length completely dominates weight

Variables included in all subsets: turn
```

As compared to the `regress`-based DA reported in the first example, this example includes a covariate that is controlled for across all model subsets.  The `All Subsets Fit Stat.     =                  0.1003` result represents the amount of the McFadden pseudo-R-square that is associated with the variable in all subsets (i.e., in `all()`).  Note that the dominance statistics reported are now residualized and reflect the removal of the all fitstats fit statistic.  Variables included in all subsets are also reported at the end of the results display (e.g., `Variables included in all subsets: turn`).

## Logistic Regression-based DA

`logit` is another command that can be accommodated in `domin`.  In this example, _rep78_ is used as a factor variable.

```
. domin foreign trunk weight, reg(logit) fitstat(e(r2_p)) sets((i.rep78))

Total of 7 regressions

General dominance statistics: Logistic regression
Number of obs             =                      59
Overall Fit Statistic     =                  0.5229

            |      Dominance      Standardized      Ranking
 foreign    |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 trunk      |         0.0834      0.1595            3 
 weight     |         0.2611      0.4994            1 
 set1       |         0.1783      0.3411            2 
-------------------------------------------------------------------------------------
Conditional dominance statistics
-------------------------------------------------------------------------------------

         #indepvars:  #indepvars:  #indepvars:
                  1            2            3
 trunk       0.2002       0.0495       0.0005
weight       0.4154       0.2272       0.1407
  set1       0.2855       0.1444       0.1051
-------------------------------------------------------------------------------------
Complete dominance designation
-------------------------------------------------------------------------------------

                    dominated?:  dominated?:  dominated?:
                         trunk       weight         set1
 dominates?:trunk            0           -1           -1
dominates?:weight            1            0            1
  dominates?:set1            1           -1            0
-------------------------------------------------------------------------------------

Strongest dominance designations

weight completely dominates trunk
set1 completely dominates trunk
weight completely dominates set1

Variables in set1: i.rep78
```

The `sets()` option incorporates independent variables as inseparable sets that are considered *an* independent variable in the DA.  The label for each set is denoted sequentially as they are included (e.g., _set1_).  In addition, variables in sets are included near the end of the output (i.e., `Variables in set1: i.rep78`). 

# Intermediate Examples

## Linear Regression-based DA with Non-additive and Non-linear Effects

The example below outlines a process in which two quadratic effects and a product are residualized and used as independent variables in a linear regression-based DA.

```
.     generate mpg2 = mpg^2

. 
.     generate headr2 = headroom^2

. 
.     generate mpg_headr = mpg*headroom

. 
.     regress mpg2 mpg

      Source |       SS           df       MS      Number of obs   =        74
-------------+----------------------------------   F(1, 72)        =   2425.89
       Model |  5640570.83         1  5640570.83   Prob > F        =    0.0000
    Residual |  167411.003        72  2325.15282   R-squared       =    0.9712
-------------+----------------------------------   Adj R-squared   =    0.9708
       Total |  5807981.84        73   79561.395   Root MSE        =     48.22

------------------------------------------------------------------------------
        mpg2 |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
         mpg |   48.04619   .9754908    49.25   0.000     46.10159     49.9908
       _cons |  -536.6594   21.51824   -24.94   0.000    -579.5552   -493.7636
------------------------------------------------------------------------------

. 
.     predict mpg2r, resid

. 
.     regress headr2 headroom

      Source |       SS           df       MS      Number of obs   =        74
-------------+----------------------------------   F(1, 72)        =   3197.55
       Model |  1947.18924         1  1947.18924   Prob > F        =    0.0000
    Residual |  43.8453928        72  .608963789   R-squared       =    0.9780
-------------+----------------------------------   Adj R-squared   =    0.9777
       Total |  1991.03463        73   27.274447   Root MSE        =    .78036

------------------------------------------------------------------------------
      headr2 |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
    headroom |    6.10485   .1079609    56.55   0.000     5.889633    6.320066
       _cons |  -8.607759   .3356446   -25.65   0.000    -9.276855   -7.938664
------------------------------------------------------------------------------

. 
.     predict headr2r, resid

. 
.     regress mpg_headr mpg headroom

      Source |       SS           df       MS      Number of obs   =        74
-------------+----------------------------------   F(2, 71)        =   1102.79
       Model |  24303.2769         2  12151.6384   Prob > F        =    0.0000
    Residual |  782.348117        71  11.0189876   R-squared       =    0.9688
-------------+----------------------------------   Adj R-squared   =    0.9679
       Total |   25085.625        73  343.638699   Root MSE        =    3.3195

------------------------------------------------------------------------------
   mpg_headr |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
         mpg |   2.837136   .0737653    38.46   0.000     2.690052     2.98422
    headroom |   20.40514   .5044586    40.45   0.000     19.39927      21.411
       _cons |  -59.75086   2.619192   -22.81   0.000    -64.97338   -54.52834
------------------------------------------------------------------------------

. 
.     predict mpg_headrr, resid

. 
.     domin price mpg headroom mpg2r headr2r mpg_headrr
Regression type not entered in reg(). 
reg(regress) assumed.

Fitstat type not entered in fitstat(). 
fitstat(e(r2)) assumed.


Total of 31 regressions

Progress in running all regression subsets
0%------50%------100%
....................
General dominance statistics: Linear regression
Number of obs             =                      74
Overall Fit Statistic     =                  0.3948

            |      Dominance      Standardized      Ranking
 price      |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 mpg        |         0.2274      0.5761            1 
 headroom   |         0.0142      0.0361            4 
 mpg2r      |         0.1044      0.2644            2 
 headr2r    |         0.0388      0.0984            3 
 mpg_headrr |         0.0099      0.0251            5 
-------------------------------------------------------------------------------------
Conditional dominance statistics
-------------------------------------------------------------------------------------

             #indepvars:  #indepvars:  #indepvars:  #indepvars:  #indepvars:
                      1            2            3            4            5
       mpg       0.2196       0.2191       0.2236       0.2323       0.2425
  headroom       0.0131       0.0091       0.0101       0.0156       0.0233
     mpg2r       0.1203       0.1104       0.1029       0.0973       0.0910
   headr2r       0.0429       0.0435       0.0419       0.0376       0.0282
mpg_headrr       0.0072       0.0109       0.0124       0.0118       0.0072
-------------------------------------------------------------------------------------
Complete dominance designation
-------------------------------------------------------------------------------------

                        dominated?:  dominated?:  dominated?:  dominated?:  dominated?:
                               mpg     headroom        mpg2r      headr2r   mpg_headrr
       dominates?:mpg            0            1            1            1            1
  dominates?:headroom           -1            0           -1           -1            0
     dominates?:mpg2r           -1            1            0            1            1
   dominates?:headr2r           -1            1           -1            0            1
dominates?:mpg_headrr           -1            0           -1           -1            0
-------------------------------------------------------------------------------------

Strongest dominance designations

mpg completely dominates headroom
mpg2r completely dominates headroom
headr2r completely dominates headroom
mpg completely dominates mpg2r
mpg completely dominates headr2r
mpg2r completely dominates headr2r
mpg completely dominates mpg_headrr
mpg2r completely dominates mpg_headrr
headr2r completely dominates mpg_headrr
headroom generally dominates mpg_headrr
```

## Linear Regression-based Relative Weights Analysis with Bootstrapped Standard Errors

`domin` can also estimate relative weights analysis (RWA) using the `epsilon` option.  RWA is a faster, approximation to DA and obviates the each subset regression by orthogonalizing independent variables using singular value decomposition (see `matrix svd`).  `epsilon`'s singular value decomposition approach is not equivalent to the all possible combinations ensemble approach but is many fold faster for models with many independent variables and tends to produce similar answers regarding relative importance.  

DA and RWA can also be `boostrap`-ped to obtain stadard errors and evaluate stability if desired.

```
.     bootstrap, reps(500): domin price mpg headroom trunk turn gear_ratio foreign length weight, epsilon
(running domin on estimation sample)

Bootstrap replications (500)
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 
..................................................    50
..................................................   100
..................................................   150
..................................................   200
..................................................   250
..................................................   300
..................................................   350
..................................................   400
..................................................   450
..................................................   500

General dominance statistics: Custom user analysis
Number of obs             =                      74
Overall Fit Statistic     =                  0.5806

            |      Dominance      Standardized      Ranking
 price      |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 mpg        |         0.0803      0.1382            4 
 headroom   |         0.0148      0.0255            8 
 trunk      |         0.0349      0.0601            7 
 turn       |         0.0533      0.0918            6 
 gear_ratio |         0.0619      0.1067            5 
 foreign    |         0.0910      0.1567            2 
 length     |         0.0824      0.1420            3 
 weight     |         0.1620      0.2791            1 
-------------------------------------------------------------------------------------

. 
.     estat bootstrap

Dominance analysis                              Number of obs     =         74
                                                Replications      =        500

------------------------------------------------------------------------------
             |    Observed               Bootstrap
       price |       Coef.       Bias    Std. Err.  [95% Conf. Interval]
-------------+----------------------------------------------------------------
         mpg |   .08025451   .0056983   .02608266    .0428729   .1343196  (BC)
    headroom |   .01480511    .005788   .00994976    .0045644   .0253761  (BC)
       trunk |   .03488086   .0052553   .01569048    .0146097   .0700839  (BC)
        turn |   .05328024   .0088603   .01346109    .0228937   .0711584  (BC)
  gear_ratio |   .06193881   .0102595   .02574597    .0342715   .1201322  (BC)
     foreign |   .09099184   .0039005   .03310621    .0428155   .1721853  (BC)
      length |   .08242763   .0060164   .01522731    .0471752   .1077507  (BC)
      weight |   .16203422   .0023967   .03416152    .0966866   .2320645  (BC)
------------------------------------------------------------------------------
(BC)   bias-corrected confidence interval
```

## Multivariate Linear Regression-based DA

`domin` comes packaged with a sub-wrapper program `mvdom` that allows for estimating DA statistics with multiple dependent variables.

The dependent variables in the below model are _price_, *gear_ratio*, _foreign_, _length_, and _weight_. 


```
. domin price mpg headroom trunk turn, reg(mvdom, dvs(gear_ratio foreign length weight)) fitstat(e(r2))

Total of 15 regressions

General dominance statistics: Multivariate regression
Number of obs             =                      74
Overall Fit Statistic     =                  0.8659

            |      Dominance      Standardized      Ranking
 price      |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 mpg        |         0.2605      0.3008            2 
 headroom   |         0.0794      0.0917            4 
 trunk      |         0.1779      0.2055            3 
 turn       |         0.3482      0.4020            1 
-------------------------------------------------------------------------------------
Conditional dominance statistics
-------------------------------------------------------------------------------------

           #indepvars:  #indepvars:  #indepvars:  #indepvars:
                    1            2            3            4
     mpg       0.6798       0.2314       0.0987       0.0319
headroom       0.2923       0.0191       0.0052       0.0009
   trunk       0.5379       0.1251       0.0339       0.0148
    turn       0.7887       0.3239       0.1782       0.1018
-------------------------------------------------------------------------------------
Complete dominance designation
-------------------------------------------------------------------------------------

                      dominated?:  dominated?:  dominated?:  dominated?:
                             mpg     headroom        trunk         turn
     dominates?:mpg            0            1            1           -1
dominates?:headroom           -1            0           -1           -1
   dominates?:trunk           -1            1            0           -1
    dominates?:turn            1            1            1            0
-------------------------------------------------------------------------------------

Strongest dominance designations

turn completely dominates mpg
mpg completely dominates headroom
trunk completely dominates headroom
turn completely dominates headroom
mpg completely dominates trunk
turn completely dominates trunk
```

## Gamma regression-based DA

`domin` allows for any scalar valued fit metric to br used as a DA-able metric.

In this example, the model deviance is used as a fit statistic.

```
. domin price mpg rep78 headroom, reg(glm, family(gamma) link(power -1)) fitstat(e(deviance)) consmodel reverse

Total of 7 regressions

General dominance statistics: Generalized linear models
Number of obs             =                      69
Overall Fit Statistic     =                  7.4580
Constant-only Fit Stat.   =                 11.6141

            |      Dominance      Standardized      Ranking
 price      |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 mpg        |        -3.7935      0.9127            1 
 rep78      |        -0.2118      0.0510            2 
 headroom   |        -0.1509      0.0363            3 
-------------------------------------------------------------------------------------
Conditional dominance statistics
-------------------------------------------------------------------------------------

           #indepvars:  #indepvars:  #indepvars:
                    1            2            3
     mpg      -3.6135      -3.8087      -3.9583
   rep78      -0.0007      -0.2270      -0.4077
headroom      -0.1876      -0.1661      -0.0989
-------------------------------------------------------------------------------------
Complete dominance designation
-------------------------------------------------------------------------------------

                      dominated?:  dominated?:  dominated?:
                             mpg        rep78     headroom
     dominates?:mpg            0            1            1
   dominates?:rep78           -1            0            0
dominates?:headroom           -1            0            0
-------------------------------------------------------------------------------------

Strongest dominance designations

mpg completely dominates rep78
mpg completely dominates headroom
rep78 generally dominates headroom
```

Important notes about this model is that a `consmodel` is estimated (e.g., `Constant-only Fit Stat.   =                 11.6141`) that is, like entries in the `all()` model, subtracted from the DA statistics.  This differs from `all()` in that it assumes only estimation model *_cons*-tants/intercepts are estimated. This DA also `reverse`s the interpretation of the fit metrics such that smaller values represent more important independent variables (i.e., larger negative values).

## Linear Mixed-effects Regression-based DA

`domin` also comes packaged with a sub-wrapper program `mixdom` that allows for estimating DA statistics with linear mixed effects models.

There can only be two levels of clustering/nesting in the data.  The analysis below uses the within-cluster r-square. 

```
.     webuse nlswork, clear
(National Longitudinal Survey.  Young Women 14-26 years of age in 1968)

. 
.     domin ln_wage tenure hours age collgrad, reg(mixdom, id(id)) fitstat(e(r2_w)) sets((i.race))

Total of 31 regressions

Progress in running all regression subsets
0%------50%------100%
....................
General dominance statistics: Mixed-effects ML regression
Number of obs             =                   28036
Overall Fit Statistic     =                  0.2688

            |      Dominance      Standardized      Ranking
 ln_wage    |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 tenure     |         0.0958      0.3564            2 
 hours      |         0.0031      0.0114            5 
 age        |         0.0406      0.1510            3 
 collgrad   |         0.1175      0.4371            1 
 set1       |         0.0119      0.0441            4 
-------------------------------------------------------------------------------------
Conditional dominance statistics
-------------------------------------------------------------------------------------

           #indepvars:  #indepvars:  #indepvars:  #indepvars:  #indepvars:
                    1            2            3            4            5
  tenure       0.1389       0.1130       0.0916       0.0744       0.0611
   hours       0.0053       0.0040       0.0029       0.0020       0.0012
     age       0.0860       0.0590       0.0364       0.0180       0.0034
collgrad       0.1461       0.1274       0.1132       0.1034       0.0975
    set1       0.0157       0.0135       0.0117       0.0100       0.0084
-------------------------------------------------------------------------------------
Complete dominance designation
-------------------------------------------------------------------------------------

                      dominated?:  dominated?:  dominated?:  dominated?:  dominated?:
                          tenure        hours          age     collgrad         set1
  dominates?:tenure            0            1            1            0            1
   dominates?:hours           -1            0           -1           -1           -1
     dominates?:age           -1            1            0           -1            0
dominates?:collgrad            0            1            1            0            1
    dominates?:set1           -1            1            0           -1            0
-------------------------------------------------------------------------------------

Strongest dominance designations

tenure completely dominates hours
age completely dominates hours
collgrad completely dominates hours
set1 completely dominates hours
tenure completely dominates age
collgrad completely dominates age
tenure completely dominates set1
collgrad completely dominates set1
collgrad conditionally dominates tenure
age generally dominates set1

Variables in set1: i.race
```

# Advanced Examples

## Multinomial Logistic Regression-based DA
### Custom Program Generation using BIC as Fit Statistic

`domin` can be adapted to dominance analyse many programs--so long as they adhere to a specific structure `domin` needs to properly evaluate the model.

In the example below, a custom program _myprog_ is defined to obtain the Bayesian Information Criterion following a `mlogit`.

```
.     program define myprog, eclass
  1. 
.     syntax varlist if , [option]
  2. 
.     tempname estlist
  3. 
.     mlogit `varlist' `if'
  4. 
.     estat ic
  5. 
.     matrix `estlist' = r(S)
  6. 
.     ereturn scalar bic = `estlist'[1,6]
  7. 
.     end

. 
.     domin race tenure hours age nev_mar, reg(myprog) fitstat(e(bic)) consmodel reverse

Total of 15 regressions

General dominance statistics: Multinomial logistic regression
Number of obs             =                   28022
Overall Fit Statistic     =              36025.3867
Constant-only Fit Stat.   =              36482.9057

            |      Dominance      Standardized      Ranking
 race       |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 tenure     |        13.2083     -0.0289            4 
 hours      |       -84.8225      0.1854            2 
 age        |         3.0220     -0.0066            3 
 nev_mar    |      -388.9268      0.8501            1 
-------------------------------------------------------------------------------------
Conditional dominance statistics
-------------------------------------------------------------------------------------

          #indepvars:  #indepvars:  #indepvars:  #indepvars:
                   1            2            3            4
 tenure      15.4426      10.2512      12.7821      14.3572
  hours    -102.4587     -91.2406     -77.0439     -68.5468
    age       9.5965       2.0195       1.4031      -0.9312
nev_mar    -399.0177    -395.1525    -384.8246    -376.7124
-------------------------------------------------------------------------------------
Complete dominance designation
-------------------------------------------------------------------------------------

                     dominated?:  dominated?:  dominated?:  dominated?:
                         tenure        hours          age      nev_mar
 dominates?:tenure            0           -1           -1           -1
  dominates?:hours            1            0            1           -1
    dominates?:age            1           -1            0           -1
dominates?:nev_mar            1            1            1            0
-------------------------------------------------------------------------------------

Strongest dominance designations

hours completely dominates tenure
age completely dominates tenure
nev_mar completely dominates tenure
nev_mar completely dominates hours
hours completely dominates age
nev_mar completely dominates age
```

## Logistic Regression-based DA: Revisited
### Multiply Imputed Estimates

`domin` can also accommodate multiply imputed data.

```
.     webuse mheart1s20, clear
(Fictional heart attack data; bmi missing)

. 
.     domin attack smokes age bmi hsgrad female, reg(logit) fitstat(e(r2_p)) mi

Total of 31 regressions

Progress in running all regression subsets
0%------50%------100%
....................
General dominance statistics: Logistic regression
Number of obs             =                     154
Overall Fit Statistic     =                  0.1056

            |      Dominance      Standardized      Ranking
 attack     |      Stat.          Domin. Stat.
------------+------------------------------------------------------------------------
 smokes     |         0.0553      0.5234            1 
 age        |         0.0248      0.2354            3 
 bmi        |         0.0250      0.2369            2 
 hsgrad     |         0.0003      0.0028            4 
 female     |         0.0002      0.0016            5 
-------------------------------------------------------------------------------------
Conditional dominance statistics
-------------------------------------------------------------------------------------

         #indepvars:  #indepvars:  #indepvars:  #indepvars:  #indepvars:
                  1            2            3            4            5
smokes       0.0548       0.0549       0.0551       0.0554       0.0560
   age       0.0231       0.0238       0.0247       0.0257       0.0270
   bmi       0.0222       0.0235       0.0249       0.0264       0.0281
hsgrad       0.0000       0.0001       0.0002       0.0004       0.0008
female       0.0000       0.0001       0.0002       0.0002       0.0003
-------------------------------------------------------------------------------------
Complete dominance designation
-------------------------------------------------------------------------------------

                    dominated?:  dominated?:  dominated?:  dominated?:  dominated?:
                        smokes          age          bmi       hsgrad       female
dominates?:smokes            0            1            1            1            1
   dominates?:age           -1            0            0            1            1
   dominates?:bmi           -1            0            0            1            1
dominates?:hsgrad           -1           -1           -1            0            0
dominates?:female           -1           -1           -1            0            0
-------------------------------------------------------------------------------------

Strongest dominance designations

smokes completely dominates age
smokes completely dominates bmi
smokes completely dominates hsgrad
age completely dominates hsgrad
bmi completely dominates hsgrad
smokes completely dominates female
age completely dominates female
bmi completely dominates female
bmi generally dominates age
hsgrad generally dominates female
```