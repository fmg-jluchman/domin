{smcl}
{* *! version 2.1.1 December 20, 2024 J. N. Luchman}{...}
{cmd:help mixdom}
{hline}{...}

{title:Title}

{pstd}
Wrapper program for {cmd:domin} to conduct linear mixed effects regression dominance analysis{p_end}

{title:Syntax}

{phang}
{cmd:mixdom} {it:depvar} {it:indepvars} {it:{help if}} {weight} {cmd:,} 
{opt id(idvar)} [{opt {ul on}re{ul off}opt(re_options)} {opt {ul on}m{ul off}opt(mixed_options)}]

{phang}{cmd:pweight}s and {cmd:fweight}s are allowed (see help {help weights:weights}).  {help fvvarlist:Factor} and 
{help tsvarlist:time series variables} are allowed.  

{title:Description}

{pstd}
{cmd:mixdom} sets the data up in a way to allow for the dominance analysis of a linear mixed effects regression by utilizing {help mixed}.
The method outlined here follows that for the within- and between-cluster Snijders and Bosker (1994) R2 metric described by Luo and Azen (2013). 

{pstd}
{cmd:mixdom} only allows 1 level of clustering in the data (i.e., 1 random effect), which must be the cluster constant/mean/intercept. Luo and 
Azen (2013) recommend that even if random coefficients are present in the data, they should be restricted to a fixed effect only in the dominance 
analysis.

{pstd}
{cmd:mixdom} is intended for use only as a wrapper program with {cmd:domin} for the dominance analysis of mixed-effects linear regression, and its syntax is designed to conform with {cmd:domin}'s expectations.  
It is not recommended for use as an estimation command outside of {cmd:domin}.

{pstd}
Note that negative R2 values indicate likely model misspecification.

{marker options}{...}
{title:Options}

{phang}{opt id()} specifies the variable on which clustering occurs and that will appear after the random effects specification (i.e., ||) in the 
{cmd:mixed} syntax.

{phang}{opt reopt()} passes options to {cmd: mixed} specific to the random intercept effect (i.e., {opt pweight()} the user would 
like to utilize during estimation.

{phang}{opt mopt()} passes options to {cmd:mixed} that the user would like to utilize during estimation.
This option was named {opt xtmopt()} in {cmd:mixdom} versions previous to 2.0.0 and is now defunct.

{title:Saved results}

{phang}{cmd:mixdom} saves the following results to {cmd: e()}:

{synoptset 16 tabbed}{...}
{p2col 5 15 19 2: scalars}{p_end}
{synopt:{cmd:e(r2_w)}}within-cluster R2{p_end}
{synopt:{cmd:e(r2_b)}}between-cluster R2{p_end}
{p2col 5 15 19 2: macros}{p_end}
{synopt:{cmd:e(title)}}"Mixed-effects ML regression"{p_end}
{p2col 5 15 19 2: functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{title:References}

{p 4 8 2}Luo, W., & Azen, R. (2013). Determining predictor importance in hierarchical linear models using dominance analysis. {it:Journal of Educational and Behavioral Statistics, 38(1)}, 3-31.{p_end}
{p 4 8 2}Snijders, T. A. B., & Bosker, R. J. (1994). Modeled variance in two-level models. {it:Sociological Methods & Research, 22(3)}, 342-363.{p_end}

{title:Author}

{p 4}Joseph N. Luchman{p_end}
{p 4}Research Fellow{p_end}
{p 4}Fors Marsh{p_end}
{p 4}Arlington, VA{p_end}
{p 4}jluchman@forsmarsh.com{p_end}

{title:Also see}

{psee}
{manhelp mixed R}.
{p_end}
