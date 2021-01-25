from sys import argv as use_stata_arguments
import itertools as it
import sfi
import statistics as stat
from math import factorial as fctl

    # ~~ Process arguments passed from stata ~~ #

Stata_Regression = use_stata_arguments[1]
Dep_Var = use_stata_arguments[2]
Indep_Vars_Unprocessed = use_stata_arguments[3]
AllSubsets_Indep_Vars = use_stata_arguments[4]
If_Conditions = use_stata_arguments[5]
Regress_Options = use_stata_arguments[6]
Fit_Statistic = use_stata_arguments[7]
Mult_Impute_Flag = use_stata_arguments[8]
Mult_Impute_Opts = use_stata_arguments[9]
FitStat_Adjustment = float(use_stata_arguments[10])
Conditional_Flag = not bool(len(use_stata_arguments[11]))
Complete_Flag = not bool(len(use_stata_arguments[12]))
Mult_Impute_File = use_stata_arguments[13]
Mult_Imputes_toUse = use_stata_arguments[14]


    # ~~ Create independent variable list ~~ #
    
Indep_Var_List = Indep_Vars_Unprocessed.strip().split("<") # parses IV sets, if any (individual IVs unprocessed)

if len(Indep_Var_List[0]) > 0:
    Indep_Var_List = ( Indep_Var_List[0].strip().split(" ") + 
        Indep_Var_List[1:len(Indep_Var_List)] ) # if there are only individual IVs, or if there are invidual IVs and IV sets, parse individual IVs (and include IV sets if any)
else:
    Indep_Var_List = Indep_Var_List[1:len(Indep_Var_List)] # otherwise, remove the "empty" space that's produced by the missing individual IVs


    # ~~ Create independent variable combination list ~~ #
    
Combination_List = list(map(it.combinations, # use map() function to apply combinations function to ...
                            list(it.repeat(Indep_Var_List, len(Indep_Var_List))), #... the IV list - which is repeated so that...
                            range(1, len(Indep_Var_List)+1))) # ... each number of combinations of a specific number of elements can be applied to get all possible combinations (note: saved as a combination object to be evaluated later and not as the list of combinations)

Total_Indep_Vars = len(Combination_List) # number of IVs in model

Total_Models_to_Estimate = 2**Total_Indep_Vars - 1 # total number of models to estimate

sfi.SFIToolkit.stata("display _newline \"{txt}Total of {res}" + str(2**len(Indep_Var_List)-1) + " {txt}regressions\"")

if Total_Indep_Vars > 4:
    sfi.SFIToolkit.stata("display _newline \"{txt}Progress in running all regression subsets\" _newline " + 
        "\"{res}0%{txt}{hline 6}{res}50%{txt}{hline 6}{res}100%\"")
    
    print(".", end="")


    # ~~ Define function to call regression model in Stata ~~ #
     
def st_model_call(Indep_Var_combination, report): # standard function to call and catch stata model
    sfi.SFIToolkit.stata( "quietly " +
                         Stata_Regression + " " + 
                         Dep_Var + " " +   
                         " ".join(Indep_Var_combination) + " " +
                         AllSubsets_Indep_Vars +
                         " if " + If_Conditions +
                          "," + Regress_Options ) 
    
    if report: print(".", end="") 
    
    return( (Indep_Var_combination,
             sfi.Scalar.getValue(Fit_Statistic)) ) 

def st_model_call_mi(Indep_Var_combination, report):  # multiple imputation-based function to call and catch stata model
    sfi.SFIToolkit.stata( "quietly mi estimate, saving(" +
                         Mult_Impute_File + ", replace) " +
                         Mult_Impute_Opts + ":" +
                         Stata_Regression + " " + 
                         Dep_Var + " " +   
                         " ".join(Indep_Var_combination) + " " +
                         AllSubsets_Indep_Vars +
                         " if " + If_Conditions +
                          "," + Regress_Options ) 
    
    Imputations = [int(imputation) for imputation in Mult_Imputes_toUse.split(" ")] # note imputations used 
    
    Mult_Impute_FitStat = 0
    
    for imputation in Imputations: 
        sfi.SFIToolkit.stata("estimates use " + Mult_Impute_File + ", number(" + str(imputation) + ")")
        
        Mult_Impute_FitStat = Mult_Impute_FitStat + sfi.Scalar.getValue(Fit_Statistic) # add to running sum of fit statistics
    
    if report: print(".", end="")
    
    Mult_Impute_FitStat = Mult_Impute_FitStat/len(Imputations) #Imputations # average the fit statistic values
    
    return( (Indep_Var_combination,
             Mult_Impute_FitStat) ) 


    # ~~ Obtain all subsets regression results ~~ #
    
Ensemble_of_Models = [] # initialize ensemble list container

"""
'Ensemble_of_Models' is structured such that:
1. Top level is results by number of IVs in the model
2. Middle level is model within a number of IVs
3. Bottom level is a specific result from 'st_model_call'
"""

ensemble_begin = 0 # note where the ensemble of models has begun for this set of IVs (relevant for tracking progress only)

ensemble_end = 0 # note where the ensemble of models has ended for this set of IVs (relevant for tracking progress only)

if Total_Indep_Vars > 4: # if at least 20 models, note which models will report a '.' when estimating
    flag_list = [int(twentieth/20*Total_Models_to_Estimate) for twentieth in range(1,21)]
    
else: flag_list = [] # if not at least 20 models, do not track  progress
    
for number_of_Indep_Vars in range(Total_Indep_Vars): # applying the modeling function across all IV combinations at a distinct number of IVs
    
    ensemble_end = (ensemble_end +
                    fctl(Total_Indep_Vars)/(fctl(number_of_Indep_Vars+1)*fctl(Total_Indep_Vars-(number_of_Indep_Vars+1)))) # ending point of ensemble -- cumulative
    
    report_model_list = [x in flag_list for x in range(int(ensemble_begin), int(ensemble_end))] # note which models will produce '.' when estimated in progress tracker
    
    if Mult_Impute_Flag: # if multiple imputation...
        Models_at_Indep_Var_number = list( map(st_model_call_mi,
                                   list(Combination_List[number_of_Indep_Vars]), report_model_list ) )
    
    else: # ... otherwise normal function
        Models_at_Indep_Var_number = list( map(st_model_call, 
                                   list(Combination_List[number_of_Indep_Vars]), report_model_list ) )
    
    Ensemble_of_Models.append(Models_at_Indep_Var_number) 
    
    ensemble_begin = ensemble_end # update where the ensemble tracker will begin for next round


    # ~~ Process all subsets - find the increments  ~~ #
    
Model_List = [[list(model) for model in Ensemble_of_Models[0]]]  # evaluate the map-ped models and record them - start with the single IV models...

for model in range(len(Model_List[0])): # ...for the single IV models...
    Model_List[0][model][1] = Model_List[0][model][1]-FitStat_Adjustment #... have to remove constant model results as well as all subets results

#print("model list:", Model_List) #//

for number_of_Indep_Vars in range(1, len(Ensemble_of_Models)): # when >1 IV in the model, processing needed...
    
    Model_Incremented = []  # initialize/reset container for finding subset
    
    Indep_Var_Set_at_1lessIndep_Var = [set(Candidate_Indep_Var_Set[0]) for Candidate_Indep_Var_Set in Ensemble_of_Models[number_of_Indep_Vars-1]] # collect all sets IVs (coerced to be a set object), specifically all sets at one less IV in the model than the current number of IVs
    
    for model in range(0, len(Ensemble_of_Models[number_of_Indep_Vars])): # loop through all models at a specific number of IVs in the model...

        Indep_Var_Set = set(Ensemble_of_Models[number_of_Indep_Vars][model][0]) # IV set for a focal model; coerced to be set object
        
        for at1less_model in range(0, len(Indep_Var_Set_at_1lessIndep_Var)): # loop through all models at one less than the specific number of IVs in the model...

            if Indep_Var_Set_at_1lessIndep_Var[at1less_model].issubset(Indep_Var_Set): # if IV set at one less is a subset of the predictors in the focal model...
                
                Model_Incremented.append( 
                [ Ensemble_of_Models[number_of_Indep_Vars][model][0], # append IV names at focal ...
                  Ensemble_of_Models[number_of_Indep_Vars-1][at1less_model][0], # ...IV names at one less...
                  Ensemble_of_Models[number_of_Indep_Vars][model][1] - Ensemble_of_Models[number_of_Indep_Vars-1][at1less_model][1] ] # ...and the increment to the fit metric
                )
                
    Model_List.append(Model_Incremented) 
        

"""
'Model_List' is structured such that:
1. Top level is results by number of IVs in the model
2. Middle level is model within a number of IVs
3. Bottom level is a specific increment's information (full_model, reduced_model, fit metric difference)
"""


    # ~~ Obtain complete and conditional dominance statistics  ~~ #
    
Conditional_Dominance = [] # conditional dominance container

if Complete_Flag: Complete_Dominance = [] # complete dominance container

for Indep_Var in range(0, len(Model_List[0])): # for each IV in the model...
    
    Conditional_atIndep_Var = [] # initialize/reset container for conditional dominance
    
    Conditional_atIndep_Var.append(Model_List[0][Indep_Var][1]) # for IV alone - copy fit statistic

    Indep_Varname = set(Model_List[0][Indep_Var][0]) # record name of focal IV; coerce to set

    if Complete_Flag: Complete_atIndep_Var = [
        [Other_Indep_Var[1] < Model_List[0][Indep_Var][1]] # compare fit statistic values (is focal IV larger than other IV?) ...
        
        for Other_Indep_Var in Model_List[0][0:len(Model_List[0])] ] #... for other IVs alone (will compare to self also)

    for number_of_Indep_Vars in range(1, len(Model_List)): # for all numbers of IVs greater than 1...
 
        Relevant_Increments = [] # initialize/reset container for collecting specific/relevant conditional dominance increments
        
        for model in range(0, len(Model_List[number_of_Indep_Vars])): # for each individual model within a specific number of IVs...

            proceed_to_record = ( Indep_Varname.issubset( set(Model_List[number_of_Indep_Vars][model][0]) ) and not # flag this entry for recording if the focal IV name is in the IV set...
                         Indep_Varname.issubset( set(Model_List[number_of_Indep_Vars][model][1]) ) ) # ...but is _not_ in the IV set less one - thus, the fit statistic here is a valid "increment" for the focal IV

            if proceed_to_record: 
                Relevant_Increments.append(Model_List[number_of_Indep_Vars][model][2]) # always collect the fit statistic for conditional dominance computations
                
                if Complete_Flag:
                    for other_model in range(0, len(Model_List[number_of_Indep_Vars])): # also proceed to collect complete dominance data using this loop comparing to all other models within this number of IVs to find relevant comparisons

                       relevant_complete = ( # a relevant complete dominance comparsion is found when ...
                            set(Model_List[number_of_Indep_Vars][model][0]).issubset( set(Model_List[number_of_Indep_Vars][other_model][0]) ) and # ...the focal full model and the full other model have the same IV set (the only way they can be a 'subset' here) ...
                            len(set(Model_List[number_of_Indep_Vars][model][1]).difference( set(Model_List[number_of_Indep_Vars][other_model][1])) ) == 1 ) #... but their reduced IV set differs by one IV (this ensures it is not trying to compare the subset to itself)
 
                       if relevant_complete: 
                            MatrixLocation_Complete = [Position_IV[0] for Position_IV in Model_List[0]].index(( # when a relevant comparison, obtain the index value for ...
                                       set(Model_List[number_of_Indep_Vars][model][1]).difference( set(Model_List[number_of_Indep_Vars][other_model][1]) ).pop(), )) #... the different element in the reduced model (to place it in the correct "row" for the dominance matrix/list)
                            
                            Complete_atIndep_Var[MatrixLocation_Complete].append( #at the correct location in the complete dominance matrix, append...
                                Model_List[number_of_Indep_Vars][other_model][2] < Model_List[number_of_Indep_Vars][model][2] ) # ...whether the other model's increment is bigger than the focal
                            
                        
        Conditional_atIndep_Var.append(stat.mean(Relevant_Increments)) # compute conditional dominance at number of IVs for specific IV and append
    
    Conditional_Dominance.append(Conditional_atIndep_Var) # append full row of IV's conditional dominance statistics
    
    if Complete_Flag: Complete_Dominance.append(Complete_atIndep_Var) # append full row of IV's complete dominance logicals/designations

if Complete_Flag:
    Complete_Dominance = [list(map(all, Indep_Var)) for Indep_Var in Complete_Dominance] # for each focal IV, make list comprehension that flags whether at each comparison (i.e., other IV) are all entries (i.e., specific comarisons between similar models) in list True?

    Complete_Dominance = [[int(IV_Other_Compare) for IV_Other_Compare in Indep_Var] for Indep_Var in Complete_Dominance] # for each IV and other comparison, change boolean to integer for use in Stata


if Conditional_Flag:
    sfi.Matrix.create('r(cdldom)', len(Ensemble_of_Models), len(Ensemble_of_Models), 0) # create conditional dominance matrix container in Stata
    sfi.Matrix.store('r(cdldom)', Conditional_Dominance) # post conditional dominance matrix

if Complete_Flag:
    sfi.Matrix.create('r(cptdom)', len(Ensemble_of_Models), len(Ensemble_of_Models), 0) # create complete dominance matrix container in Stata
    sfi.Matrix.store('r(cptdom)', Complete_Dominance) # post complete dominance matrix


    ## ~~ Compute general dominance and fit statistic  ~~ ##
    
General_Dominance = list(map(stat.mean, Conditional_Dominance)) # average conditional dominance statistics to produce general dominance

FitStat = stat.fsum(General_Dominance) + FitStat_Adjustment # adjust overall fit statistic by replacing all subsets component and constant model component
sfi.Scalar.setValue('r(fs)', FitStat) # post overall fitstat

sfi.Matrix.create('r(domwgts)', 1, len(Ensemble_of_Models), 0) # create general dominance matrix container in Stata
sfi.Matrix.store('r(domwgts)', General_Dominance) # post general dominance statistics

sfi.Matrix.create('r(sdomwgts)', 1, len(Ensemble_of_Models), 0) # create standardized general dominance matrix container in Stata
sfi.Matrix.store('r(sdomwgts)', list(map(lambda Gen_Dom_Stat: Gen_Dom_Stat/FitStat, General_Dominance)) ) # produce standardized general dominance

General_Dominance_Ranks = [sorted(General_Dominance, reverse = True).index(Indep_Var)+1 for Indep_Var in General_Dominance] # rank general dominance statistics
sfi.Matrix.create('r(ranks)', 1, len(Ensemble_of_Models), 0) # create general dominance ranking matrix container in Stata
sfi.Matrix.store('r(ranks)', General_Dominance_Ranks) # post general dominance statistic rankings
