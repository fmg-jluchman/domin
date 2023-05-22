*! dominance.mata version 0.1.0  xx/xx/202x Joseph N. Luchman

version 12.1

**# domin_specs structure
	
	//each of these elements are called by the regression model

mata:

mata set matastrict on

struct domin_specs {
	
	string scalar mi, reg, dv, all, weight, exp, touse, regopts 
	
}

end

**# dominance() function

	//Mata function to compute all combinations of predictors or predictor sets run all subsets regression, and compute all dominance criteria
	
mata: 

mata set matastrict on

void dominance(struct domin_specs scalar model_specs, pointer scalar model_call,  
	string scalar cdlcompu, string scalar cptcompu, string scalar iv_string, 
	real scalar all_subsets_fitstat, real scalar constant_model_fitstat)
{
	/*# object declarations*/
	real matrix IV_indicator_matrix, wset_indicator_matrix,
		basic_wsets, expd_wsets, expd_combn, combin_at_order, all_of_wset,
		fitstat_increments, permutes_by_increment,
		weighted_log_increments, 
		conditional_dominance, cdl_increments,
		complete_dominance

	real rowvector fitstat_vector, orders, 
		IVs_in_row, IV_rest_combs,
		which_focal_IV, which_other_IVs,
		find_focal_pattern, find_other_pattern,
		general_dominance, which_basic, which_wset, 
		row_rownumbers, col_rownumbers, row_IVsort, col_IVsort,
		row_IVsort_numbers, col_IVsort_numbers

	real colvector wset_lengths, wset_reps, 
		rows_have_IV, IVs_at_order, IVs_to_end, not_wset,
		other_column_indexes, cdl_model_sizes, cdl_sums,
		cdl_wset_sizes, cdl_wset_denom, 
		cpt_setup_vec, cpt_iterator, cpt_IVs, cpt_wsets, sort_omit_xy
		
	string colvector IVs, wset_IVs, wset_IVs_expd
	
	real scalar x, y, number_of_IVs, posit, display,
		expect_size, log_permutes, fitstat

	string scalar IVs_in_model
	
	transmorphic t, cdl_wset_panel, cpt_permute_container
	
	/* ~~ step 1: parse the IV inputs ~~ */
	
	// This section sets up parsing rules to separate 'iv_string' into tokens that differentiates between ...
		// ... basic sets in '<>' characters from those in wsets in '{}' characters.
		// The ending 'IVs' matrix is a string row vector of all the IVs, basic sets, and wsets.
		
	t = tokeninit((" "), (" "),("<>", "{}")) 
	
	tokenset(t, iv_string) 
	
	IVs = tokengetall(t)' 
	
	"IVs as parsed" // ~~
	IVs // ~~
	
	/* ~~ step 2: parse and process wsets and basic sets ~~ */
	
	// This section first identifies whether there are wsets present in the 'IVs' vector.
		// If there are wsets, they are separated into a different vector: 'wset_IVs' and their enclosing characters '{}' are removed.
		// The remaining IVs and basic sets remain in the 'IVs' vector. 
		// If there are remaining IVs, each element is checked to see whether it is a basic set and their enclosing characters '<>' are removed.
	
	/* ~~~ substep 2a: process wsets ~~~ */
	if ( any(substr(IVs, 1, 1):=="{") ) { 
		
		wset_IVs = select(IVs, (substr(IVs, 1, 1):=="{") ) 
		
			for (x = 1; x <= rows(wset_IVs); x++) { 
				
				wset_IVs[x] = substr(wset_IVs[x], 1, strlen(wset_IVs[x]) - 1 ) 
					
				wset_IVs[x] = substr(wset_IVs[x], 2, strlen(wset_IVs[x]) ) 
				
			}
		
		IVs = select(IVs, (substr(IVs, 1, 1):!="{") ) 
		
	}
	
	/* ~~~ substep 2b: process IVs and basic sets (if any) ~~~ */
	if (length(IVs)) { 
	
		if ( any(substr(IVs, 1, 1):=="<") ) { 
		
			for (x = 1; x <= rows(IVs); x++) { 
			
				if (substr(IVs[x], 1, 1) == "<") { 
				
					IVs[x] = substr(IVs[x], 1, strlen(IVs[x]) - 1 ) 
					
					IVs[x] = substr(IVs[x], 2, strlen(IVs[x]) )
					
				}
				
			}
		
		}
	
	}
	
	/* ~~ step 3: compute all combinations of IVs and basic sets ~~*/
	
	// This section computes all combinations of IVs and basic sets using the -domin_combn()- function.
		// If there are no elements in 'IVs', the 'IV_indicator_matrix' matrix is empty.
	
	number_of_IVs = rows(IVs) 
	
	/*??*/ //if (number_of_IVs > 12) printf("\n{txt}Computing all independent variable combination sub-models\n")

	IV_indicator_matrix = domin_combn(number_of_IVs)
	
	/* ~~ step 4: set-up wsets for combination construction ~~ */
	
	// This section collects the number of elements in each wset, if any, in the 'wset_lengths' row vector. 
		// This section also maps wsets to IVs using the 'wset_reps' row vector.
		// When there are no wsets, the both vectors record that the wsets are size 1 and each IV is its own wset.
	
	if ( rows(wset_IVs) ) { 
		
		/* ~~~ substep 4a: if wsets, determine number of IV elements in each wset ~~~ */
		wset_lengths = J(rows(wset_IVs), 1, .)
		
		for (x = 1; x <= rows(wset_IVs); x++) { 
			
			wset_lengths[x] = length( tokens(wset_IVs[x]) ) 
			
		}
		
		/* ~~~ substep 4b: if wsets, map wset membership to IVs ~~~ */
		wset_reps = J(sum(wset_lengths), 1, .) 
			
		posit = 1 
		
		for (x = 1; x <= length(wset_lengths); x++) { 
			
			for (y = 1; y <= wset_lengths[x]; y++) { 
				
				wset_reps[posit] = x 
				
				posit++ 
				
			}
			
		}
		
		/* ~~ step 4+: construct all combinations of wsets ~~ */
	
		// This section is only relevant when there are wsets. 
			// If there are no wsets, the 'IV_indicator_matrix' as constructed in step 3 is all that is required.
			// If there are wsets, a container matrix is created, 'wset_indicator_matrix', that is filled in with all relevant combinations.
			// The first set of combinations, 'basic_wsets', computes all combinations as though they were basic sets --  all IV elements in or out.
			// Following including 'basic_wsets', the all combinations of wsets are computed.
			// For each wset, the 'expect_size' scalar is computed indicating the number of combinations of IVs allowing IVs to vary within set but holding constant across sets is computed.
			// The combinations are computed (using -domin_combn-) and spread across the levels of the other wsets (again that are constant within set here).
			// The combinations of all wsets in 'wset_indicator_matrix' are then spread across all IV and basic set combinations, updating 'IV_indicator_matrix'.
			// The 'IVs', 'wset_lengths', 'wset_reps', and 'number_of_IVs' metadata are also updated.
		
		wset_indicator_matrix = 
			J( sum(wset_lengths), 
				sum( ((2:^wset_lengths):-2):*(2^(rows(wset_lengths)-1)) ) + 
					2^rows(wset_lengths), 
				. ) 
		
		/* ~~~ substep 4+a: wsets-as-basic sets ~~~ */
		basic_wsets = domin_combn( rows(wset_lengths) ) 
		
		basic_wsets = basic_wsets[wset_reps, ] 
		
		wset_indicator_matrix[, 1..cols(basic_wsets)] = basic_wsets 
		
		posit = cols(basic_wsets) 
		
		/* ~~~ substep 4+b: wsets across one another ~~~ */
		for (x = 1; x <= rows(wset_IVs); x++) { 
		
			expect_size = 
				(2^wset_lengths[x] - 2)*( 2^(length(wset_lengths) - 1) ) 
		
			expd_combn = domin_combn( length( tokens(wset_IVs[x]) ) ) 
			
			expd_combn = expd_combn[,2..cols(expd_combn)-1] 
			
			if (cols(expd_combn) != expect_size) { 
				
				expd_combn = J(1, expect_size/cols(expd_combn), expd_combn) 
				
			}
			
			expd_wsets = 
				select(basic_wsets, 
					( colsum( select(basic_wsets, wset_reps:==x) ):==0 ) ) 
			
			
			expd_wsets = 
				( expd_wsets # J(1, expect_size/cols(expd_wsets), 1) ) 
			
			expd_wsets[ select( 1::length(wset_reps), wset_reps:==x ) , ] = 
				expd_combn 
			
			wset_indicator_matrix[,posit+1..posit+cols(expd_wsets)] = 
				expd_wsets 
			
			posit = posit + cols(expd_wsets) 
		
		}
		
		/* ~~~ substep 4+c: wset combinations across IVs and basic sets ~~~ */
		IV_indicator_matrix = 
			IV_indicator_matrix # J(1, cols(wset_indicator_matrix), 1) 
		
		wset_indicator_matrix = 
			J(1, cols(IV_indicator_matrix)/cols(wset_indicator_matrix), 
				wset_indicator_matrix)
				
		IV_indicator_matrix = 
			IV_indicator_matrix \ wset_indicator_matrix
		
		/* ~~~ substep 4+d: update IV names with IVs in each wset ~~~ */
		wset_IVs_expd = J( sum(wset_lengths), 1, "") 
		
		posit = 1		
		
		for (x = 1; x <= length(wset_IVs); x++)	{
			
			wset_IVs_expd[posit..(posit+wset_lengths[x]-1)] = 
				tokens(wset_IVs[x])'
			
			posit = posit + wset_lengths[x]
			
		}
		
		/* ~~~ substep 4+e: include IVs and basic sets in wset metadata vectors ~~~ */
		if (length(IVs)) {
			
			wset_reps = 
				( 1::length(IVs) ) \ 
					wset_reps:+length(IVs) 
					
			wset_lengths = 
				J(length(IVs), 1, 1) \ 
					wset_lengths
			
		}
		
		/* ~~~ substep 4+f: include IVs and basic sets in wset metadata vectors ~~~ */
		IVs = IVs \ wset_IVs_expd
		
		number_of_IVs = number_of_IVs + rows(wset_IVs_expd)
		
	}
	
	/* ~~~ substep 4c: if no wsets, create wset metadata vectors with IVs and basic sets ~~~ */
	else {
		
		wset_reps = 1::length(IVs)
		
		wset_lengths = J(length(IVs), 1, 1)
		
	}
	
	"wset_reps" // ~~
	wset_reps // ~~
	"wset_lengths" // ~~
	wset_lengths // ~~
	
	/* ~~ step 5: orders and sort ~~ */
	
	// This section determines the number of IVs in each model/row in the 'orders' column vector.
		// The 'IV_indicator_matrix' with all models to be estimated from the data are bound to the 'orders' vector, transposed, and sorted such that fewer IVs in the model are at lower row index numbers.
	
	orders = colsum(IV_indicator_matrix)
	
	IV_indicator_matrix = sort((orders \ IV_indicator_matrix)', 1..rows(IV_indicator_matrix)+1) // ~~
	
	"final, sorted indicator matrix" // ~~
	IV_indicator_matrix // ~~
	
	/* ~~ step 6: estimate all fit statistics ~~ */
	
	// This section loops over all combinations of IVs defined in 'IV_indicator_matrix' and submits them, as names, to the focal model.
		// A key activity below is in displaying progress toward completing all model runs when there are 6 or more columns (i.e., 5 IVs) in 'IV_indicator_matrix'.
		// Below a container row vector, 'fitstat_vector', is created to hold all fit statistic values and the already known value with no IVs other than those in the -all- submodels value is placed into index value 1.
		// The 'IVs_in_model' string scalar is created using Mata's string selection-via-multiplication method and submitted as an argument to '*model_call' a pointer to -domin_call- defined below.
		// The value returned by -domin_call- is placed directly into each element of 'fitstat_vector'.
	
	/* ~~~ substep 5a: establish display routine ~~~ */
	display = 1 
	
	printf("\n{txt}Total of {res}%f {txt}sub-models\n", rows(IV_indicator_matrix) )
	
	if (rows(IV_indicator_matrix) > 6) {
	
		printf("\n{txt}Progress in running all sub-models\n{res}0%%{txt}{hline 6}{res}50%%{txt}{hline 6}{res}100%%\n")
		
		printf(".")
		
		displayflush()
		
	}
	
	/* ~~~ substep 5b:  ~~~ */
	fitstat_vector = J(1, rows(IV_indicator_matrix), .) //pre-allocate container vector that will contain fitstats across all models
	
	fitstat_vector[1] = constant_model_fitstat + all_subsets_fitstat 
	
	for (x = 2; x <= rows(IV_indicator_matrix); x++) { 
	
		if (rows(IV_indicator_matrix) > 4) {
	
			if (floor(x/rows(IV_indicator_matrix)*20) > display) {
			
				printf(".")
				
				displayflush()
				
				display++	
				
			}
			
		}
		
		IVs_in_model = 
			invtokens( IV_indicator_matrix[x, 2..cols(IV_indicator_matrix)]:*(IVs') ) 
		
		fitstat_vector[x] = (*model_call)(IVs_in_model, 
			all_subsets_fitstat, constant_model_fitstat, model_specs) 

	}
	
	"fitstat vector" // ~~
	fitstat_vector // ~~
	
	/* ~~ step 6: obtain all fit statistic increments ~~ */
	
	// This section computes the increment the focal column's IV has from all other column IVs included in the same row/model of 'IV_indicator_matrix'.
		// All IVs and basic sets do not have valid increments that include partial sets of wset elements. The 'all_subsets_fitstat' matrix flags which rows have either all or none of the elements of each wset to adjust 'IV_indicator_matrix'.
		// The 'all_of_wset' matrix is used to remove indications/'1's from 'IV_indicator_matrix' when the column IV is an IV or basic set during increment computation with 'fitstat_increments'.
		// Increment computation first obtains the location of all indications/'1's from 'IV_indicator_matrix' for each column IV. 
		// For each row in which the focal column IV is present in the model, the pattern of other column IVs is obtained and compated to all other rows in which the focal column IV is not included.
		// The row indexes for the row with the focal column IV given a specific pattern of the other column IVs is matched with the row indexes for the row without the focal column IV given the same pattern of the other column IVs.
		// The values in the matched row indexes are subtracted and recorded in 'fitstat_increments' in the same row as the original column IV's location and in the focal IV's column.
	
	/* ~~~ substep 6a: identify which wsets contain all or none of within set IV elements ~~~ */
	all_of_wset = 
		J(rows(IV_indicator_matrix), length(wset_lengths), 1) 
		
	for (x = 1; x <= length(wset_lengths); x++) {
		
		if (wset_lengths[x] == 1) continue 
		
		all_of_wset[,x] = 
			( rowsum( 
				IV_indicator_matrix[, 
					select(1..cols(IV_indicator_matrix), 
						(0, wset_reps'):==x ) ] ):==wset_lengths[x] ):+( 
				rowsum( 
					IV_indicator_matrix[, 
						select(1..cols(IV_indicator_matrix), 
							(0 , wset_reps'):==x ) ] ):==0 )
		
	}
	
	/* ~~~ substep 6b: compute fit statistic increments ~~~ */
	fitstat_increments =
		J(rows(IV_indicator_matrix), cols(IV_indicator_matrix)-1, 0) 
	
	for (x = 2; x <= cols(IV_indicator_matrix); x++) {
		
		if (wset_lengths[ wset_reps[x-1] ] == 1) { 
			
			which_focal_IV = 
				selectindex( IV_indicator_matrix[,x]:*(
					rowsum(all_of_wset):==cols(all_of_wset) ) ) 
			
			which_other_IVs = 
				selectindex( !IV_indicator_matrix[,x]:*(
					rowsum(all_of_wset):==cols(all_of_wset) ) ) 
			
		}
		
		else { 
			
			which_focal_IV = selectindex( IV_indicator_matrix[,x] ) 
			
			which_other_IVs = selectindex( !IV_indicator_matrix[,x] ) 
			
		}
		
		other_column_indexes = 
			select( (2..cols(IV_indicator_matrix)), 
				(2..cols(IV_indicator_matrix)):!=x ) 
		
		find_focal_pattern = 
			select(
				J( length(which_other_IVs), 
					1, 
					which_focal_IV), ///
				J( length(which_other_IVs), 
					1, 
					rowsum( 
						IV_indicator_matrix[which_focal_IV, other_column_indexes]:*10:^(cols(IV_indicator_matrix)-3..0) ) ):==( 
							rowsum( IV_indicator_matrix[which_other_IVs, other_column_indexes]:*10:^(cols(IV_indicator_matrix)-3..0) ) # 
								J( rows(which_focal_IV), 1, 1) ) )
					
		find_other_pattern = 
			select(
				J( length(which_focal_IV), 
					1, 
					which_other_IVs), ///
				J( length(which_focal_IV), 
					1, 
					rowsum( 
						IV_indicator_matrix[which_other_IVs, other_column_indexes]:*10:^(cols(IV_indicator_matrix)-3..0) ) ):==( 
							rowsum( IV_indicator_matrix[which_focal_IV, other_column_indexes]:*10:^(cols(IV_indicator_matrix)-3..0) ) # 
								J( rows(which_other_IVs), 1, 1) ) )
				
		fitstat_increments[find_focal_pattern, x-1] = 
			( fitstat_vector[find_focal_pattern]:-fitstat_vector[find_other_pattern] )'
		
	}
	
	"increments"
	fitstat_increments // ~~
	
	/* ~~ step 7: increment weighting ~~ */
	
	// This section weights the increments by the number of permutations in which that increment appears.
		// 'log_permutes' establishes the total number of permutations given the number of IVs, basic sets, and wsets. This value is log transformed for better storage performance and convenience in computation.
		// 'permutes_by_increment' collects the number of (log) permutations associated with each increment.
		// Note the differences between IVs and basic sets (i.e., 'wset_lengths' of 1) as compated to wsets in terms of how permutations are registered.
		// Note also that the increment log permutations are registered in the same locations as their increments.
	
	
	log_permutes = 
		lnfactorial( length(wset_lengths) ) + sum( lnfactorial(wset_lengths) )
	
	permutes_by_increment = 
		J(rows(fitstat_increments), cols(fitstat_increments), .)
	
	for (x = 1; x <= cols(permutes_by_increment); x++) { 
		
		rows_have_IV = selectindex( fitstat_increments[,x]:>0 ) 
			
		IV_rest_combs = selectindex( (1..length(wset_reps)):!=x )
		
		not_wset = selectindex( (1..length(wset_lengths)):!=wset_reps[x]  )
		
		for (y = 1; y <= length(rows_have_IV); y++) {
			
			if ( wset_lengths[ wset_reps[x] ] == 1 ) 
				IVs_in_row = fitstat_increments[ rows_have_IV[y], IV_rest_combs ]:>0
			
			else 
				IVs_in_row = IV_indicator_matrix[ rows_have_IV[y], IV_rest_combs:+1 ]
			
			IVs_at_order = 
				panelsum( IVs_in_row', 
					panelsetup( wset_reps[ IV_rest_combs ] , 1) ) 
			
			IVs_to_end = 
				panelsum( J(length(IVs_in_row), 1, 1), 
					panelsetup( wset_reps[ IV_rest_combs ] , 1) ):-IVs_at_order
					
			if ( wset_lengths[ wset_reps[x] ] == 1 ) 
				permutes_by_increment[rows_have_IV[y], x] = 
					lnfactorial( length( select(IVs_at_order, IVs_at_order:>0) ) )  + 
					sum( lnfactorial(IVs_at_order) ) + 
					lnfactorial( length( select(IVs_to_end, IVs_to_end:>0) ) ) + 
					sum( lnfactorial(IVs_to_end) )
			
			else 
				permutes_by_increment[rows_have_IV[y], x] = 
					lnfactorial( sum(IVs_at_order[ not_wset ]:>0) ) + 
					sum( lnfactorial( IVs_at_order ) ) + 
					lnfactorial( sum(IVs_to_end[ not_wset ]:>0) ) + 
					sum( lnfactorial( IVs_to_end ) ) 
			
		}
			
	}
	
	// !! This will result in errors if the increments are negative !! <- rescale all values such that the min value is '+' then de-scale back??
	weighted_log_increments = ln(fitstat_increments) + permutes_by_increment
	
	"permutes_by_increment" // ~~
	exp( permutes_by_increment ) // ~~
	"permutes_by_increment - summed" // ~~
	colsum( exp( permutes_by_increment ) ) // ~~
	"number of orderings" // ~~
	exp( log_permutes ) // ~~
	"combined" // ~~
	exp( weighted_log_increments ) // ~~
	
	/* ~~ step 8: general dominance ~~ */
	
	// This section computes general dominance statistics (i.e., Shapley or Owen values)
		// The computations are implemented by summing 'weighted_log_increments' minus the 'log_permutes'. This value is also fed back to the Stata environment.
		// The ranks, summed overall fit statistic, and standardized/normed general dominance statistics are computed and fed back to the Stata environment.
	
	general_dominance = 
		colsum( 
			exp( weighted_log_increments:-J( 
				1, cols(weighted_log_increments), log_permutes ) ) )
				
	fitstat = rowsum(general_dominance) + all_subsets_fitstat + constant_model_fitstat
	
	"general dominance" // ~~
	general_dominance // ~~
	"fitstat" // ~~
	fitstat // ~~
	
	st_matrix("r(domwgts)", general_dominance) 

	st_matrix("r(sdomwgts)", general_dominance:/fitstat) 
	
	st_matrix("r(ranks)", mm_ranks(general_dominance'*-1, 1, 1)') 

	st_numscalar("r(fs)", fitstat) 
	
	/* ~~ step 9: standard conditional dominance ~~ */
	
	// This section prepares the results the computation of conditional domiance statistics - a version of Shapley/Owen values that works by evaluating increments by number of IVs in the model.
		// The section begins by counting the number of IVs, basic sets, and wsets in the model. 
		// Importantly, wsets are treated 'as though' they are basic sets in terms of aggregating results; the 'cdl_model_sizes' vector reflects this updated model size value.
		// All IVs and basic sets are identified using 'which_basic' that first identifies which elements have values of '1' in 'wset_lengths' which is then mapped to 'wset_reps'.
		// If there are IVs or basic sets, 'cdl_increments' is generated reflecting only those columns of 'fitstat_increments' associated with the IVs and basic sets which are subsequently filtered to those rows with non-0 entries.
		// Finally, the results are panelsum()-med by model size and are divided by the number of individual model combinations at each updated model size to get the average.
	
	if (strlen(cdlcompu) == 0) {
		
		conditional_dominance = J(length(wset_reps), length(wset_lengths), .)

		cdl_model_sizes = 	
			colsum(
				floor( 
					panelsum(
						IV_indicator_matrix[,2..cols(IV_indicator_matrix)]',
						panelsetup( wset_reps, 1 ) ):/wset_lengths ) )'
		
		which_basic = selectindex(wset_lengths:==1)'
		
		if ( length(which_basic) ) {
		
			cdl_increments =  
				(cdl_model_sizes, fitstat_increments[,which_basic] )
						
			cdl_increments = 
				select(cdl_increments, rowsum(cdl_increments[, 2..cols(cdl_increments)]):>0)
				
			
			cdl_sums = 
				(panelsum(cdl_increments[,2..cols(cdl_increments)], 
					panelsetup(cdl_increments, 1) ):/(
						panelsum(cdl_increments[,2..cols(cdl_increments)]:>0, 
							panelsetup(cdl_increments, 1) ) ))
					
			conditional_dominance[selectindex(which_basic),] = cdl_sums'
		
		}
		
		/* ~~ step 9+: wset-based conditional dominance ~~ */
	
		// This section is only relevant when there are wsets. 
			// If there are no wsets, the 'conditional_dominance' matrix will be filled in the last step and the loop below will terminate immediately as the size of 'which_wset' will be 0.
			// If there are wsets, a loop is initiated which considers each wset separately for conditional dominance statistic computation.
			// Within the loop, the first update to 'which_wset' identifies the colums of 'fitstat_increments' that are associated with the focal wset.
			// The wsets use a different model size computation where the values from the overall update to model size is incremented by 1--except for those models in which all of the wset IVs are included.
			// The model sizes and fitstat increments are then selected and sorted in preparation for panelsum()-ming.
			// Prior to using the panelsum(), the data are panelsetup()-ed and the individual panels (i.e., model sizes) are looped through in a second loop.
			// The second loop effectively initiates a full Shapley value decompositon for each model size. 
				// Specifically, the submatrix associated with just one model size is extracted, from which the values are sorted in preparation for a panelsum().
				// The panelsum() is used only to get the number of combinations observed including the focal IV at the within wset by within-set model size to use as a denominator.
				// 'cdl_wset_denom' is then 'spread' using repeated index values to all the observed increments and column summed to obtain the final wset conditional dominance statistic at the focal model size.
			// 'which_wset' is refreshed prior to re-initiating the loop for the next wset (if any).
			
		which_wset = selectindex(wset_lengths:>1)'
		
		for (x = 1; x <= length(which_wset); x++) {
				
			which_wset = 
				selectindex( (wset_reps'):==which_wset[x] )
			
			cdl_wset_sizes = 
				(cdl_model_sizes:+1):-(
					rowsum(
						IV_indicator_matrix[,which_wset:+1]):==length(
							which_wset) )
					
			cdl_increments =  
				(cdl_wset_sizes, fitstat_increments[,which_wset] )
					
			cdl_increments = 
				select(cdl_increments, rowsum(cdl_increments[, 2..cols(cdl_increments)]):>0)
						
			cdl_increments = 
				sort(cdl_increments, 1..cols(cdl_increments))
			
			cdl_wset_panel = 
				panelsetup(cdl_increments, 1)
				
			for (y = 1; y <= max(cdl_increments[,1]); y++) {
				
				cdl_sums = 
					panelsubmatrix(
						cdl_increments[,2..cols(cdl_increments)], y, 
							cdl_wset_panel)
				
				cdl_sums = 
					( sort( 
						(rowsum(cdl_sums:>0), cdl_sums), 
							1..cols(cdl_sums)+1 ) )[,2..cols(cdl_sums)+1]				
				
				cdl_wset_denom = 
					panelsum( cdl_sums[,1]:>0 , 
						panelsetup( rowsum(cdl_sums:>0), 1 ) 
							):*cols(cdl_sums)
				
				cdl_sums = 
					colsum( 
						cdl_sums:/cdl_wset_denom[rowsum(cdl_sums:>0)] )
				
				conditional_dominance[which_wset, y] = cdl_sums'
			
			}
			
			which_wset = selectindex(wset_lengths:>1)'
				
		}

		"conditional_dominance" // ~~
		conditional_dominance // ~~
		
		st_matrix("r(cdldom)", conditional_dominance) 
	
	}
	
	/* ~~ step 10: complete dominance ~~ */
	
	// This section computes complete dominance statistics.
		// These values are an experimental statistic representing the proportion of models that the focal IV produces a bigger value than the comparison IV.
		// Note that any comparison between an IV in a wset with an IV in another wset as well as any non-wset IV or basic set is currently considered to be invalid.
		
		// !! ~~ redo - ensure these values are, like the gen + cond doms, subdivisions of R2 by sets to allow for all IVs to be compared.
		
		// !! ~~ perspective taken below is that it is less expensive to get row indexes, filter, sort, filter than to filter the matrix itself and store it's values - especially as this scales
	
	if (strlen(cptcompu) == 0) {
		
		complete_dominance = 
			J(length(wset_reps), length(wset_reps), .) 
			
		cpt_setup_vec = 
			( J(2, 1, 1) \ J( length(wset_lengths)-2, 1, 0) )
	
		cpt_permute_container = cvpermutesetup(cpt_setup_vec) 
		
		// selects wsets (invididual ivs are a kind of wset with 1 value/ basic set is wset with 1 value and multiple ivs) to compare
		cpt_iterator = selectindex( cvpermute(cpt_permute_container) )
		
		while ( length(cpt_iterator) ) {
			
			cpt_iterator // ~~
			
			// creates matrix of 'wset_reps' repeated to find which IV-cols are associated with the wsets - accommodates wsets of size > 1
			cpt_IVs = 
				selectindex( 
					rowsum( 
						J(1, length(cpt_iterator), wset_reps
							):==cpt_iterator' ) )
				
			cpt_IVs // ~~
			
			// which rows are associated with observed increments of the first wset (the iv in the row of cpt)
			row_rownumbers = 
				selectindex( 
					rowsum(
						fitstat_increments[, 
							cpt_IVs[
								1..wset_lengths[ cpt_iterator[1] ] ] ]:>0):==wset_lengths[ 
									cpt_iterator[1] ] )
									
			row_rownumbers // ~~

			// which rows are associated with observed increments of the second wset (the iv in the col of cpt)
			col_rownumbers = 
				selectindex(
					rowsum(
						fitstat_increments[, 
							cpt_IVs[ 
								wset_lengths[ 
									cpt_iterator[1] ]+1..wset_lengths[ 
										cpt_iterator[1] ]+wset_lengths[ 
											cpt_iterator[2] ] ] ]:>0):==wset_lengths[ 
												cpt_iterator[2] ] )
			
			col_rownumbers // ~~
			
			for (x = cpt_IVs[1]; x <= cpt_IVs[ wset_lengths[ cpt_iterator[1] ] ]; x++) { // for all individual IV-cols in row wset
					
				for (y = cpt_IVs[ wset_lengths[ cpt_iterator[1] ]+1 ]; 
					y <= cpt_IVs[ wset_lengths[ cpt_iterator[1] ]+wset_lengths[ cpt_iterator[2] ] ]; y++) {	// for all individual IV-cols in col wset
					
					// which IV-cols are not in the comparison; used to sort 
					sort_omit_xy = 
						selectindex(
							((1..cols(fitstat_increments)):!=x):&(
								(1..cols(fitstat_increments)):!=y) )
					
					sort_omit_xy // ~~
						
					row_IVsort = 
						order(
							fitstat_increments[row_rownumbers, ]:>0, 
								(sort_omit_xy, y))
								
					row_IVsort // ~~
								
					row_IVsort_numbers = 
						selectindex( 
							( (fitstat_increments[
								row_rownumbers, ]:>0)[row_IVsort, y]):==0 )
						
					row_IVsort_numbers // ~~
					
					col_IVsort = 
						order(
							fitstat_increments[col_rownumbers, ]:>0, 
								(sort_omit_xy, x))
								
					col_IVsort // ~~
								
					col_IVsort_numbers = 
						selectindex( 
							( (fitstat_increments[
								col_rownumbers, ]:>0)[col_IVsort, x]):==0 )
						
					col_IVsort_numbers // ~~
								
					( ( fitstat_increments[row_rownumbers,] 
						)[row_IVsort,] 
							)[row_IVsort_numbers, ] // ~~
					
					( ( fitstat_increments[col_rownumbers,] 
						)[col_IVsort,] 
							)[col_IVsort_numbers, ] // ~~
					
					complete_dominance[x, y] = 
						sum(
							( ( ( fitstat_increments[row_rownumbers,] 
								)[row_IVsort,] 
									)[row_IVsort_numbers, x]):>(
										( ( fitstat_increments[col_rownumbers,] 
											)[col_IVsort,] 
												)[col_IVsort_numbers, y] ) )
																
					complete_dominance[x, y] = 
						complete_dominance[x, y]/(2^(length(wset_lengths)-2))
											
					complete_dominance[y, x] = 
						sum(
							( ( ( fitstat_increments[row_rownumbers,] 
								)[row_IVsort,] 
									)[row_IVsort_numbers, x]):<(
										( ( fitstat_increments[col_rownumbers,] 
											)[col_IVsort,] 
												)[col_IVsort_numbers, y] ) )
											
					complete_dominance[y, x] = 
						complete_dominance[y, x]/(2^(length(wset_lengths)-2))
						
				}
					
			}

			cpt_iterator = selectindex( cvpermute(cpt_permute_container) )
		
		}
		
		if ( any(wset_lengths:>1) ) {
			
			"within!" // ~~
			
			cpt_wsets = selectindex(wset_lengths:>1)
			
			cpt_wsets // ~~
			
			for (x = 1; x <= length(cpt_wsets); x++) {
			
				cpt_setup_vec = 
					( J(2, 1, 1) \ J( wset_lengths[ cpt_wsets[x] ]-2, 1, 0) ) 
			
				cpt_permute_container = cvpermutesetup(cpt_setup_vec) 
				
				cpt_iterator = cvpermute(cpt_permute_container)
				
				while ( length(cpt_iterator) ) {
					
					cpt_iterator // ~~
					
					cpt_IVs = 
						select(
							selectindex(wset_reps:==cpt_wsets[x]), 
								cpt_iterator )
								
					cpt_IVs // ~~
					
					// need to sort these values like above to ensure they align
					row_rownumbers = 
						selectindex( 
								(IV_indicator_matrix[, cpt_IVs[1]+1]:==1 ):&(
									IV_indicator_matrix[, cpt_IVs[2]+1]:==0 ) )
								
					col_rownumbers = 
						selectindex( 
								(IV_indicator_matrix[, cpt_IVs[1]+1]:==0 ):&(
									IV_indicator_matrix[, cpt_IVs[2]+1]:==1 ) )
								
					IV_indicator_matrix[row_rownumbers,], IV_indicator_matrix[col_rownumbers,] // ~~
					
					fitstat_increments[row_rownumbers,], fitstat_increments[col_rownumbers,] // ~~
					
					complete_dominance[cpt_IVs[1], cpt_IVs[2]] = 
						sum(
							fitstat_increments[
								row_rownumbers, cpt_IVs[1]]:>fitstat_increments[
									col_rownumbers, cpt_IVs[2]]
									)/rows(row_rownumbers)
									
					complete_dominance[cpt_IVs[2], cpt_IVs[1]] = 
						sum(
							fitstat_increments[
								row_rownumbers, cpt_IVs[1]]:<fitstat_increments[
									col_rownumbers, cpt_IVs[2]]
									)/rows(row_rownumbers)
								
					cpt_iterator = cvpermute(cpt_permute_container)
					
				}
							
			}
				
		}
		
		"complete_dominance" // ~~
		complete_dominance // ~~
	
		st_matrix("r(cptdom)", complete_dominance) 
	
	}
	
}

end

**# Mata function to compute combinations
version 12.1

mata: 

mata set matastrict on

real matrix domin_combn(real scalar number_of_IVs) {
	
	real scalar x
	
	real matrix IV_indicator_matrix, binary_pattern
	
	IV_indicator_matrix = J(number_of_IVs, 2^number_of_IVs, .)

	for (x = 1; x <= rows(IV_indicator_matrix); x++) {	//fills in the IV indicators matrix - for each row...
		
		binary_pattern = J(1, 2^(x-1), 0), J(1, 2^(x-1), 1)	//...make a binary matrix that grows exponentially; (0,1), then (0,0,1,1), then (0,0,0,0,1,1,1,1) growing in size until it fills the entire set of columns with sequences of 0's and 1's...
			
		IV_indicator_matrix[x, .] = J(1, 2^(number_of_IVs-x), binary_pattern)	//...spread the binary pattern across all rows forcing it to repeat when not long enough to fill all columns - net effect is staggering all binary patters across rows to obtain all subsets in the final matrix
			
	}
	
	return(IV_indicator_matrix)
		
}
	
end

**# Mata function to execute 'domin-flavored' models
version 12.1

mata:

	mata set matastrict on

	real scalar domin_call(string scalar IVs_in_model,  
		real scalar all_subsets_fitstat, real scalar constant_model_fitstat, ///
		struct domin_specs scalar model_specs) 
	{ 

		real scalar fitstat

		if (strlen(model_specs.mi) == 0) { //if not multiple imputation, then regular regression

			stata(model_specs.reg + " " + model_specs.dv + " " + ///
				model_specs.all + " " + IVs_in_model + " [" + ///
				model_specs.weight + model_specs.exp + "] if " + 
				model_specs.touse + ", " + model_specs.regopts, 1) //conduct regression

			fitstat = st_numscalar(st_local("fitstat")) - all_subsets_fitstat - constant_model_fitstat //record fitstat omitting constant and all subsets values; note that the fitstat to be pulled from Stata is stored as the Stata local "fitstat"

		}

		else { //otherwise, regression with "mi estimate:"

			stata("mi estimate, saving(\`mifile', replace) \`miopt': " + ///
				model_specs.reg + " " + model_specs.dv + " " + ///
				model_specs.all + " " + IVs_in_model + " [" + ///
				model_specs.weight + model_specs.exp + "] if " + 
				model_specs.touse + ", " + model_specs.regopts, 1) //conduct regression with "mi estimate:"

			stata("mi_dom, name(\`mifile') fitstat(\`fitstat') list(\`=e(m_est_mi)')", 1) //use built-in program to obtain average fitstat across imputations

			fitstat = st_numscalar("r(passstat)") - all_subsets_fitstat - constant_model_fitstat //record fitstat omitting constant and "all" subsets values with "mi estimate:"

		}

		return(fitstat)

	}
	
end

**# Mata function to execute epsilon-based relative importance
version 12.1

mata: 

mata set matastrict on

void eps_ri(string scalar varlist, string scalar reg, string scalar touse, string scalar regopts) 
{
	/*object declarations*/
	real matrix X, L, R, Lm, glmdat, glmdats, L2, R2, orth

	real rowvector V, Bt, V2, glmwgts
	
	string rowvector orthnames
	
	real scalar sd_yhat, cor_yhat
	
	string scalar predmu
	
	/*begin processing*/
	X = correlation(st_data(., tokens(varlist), st_varindex(touse))) //obtain correlations
	
	L = R = X[2..rows(X), 2..cols(X)] //set-up for svd()
	
	V = J(1, cols(X)-1, .) //placeholder for eigenvalues
	
	svd(X[2..rows(X), 2..cols(X)], L, V, R) //conduct singular value decomposition
	
	Lm = (L*diag(sqrt(V))*R) //process orthogonalized predictors
	
	if (reg == "regress") Bt = invsym(Lm)*X[2..rows(X), 1] //obtain adjusted regression weights
	
	else if (reg == "glm") { //if glm-based...
	
		glmdat = st_data(., tokens(varlist), st_varindex(touse)) //pull data into Mata
		
		L2 = R2 = glmdat[., 2..cols(glmdat)] //set-up for svd()
		
		V2 = V //placeholder for eigenvalues
		
		glmdats = (glmdat[., 2..cols(glmdat)]:-mean(glmdat[., 2..cols(glmdat)])):/sqrt(diagonal(variance(glmdat[., 2..cols(glmdat)])))' //standardize the input data
	
		svd(glmdats, L2, V2, R2) //conduct singular value decomposition on full data
		
		orth = L2*R2 //produce the re-constructed orthogonal predictors for use in regression
		
		orth = (orth:-mean(orth)):/sqrt(diagonal(variance(orth)))' //standardize the orthogonal predictors
		
		orthnames = st_tempname(cols(orth))
		
		st_addvar("double", orthnames) //generate some tempvars for Stata
		
		st_store(., orthnames, st_varindex(touse), orth) //put the orthogonalized variables in Stata
		
		stata("capture " + reg + " " + tokens(varlist)[1] + " " + invtokens(orthnames) + " if " + touse + ", " + regopts) //conduct the analysis

		if (st_numscalar("c(rc)")) {
		
			display("{err}{cmd:" + reg + "} failed when executing {cmd:epsilon}.")
		
			exit(st_numscalar("c(rc)"))
			
		}
		
		glmwgts = st_matrix("e(b)") //record the regression weights to standardize
		
		predmu = st_tempname() //generate some more tempvars for Stata
		
		sd_yhat = sqrt(variance(orth*glmwgts[1, 1..cols(glmwgts)-1]')) //SD of linear predictor
		
		stata("quietly predict double " + predmu + " if " + touse + ", mu") //translated with link function
		
		cor_yhat = correlation((glmdat[., 1], st_data(., st_varindex(predmu), st_varindex(touse))))
		
		Bt = (glmwgts[1, 1..cols(glmwgts)-1]:*((cor_yhat[2, 1])/(sd_yhat)))'

	}
	
	else { //asked for invalid reg
	
		display("{err}{opt reg(" + reg + ")} invalid with {opt epsilon}.")
	
		exit(198)
		
	}
	
	Bt = Bt:^2 //square values of regression weights
	
	Lm = Lm:^2 //square values of orthogonalized predictors

	st_matrix("r(domwgts)", (Lm*Bt)')	//produce proportion of variance explained and put into Stata
	
	st_numscalar("r(fs)", sum(Lm*Bt))	//sum relative weights to obtain R2
	
}

end
