# Calculates a correlation coefficient and p-value for variables with repeated observations. 
# The correlation coefficient is calculated following the guidelines introduced by 
# Bland & Altman (1995) "Calculating correlation coefficients with repeated observations: 
# Part 1--correlation within subjects." Specifically, this method of correlation infers whether 
# an increase in one variable within an individual is associated with an increase in another 
# variable. Users are encouraged to refer to Bland & Altman's (1995) BMJ paper for details.
#
# 16 Oct 2013 - written by Brian Hurley
# 07 Feb 2014 - included degrees of freedom in output (BH)

# INPUT
# 	var1, var2 = vectors containing values of the first and second continuous variables that 
#		         will be tested for correlation
# 	subject_id = vector of subject identifiers (values can be numeric or strings)
# 	NOTE: it is important that the subject variable be entered last

# OUTPUT
#	corr.coef = Pearson correlation coefficient
# 	corr.pval = p-value indicating significance of correlation coefficient
#	corr.df = degrees of freedom

rm_corr_test = function(var1,var2,subject_id) {
	
	# make sure that subject_id is a factor
	if(!is.factor(subject_id)) {as.factor(subject_id)};
	
	# get degrees of freedom (N-2)
	corr.df = length(unique(subject_id))-2;
	
	# model data as anova and regression objects
	mod.aov = aov(var2~subject_id+var1);
	# regression
	mod.lm = lm(var2~subject_id+var1);
	
	# extract sum of squares
	SS.pred = (summary(mod.aov))[[1]]$'Sum Sq'[2];
	SS.resid = (summary(mod.aov))[[1]]$'Sum Sq'[3];
	
	# calculate the magnitude of the correlation coefficient within a given subject by removing between-subject variation
	# note: the ratio below is akin to the coefficient of determination. Thus, taking the sqaure root is akin to the Pearson 
	# correlation coefficient.
	abs_corr = sqrt(SS.pred/(SS.pred+SS.resid));
	
	# obtain the sign and p-value of the correlation by looking at those of the regression predictor
	n_coeffs = length(summary(mod.lm)$coefficients[,1]);
	if(summary(mod.lm)$coefficients[n_coeffs,1]<0) {corr.coeff = -abs_corr} else {
		corr.coeff = abs_corr};
	corr.pval = summary(mod.lm)$coefficients[n_coeffs,4];
	
	# output results
	return(c(corr.coeff,corr.df,corr.pval));
}