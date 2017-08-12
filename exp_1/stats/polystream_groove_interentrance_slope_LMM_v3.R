# statistical analyses of polystream_groove slider data 
# fit model to inter-entrance slider slopes
# 07.12.2013 BH
################
# Read in & prep data

# store current directory
initial.dir = getwd()
# change to directory where data reside
data_dir = '/data/polystream/groove_rating/tables';
setwd(data_dir);

# Read in data
data_fname = 'timevar_sliderslope_data.txt';
cat("\n\n\nData imported from the following file:",sprintf("%s/%s",data_dir,data_fname),"\n\n\n")
sliderData = read.table(data_fname,header=T,sep='\t',na.strings=c('.'));

################
# Fit linear mixed-effects models
################

# RANDOM EFFECTS

# full model
slidermod.full = lmer(interEntrance_slope ~ 
		entrance_type*entrance_interval*musically_trained + (1|excerpt) + (1|subject_id),
		data = sliderData,REML=F,na.action = na.omit);
(summary(slidermod.full));

# full model w/ no training interactions
slidermod.noTrainingInts = lmer(interEntrance_slope ~ 
		entrance_type*entrance_interval + musically_trained + (1|excerpt) + (1|subject_id),
		data = sliderData,REML=F,na.action = na.omit);
(summary(slidermod.noTrainingInts));

cat("\nTest of fit goodness with excerpt as an additional random intercept beyond subject_id\n")
# test the comparitive fit of random effects structure with (a) stimulus root ('excerpt') & subject id and (b) subject id only
slidermod.noStimRoot = update(slidermod.full,.~. - (1|excerpt));
print(anova(slidermod.noStimRoot,slidermod.full)); # fits signfiicantly better (p<.01) with excerpt & subject id as random intercepts

# EVALUATE FIXED EFFECTS
cat("\n\n\nEvaluating fixed effects . . .\n")
print(summary(slidermod.full))

# test significance of 3-way interaction
cat("\nSignificance test of entrance_type:entrance_interval:musically_trained interaction\n")
slidermod.no3wayInt = update(slidermod.full, .~. -entrance_type:entrance_interval:musically_trained)
print(anova(slidermod.no3wayInt,slidermod.full)) # 3way interaction not significant

# assess significance of interactions
cat("\nSignificance test of entrance_type:musically_trained\n")
slidermod.noEntrtypeTraining_int = update(slidermod.no3wayInt, .~. -entrance_type:musically_trained);
print(anova(slidermod.noEntrtypeTraining_int,slidermod.no3wayInt)) # not significant (p = 0.4909)

cat("\nSignificance test of entrance_interval:musically_trained\n")
slidermod.noEntrintervalTraining_int = update(slidermod.noEntrtypeTraining_int, .~. -entrance_interval:musically_trained);
print(anova(slidermod.noEntrintervalTraining_int,slidermod.noEntrtypeTraining_int)) # interaction is significant (p<.05)

cat("\nSignificance test of entrance_type:entrance_interval interaction\n")
slidermod.noEntrtypeEnterinterval_int = update(slidermod.noEntrtypeTraining_int, .~. -entrance_type:entrance_interval);
print(anova(slidermod.noEntrtypeEnterinterval_int,slidermod.noEntrtypeTraining_int)) # interaction is significant (p<.01));

# test significance of music training factor
cat("\nSignificance test of musical training term\n")
slidermod.noTraining = update(slidermod.noTrainingInts, .~. - musically_trained);
print(anova(slidermod.noTraining,slidermod.noTrainingInts)); # not significant (p = 0.1657).



# MULTIPLE COMPARISONS
cat("\n\n\nEvaluating contrasts . . .\n\n")

# test contrasts & polynomial trends for each entrance types and levels of training
print(lsmeans(slidermod.noEntrtypeTraining_int,list(pairwise~entrance_type|entrance_interval,
	poly~entrance_interval|entrance_type,
	pairwise~musically_trained|entrance_interval,
	poly~entrance_interval|musically_trained),adjust="fdr"))
	

# plot music training X entrance interval interaction
musTrained_data = subset(sliderData, (musically_trained %in% "trained"))
musNontrained_data = subset(sliderData, (musically_trained %in% "non-trained"))
entrInterval_levels = unique(sliderData$entrance_interval)
num_intervals = length(entrInterval_levels)

training_interval_mat = matrix(data=NA,nrow=2,ncol=4)
rownames(training_interval_mat) = c("Musically Trained","Not Musically Trained")
colnames(training_interval_mat) = c("1 to 2","2 to 3","3 to 4","Last to End")

# trained means
for (i_entr in 1:num_intervals) {
	this_interval_data = subset(musTrained_data, (entrance_interval %in% entrInterval_levels[i_entr]))
	training_interval_mat[1,i_entr] = mean(this_interval_data$interEntrance_slope,na.rm = T)
}

# non-trained means
for (i_entr in 1:num_intervals) {
	this_interval_data = subset(musNontrained_data, (entrance_interval %in% entrInterval_levels[i_entr]))
	training_interval_mat[2,i_entr] = mean(this_interval_data$interEntrance_slope,na.rm = T)
}

barplot(training_interval_mat,beside=T,legend=T,ylab="Inter-Entrance Slope",xlab="Entrance Intervals")

# return to original directory
setwd(initial.dir)