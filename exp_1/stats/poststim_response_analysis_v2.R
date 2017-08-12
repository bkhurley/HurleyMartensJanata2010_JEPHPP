# analysis of subjective responses to post-stimulus questions from the polystream_groove experiment
# calculates correlations between subjective responses & mean slider rating (averaged across last few s of trials)
# fits separate LMMs to each subjective response as a function of (1) entrance type & (2) instrument density, and compares among conditions
# 09.17.2013 BH
################

# store current directory
initial.dir = getwd()

# change to directory where data reside
data_dir = '/data/polystream/groove_rating/tables'
setwd(data_dir);

# Read in mean slider ratings & replace missing data with NA
data_fname = 'slider_groove_data.txt';
cat("\n\n\nData imported from the following file:",sprintf("%s/%s",data_dir,data_fname),"\n\n\n")
polydata = read.table(data_fname,header=T,sep='\t',na.strings=c('.'));

# get correlations between subjective responses and slider means
cat("\n\nCorrelations between subjective responses and slider means\n")
# use multiple regression to get correlation of post-stim responses with slider ratings while controlling for differences among subjects (following the guidelines laid out by Bland & Altman (1995) for calculating correlation coefficients with repeated observations). The code for this is repeated for each post-stim response, and should really be abstracted to a function.

# ENJOYMENT
enjoy_meanslider.rmcorrtest = with(polydata,rm_corr_test(enjoyment,late_mean,subject_id))
	cat("\nEnjoyment vs. Mean Slider rating\n",
		sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
			enjoy_meanslider.rmcorrtest[1],enjoy_meanslider.rmcorrtest[2],enjoy_meanslider.rmcorrtest[3]))

# URGE TO MOVE
urgemv_meanslider.rmcorrtest = with(polydata,rm_corr_test(urge_to_move,late_mean,subject_id))
cat("\nUrge to Move vs. Mean Slider rating\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		urgemv_meanslider.rmcorrtest[1],urgemv_meanslider.rmcorrtest[2],urgemv_meanslider.rmcorrtest[3]))

# LIKE TO CONTINUE
likecont_meanslider.rmcorrtest = with(polydata,rm_corr_test(like_to_continue,late_mean,subject_id))
cat("\nLike to Continue vs. Mean Slider rating\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		likecont_meanslider.rmcorrtest[1],likecont_meanslider.rmcorrtest[2],likecont_meanslider.rmcorrtest[3]))


# ignore instrument density condition responses
data_entrtype = subset(polydata, !(entrance_type %in% ""));
summary(data_entrtype)

# testing paired differences w/ lmer & glht to account for inter-stim variability, rather than collapsing across stims for a traditional paired t-test

cat("\n\n\nComparisons between subjective responses for simultaneous and staggered conditions, as modeled with lmer and compared with lsmeans . . .\n")
# model subjective responses as function of entrance type & compare simult - stag conditions
# enjoyment
enjoy_entrtype.lmm = lmer(enjoyment ~ entrance_type + (1|root_excerpt) + (1|subject_id),
	data = data_entrtype,REML=T,na.action=na.omit);
cat("\nEnjoyment\n")
print(lsmeans(enjoy_entrtype.lmm,pairwise~entrance_type,adjust="fdr")) # no need for adjustment since only 1 comparison

# urge to move
urgemove_entrtype.lmm = lmer(urge_to_move ~ entrance_type + (1|root_excerpt) + (1|subject_id),
	data = data_entrtype,REML=T,na.action=na.omit);
cat("\nUrge to move\n")
print(lsmeans(urgemove_entrtype.lmm,pairwise~entrance_type,adjust="fdr"))

# like to continue
likecontin.lmm = lmer(like_to_continue ~ entrance_type + (1|root_excerpt) + (1|subject_id),
	data = data_entrtype,REML=T,na.action=na.omit);
cat("\nLike to continue\n")
print(lsmeans(likecontin.lmm,pairwise~entrance_type,adjust="fdr"))

####################################################################
# now model subjective responses as a function of instrument density
cat("\n\n\nEffects of instrument density on subjective responses and comparisons\n")
# ignore instrument density condition responses
data_instdens = subset(polydata, !(instrument_density %in% ""));
summary(data_instdens);

# enjoyment
enjoy_instdens.full = lmer(enjoyment ~ instrument_density + (1|root_excerpt) + (1|subject_id),
	data = data_instdens,REML=T,na.action=na.omit);
enjoy_instdens.restrict = update(enjoy_instdens.full, .~. -instrument_density);
cat("\nEnjoyment\n")
print(anova(enjoy_instdens.restrict,enjoy_instdens.full));
print(lsmeans(enjoy_instdens.full,pairwise~instrument_density,adjust="fdr"))

# urge to move
urgemove_instdens.full = lmer(urge_to_move ~ instrument_density + (1|root_excerpt) + (1|subject_id),
	data = data_instdens,REML=T,na.action=na.omit);
urgemove_instdens.restrict = update(urgemove_instdens.full,.~. - instrument_density);
cat("\nUrge to move\n")
print(anova(urgemove_instdens.restrict,urgemove_instdens.full));
print(lsmeans(urgemove_instdens.full,pairwise~instrument_density,adjust="fdr"))

# like to continue
likecontin_instdens.full = lmer(like_to_continue ~ instrument_density + (1|root_excerpt) + (1|subject_id),
	data = data_instdens,REML=T,na.action=na.omit);
likecontin_instdens.restrict = update(likecontin_instdens.full,.~. -instrument_density);
cat("\nLike to continue\n")
print(anova(likecontin_instdens.restrict,likecontin_instdens.full));
print(lsmeans(likecontin_instdens.full,pairwise~instrument_density,adjust="fdr"))

# return to original directory
setwd(initial.dir);