# analysis of subjective responses to post-stimulus questions from the polystream_move experiment
# calculates correlations between subjective responses & movement responses
# fits separate LMMs to each subjective response as a function of entrance type and compares among conditions
# 09.30.2013 BH
################

# store current directory
initial.dir = getwd()

# change to directory where data reside
data_dir = '/data/polystream_move/tables/motion_analysis'
setwd(data_dir);

# Read in data & replace missing data with NA
data_fname = 'tapping_mocap_dataTbl.txt';
cat("\n\n\nData imported from the following file:",sprintf("%s/%s",data_dir,data_fname),"\n\n\n")
polymovedata = read.table(data_fname,header=T,sep='\t',na.strings=c('.'));

# ignore chameleon trials
polymovedata_loops = subset(polymovedata, !(stimulus_name %in% "chameleon.mp3"));

# for mocap analyses, ignore subjects for which we have no mocap responses
polymoveMCdata = subset(polymovedata_loops, !(subject_id %in% 
	c("08cem92181","04heg89151","09kir90191","10jzc89111")))
	
# test whether subjects' self-reported self-consciousness was any different between 1st & 2nd half of study
cat("\n\n\nDifference between self-consciousness ratings for first and second half of experiment\n")
# need to collapse responses within subjects (only 1 response per subject given, but response is repeated per trial in data table) 
collapsed_df = with(polymovedata_loops,aggregate(polymovedata_loops,
	list(subject=subject_id),mean));
# t-test for self consciousness ratings in 1st and 2nd half
selfconsc_ttest = 
with(collapsed_df,t.test(selfconscious_beginning,selfconscious_nearEnd,paired=T));
print(selfconsc_ttest);

# test whether subjects' selectivity group membership was associated with their familiarity with the term "groove" as applied to music
subSelHrdgrv_data = unique(polymovedata_loops[c("subject_id","tapping_selectivity","ever_heard_groove")])
subSelHrdgrv_data$tapping_selectivity2 = relevel(subSelHrdgrv_data$tapping_selectivity,ref="non_tapper")
cat(sprintf("%d out of %d subjects reported having heard the term groove applied to music",
	length(subSelHrdgrv_data$ever_heard_groove[subSelHrdgrv_data$ever_heard_groove=="yes"]),
	length(subSelHrdgrv_data$ever_heard_groove)))

cat("\n\nAssociation between familiarity with ther term Groove and selectivity group membership\n")
select_logtest = multinom(tapping_selectivity ~ ever_heard_groove, data=subSelHrdgrv_data)
select_logtest2 = multinom(tapping_selectivity2 ~ ever_heard_groove, data=subSelHrdgrv_data)
# logtest stats
print(summary(select_logtest))
z = summary(select_logtest)$coefficients/summary(select_logtest)$standard.errors
cat("\nz-scores associated with select_logtest (high selectivity as reference level):\n\n")
print(z)
cat("\np-values associated with select_logtest:\n\n")
print((1-pnorm(abs(z),0,1))*2)
# logtest2 stats
print(summary(select_logtest2))
z = summary(select_logtest2)$coefficients/summary(select_logtest2 )$standard.errors
cat("\nz-scores associated with select_logtest2 (non-tapper as reference level):\n\n")
print(z)
cat("\np-values associated with select_logtest2:\n\n")
print((1-pnorm(abs(z),0,1))*2)


# test differences in head movement energy among selectivity groups
cat("\n\nDifference in head movement energy among selectivity groups\n")
mcRes.sel = lmer(mean_resEnergy ~ tapping_selectivity + (1|stimulus_root) + (1|subject_id),
		data=polymoveMCdata, REML=F, na.action=na.omit)
mcRes.noSel = update(mcRes.sel, . ~ . - tapping_selectivity)
print(anova(mcRes.noSel,mcRes.sel))
print(lsmeans(mcRes.sel,pairwise~tapping_selectivity,adjust="fdr"))

##############################################################
# get correlations between subjective responses and movement responses
cat("\n\nCorrelations between subjective responses and movement responses\n\n")
# using multiple regression to get correlation of post-stim responses with movement responses while controlling for differences among subjects (following the guidelines laid out by Bland & Altman (1995) for calculating correlation coefficients with repeated observations). The code for this is repeated for each post-stim response, and should really be abstracted as a function.
###############
## TAPPING RATE
cat("\nTAPPING RATE\n\n")

# ENJOYMENT
enjoy.rmcorrtest = with(polymovedata_loops,rm_corr_test(enjoyed_music,tapping_rate,subject_id))
cat("\nEnjoyment vs. Mean Tapping Rate\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		enjoy.rmcorrtest[1],enjoy.rmcorrtest[2],enjoy.rmcorrtest[3]));

# MUSIC GROOVED
musgrooved.rmcorrtest = with(polymovedata_loops,rm_corr_test(music_grooved,tapping_rate,subject_id))
cat("\nMusic Grooved vs. Mean Tapping Rate\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		musgrooved.rmcorrtest[1],musgrooved.rmcorrtest[2],musgrooved.rmcorrtest[3]));

# LIKE TO CONTINUE
likecont.rmcorrtest = with(polymovedata_loops,rm_corr_test(like_music_to_continue,tapping_rate,subject_id))
cat("\nWanted to Continue vs. Mean Tapping Rate\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		likecont.rmcorrtest[1],likecont.rmcorrtest[2],likecont.rmcorrtest[3]));

######################
## TAPPING ENTRAINMENT
cat("\nTAPPING ENTRAINMENT\n\n")
# ENJOYMENT
enjoy_tapEnt.rmcorrtest = with(polymovedata_loops,rm_corr_test(enjoyed_music,mean_tapping_entrainmentRatio,subject_id))
cat("\nEnjoyment vs. Mean Tapping Entrainment\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		enjoy_tapEnt.rmcorrtest[1],enjoy_tapEnt.rmcorrtest[2],enjoy_tapEnt.rmcorrtest[3]));

# MUSIC GROOVED
musgrv_tapEnt.rmcorrtest = with(polymovedata_loops,rm_corr_test(music_grooved,mean_tapping_entrainmentRatio,subject_id))
cat("\nMusic Grooved vs. Mean Tapping Entrainment\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		musgrv_tapEnt.rmcorrtest[1],musgrv_tapEnt.rmcorrtest[2],musgrv_tapEnt.rmcorrtest[3]))

# LIKE TO CONTINUE
likecont_tapEnt.rmcorrtest = with(polymovedata_loops,rm_corr_test(like_music_to_continue,mean_tapping_entrainmentRatio,subject_id))
cat("\nWanted to Continue vs. Mean Tapping Entrainment\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		likecont_tapEnt.rmcorrtest[1],likecont_tapEnt.rmcorrtest[2],likecont_tapEnt.rmcorrtest[3]))

############
## MOCAP RMS
cat("\nMOCAP RMS\n\n")
# ENJOYMENT
enjoyResNrg.rmcorrtest = with(polymoveMCdata,rm_corr_test(enjoyed_music,mean_resEnergy,subject_id))
cat("\nEnjoyment vs. Mean Mocap Energy\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		enjoyResNrg.rmcorrtest[1],enjoyResNrg.rmcorrtest[2],enjoyResNrg.rmcorrtest[3]))

# MUSIC GROOVED
musgrvResNrg.rmcorrtest = with(polymoveMCdata,rm_corr_test(music_grooved,mean_resEnergy,subject_id))
cat("\nMusic Grooved vs. Mean Mocap Energy\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		musgrvResNrg.rmcorrtest[1],musgrvResNrg.rmcorrtest[2],musgrvResNrg.rmcorrtest[3]))

# LIKE TO CONTINUE
likecontResNrg.rmcorrtest = with(polymoveMCdata,rm_corr_test(like_music_to_continue,mean_resEnergy,subject_id))
cat("\nWanted to Continue vs. Mean Mocap Energy\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		likecontResNrg.rmcorrtest[1],likecontResNrg.rmcorrtest[2],likecontResNrg.rmcorrtest[3]))

####################
## MOCAP ENTRAINMENT

cat("\nMOCAP ENTRAINMENT\n\n")
# ENJOYMENT
enjoyMCEnt.rmcorrtest = with(polymoveMCdata,rm_corr_test(enjoyed_music,mean_entrainmentRatio,subject_id))
cat("\nEnjoyment vs. Mean Mocap Entrainment\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		enjoyMCEnt.rmcorrtest[1],enjoyMCEnt.rmcorrtest[2],enjoyMCEnt.rmcorrtest[3]));

# MUSIC GROOVED
musgrvMCEnt.rmcorrtest = with(polymoveMCdata,rm_corr_test(music_grooved,mean_entrainmentRatio,subject_id))
cat("\nMusic Grooved vs. Mean Mocap Entrainment\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		musgrvMCEnt.rmcorrtest[1],musgrvMCEnt.rmcorrtest[2],musgrvMCEnt.rmcorrtest[3]));

# LIKE TO CONTINUE
likecontMCEnt.rmcorrtest = with(polymoveMCdata,rm_corr_test(like_music_to_continue,mean_entrainmentRatio,subject_id))
cat("\nWanted to Continue vs. Mean Mocap Entrainment\n",
	sprintf("Correlation coefficient: %1.2f		df = %d		p-value: %1.5f\n",
		likecontMCEnt.rmcorrtest[1],likecontMCEnt.rmcorrtest[2],likecontMCEnt.rmcorrtest[3]));

##########################################################
# compare differences in subjective responses for stimultaneous & staggered conditions

cat("\n\n\nComparisons between subjective responses for simultaneous and staggered conditions, as modeled with lmer and compared with glht . . .\n")
# model subjective responses as function of entrance type & compare simult - stag conditions

# ENJOYMENT
enjoy_entrtype.lmm = lmer(enjoyed_music ~ entrance_type_stagCollapsed + (1|stimulus_root) + (1|subject_id),
	data = polymovedata_loops,REML=T,na.action=na.omit);
cat("\nEnjoyment\n")
print(lsmeans(enjoy_entrtype.lmm,pairwise~entrance_type_stagCollapsed,adjust="fdr"))

# MUSIC GROOVED
musgrv_entrtype.lmm = lmer(music_grooved ~ entrance_type_stagCollapsed + (1|stimulus_root) + (1|subject_id),
	data = polymovedata_loops,REML=T,na.action=na.omit);
cat("\nMusic Grooved\n")
print(lsmeans(musgrv_entrtype.lmm,pairwise~entrance_type_stagCollapsed,adjust="fdr"))

# WANTED MUSIC TO CONTINUE
wantcont_entrtype.lmm = lmer(like_music_to_continue ~ entrance_type_stagCollapsed + (1|stimulus_root) + (1|subject_id),
	data = polymovedata_loops,REML=T,na.action=na.omit);
cat("\nWanted Music to Continue\n")
print(lsmeans(wantcont_entrtype.lmm,pairwise~entrance_type_stagCollapsed,adjust="fdr"))

# return to original directory
setwd(initial.dir);