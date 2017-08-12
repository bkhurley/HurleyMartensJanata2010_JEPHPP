# statistical analyses of polystream_move tapping entrainment 
# fit models to post-entrance tapping entrainment
# 08.02.2013 BH
################
# Read in & prep data

# store current directory
initial.dir = getwd()
# change to directory where data reside
data_dir = '/data/polystream_move/tables/motion_analysis';
setwd(data_dir);

# Read in data
data_fname = 'postEntrance_tapping_mocap_dataTbl.txt';
cat("\n\n\nData imported from the following file:",sprintf("%s/%s",data_dir,data_fname),"\n\n\n")
postEntrTap = read.table(data_fname,header=T,sep='\t',na.strings=c('.'));

# set entrance_number to be a factor
postEntrTap$entrance_number = factor(postEntrTap$entrance_number);

# ignore chameleon trials for this set of analyses
postEntrTap_loops = subset(postEntrTap, !(stimulus_name %in% "chameleon.mp3"));

#######################
# postEntrance_tapEntrainmentRatio mixed-models
#

## FIXED EFFECTS
cat("\n\n\nEvaluating fixed effects . . .\n")
tapEnt.full = lmer(postEntrance_tapping_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number*tapping_selectivity + experiment_block,
		data=postEntrTap_loops, REML=FALSE, na.action=na.omit);

# # test significance of music training covariate
# cat("\nSignificance test of musical training covariate term\n")
# tapEnt.NoTraining = update(tapEnt.full, .~. - years_training);
# print(anova(tapEnt.NoTraining,tapEnt.full)); # years_training does not significantly contribute to model fit

cat("\n\nLMM model for HIGH SELECTIVE group\n\n")
postEntrTap.highSel = subset(postEntrTap_loops, (tapping_selectivity %in% "high_selectivity"));
tapEnt.highSel.int = lmer(postEntrance_tapping_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number + experiment_block,
		data=postEntrTap.highSel, REML=FALSE, na.action=na.omit);
tapEnt.highSel.noInt = update(tapEnt.highSel.int, .~. -entrance_type_stagCollapsed:entrance_number);
print(anova(tapEnt.highSel.noInt,tapEnt.highSel.int)); # 2-way int not significant

tapEnt.highSel.noEntrType = update(tapEnt.highSel.noInt, .~. -entrance_type_stagCollapsed);
tapEnt.highSel.noEntrNum = update(tapEnt.highSel.noInt, .~. -entrance_number);
# test main effect of entrance type
print(anova(tapEnt.highSel.noEntrType,tapEnt.highSel.noInt));
# test main effect of entrance number
print(anova(tapEnt.highSel.noEntrNum,tapEnt.highSel.noInt));
cat("\n\nPolynomial contrast of entrance_number for Highly Selective group\n\n")
print(lsmeans(tapEnt.highSel.noInt,poly~entrance_number,adjust="fdr"))

##################

cat("\n\nLMM model for LOW SELECTIVE group\n\n")
postEntrTap.lowSel = subset(postEntrTap_loops, (tapping_selectivity %in% "low_selectivity"));
tapEnt.lowSel.int = lmer(postEntrance_tapping_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number + experiment_block,
		data=postEntrTap.lowSel, REML=FALSE, na.action=na.omit);
tapEnt.lowSel.noInt = update(tapEnt.lowSel.int, .~. -entrance_type_stagCollapsed:entrance_number);
print(anova(tapEnt.lowSel.noInt,tapEnt.lowSel.int)); # 2way int not significant
# test main effects
tapEnt.lowSel.noEntrType = update(tapEnt.lowSel.noInt, .~. -entrance_type_stagCollapsed);
tapEnt.lowSel.noEntrNum = update(tapEnt.lowSel.noInt, .~. -entrance_number);
print(anova(tapEnt.lowSel.noEntrType,tapEnt.lowSel.noInt));
print(anova(tapEnt.lowSel.noEntrNum,tapEnt.lowSel.noInt));
cat("\n\nPolynomial constrast across entrance numbers for Low Selective group\n\n")
# no need to waste memory on entrance_type comparison since only 2 levels
print(lsmeans(tapEnt.lowSel.noInt,poly~entrance_number,adjust="fdr"))

##################

cat("\n\nLMM model for NON-TAPPER group\n\n")
postEntrTap.nonTap = subset(postEntrTap_loops, (tapping_selectivity %in% "non_tapper"));
tapEnt.nonTap.int = lmer(postEntrance_tapping_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number + experiment_block,
		data=postEntrTap.nonTap, REML=FALSE, na.action=na.omit);
tapEnt.nonTap.noInt = update(tapEnt.nonTap.int, .~. -entrance_type_stagCollapsed:entrance_number);
print(anova(tapEnt.nonTap.noInt,tapEnt.nonTap.int)); # significant 2way int
cat("\n\nMultiple comparisons among entrance_types:entrance_number interaction conditions for Non-Tapping group\n\n")
print(lsmeans(tapEnt.nonTap.int,list(pairwise~entrance_type_stagCollapsed|entrance_number,
	poly~entrance_number|entrance_type_stagCollapsed),adjust="fdr"))

# return to original directory
setwd(initial.dir)