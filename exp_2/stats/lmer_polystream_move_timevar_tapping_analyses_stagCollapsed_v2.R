# statistical analyses of polystream_move tapping data 
# fit models to time-varying post-entrance tap rate
# 05.01.2013 BH
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
tapByEntrance = read.table(data_fname,header=T,sep='\t',na.strings=c('.'));

# set entrance_number to be a factor
tapByEntrance$entrance_number = factor(tapByEntrance$entrance_number);

# ignore chameleon trials for this analysis (i.e., AppleLoop stimuli only)
tapByEntrance_loops = subset(tapByEntrance, !(stimulus_name %in% "chameleon.mp3"));
# for mocap analyses, ignore subjects for which we have no mocap responses
MCbyEntranceData = subset(tapByEntrance_loops, !(subject_id %in% 
	c("08cem92181","04heg89151","09kir90191","10jzc89111")))

cat("\n################################################\n")

#######################
# postEntrance_tappingRate mixed-effects models
#

taprate.full = lmer(postEntrance_tappingRate ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number*tapping_selectivity + experiment_block,
		data=tapByEntrance_loops, REML=FALSE, na.action=na.omit);

## FIXED EFFECTS
cat("\n\n\nEvaluating fixed effects . . .\n")

print(summary(taprate.full))

# # test significance of music training covariate
# cat("\nSignificance test of musical training covariate term\n")
# taprate.NoTraining = update(taprate.full, .~. - years_training);
# print(anova(taprate.NoTraining,taprate.full)); # years training not significant

# # test significance of experiment block covariate
# cat("\nSignificance test of experiment block covariate term\n")
# taprate.noBlock = update(taprate.full,. ~ . -experiment_block);
# print(anova(taprate.noBlock,taprate.full)); # block highly significant (p<.001)

# # test significance of three-way interaction
# cat("\nSignificance test of entrance_type*entrance_number*tapping_selectivity interaction\n")
# taprate.no3way = update(taprate.full,.~. - entrance_type_stagCollapsed:entrance_number:tapping_selectivity);
# print(anova(taprate.no3way,taprate.full)); # 3way interaction is significant (p<.01)

cat("\n\n\nTest significance of entrance_type*entrance_number interaction with separate models for each level of 
	selectivity . . .\n")

cat("\nHIGH SELECTIVITY\n")
tapByEntrance.highSel = subset(tapByEntrance_loops, (tapping_selectivity %in% "high_selectivity"));
taprate.highSel_int = lmer(postEntrance_tappingRate ~ (1 | stimulus_root) + (1 | subject_id) + 
    entrance_type_stagCollapsed*entrance_number + experiment_block, 
    data=tapByEntrance.highSel,REML=F,na.action=na.omit);
taprate.highSel_noInt = update(taprate.highSel_int, .~. - entrance_type_stagCollapsed:entrance_number);
print(anova(taprate.highSel_noInt,taprate.highSel_int)); # entrance_type:entrance_number significant (p<.001) for high selectivity
cat("\nMultple comparisons of entrance_type:entrance_number interaction conditions\n")
print(lsmeans(taprate.highSel_int,list(pairwise~entrance_type_stagCollapsed|entrance_number,
	poly~entrance_number|entrance_type_stagCollapsed),adjust="fdr"))

cat("\nLOW SELECTIVITY\n")
tapByEntrance.lowSel = subset(tapByEntrance_loops, (tapping_selectivity %in% "low_selectivity"));
taprate.lowSel_int = lmer(postEntrance_tappingRate ~ (1 | stimulus_root) + (1 | subject_id) + 
    entrance_type_stagCollapsed*entrance_number + experiment_block, 
    data=tapByEntrance.lowSel,REML=F,na.action=na.omit);
taprate.lowSel_noInt = update(taprate.lowSel_int, .~. - entrance_type_stagCollapsed:entrance_number);
print(anova(taprate.lowSel_noInt,taprate.lowSel_int)); # entrance_type:entrance_number significant (p<.001) for low selectivity
cat("\nMultple comparisons of entrance_type within each level of entrance_number\n")
print(lsmeans(taprate.lowSel_int,list(pairwise~entrance_type_stagCollapsed|entrance_number,
	poly~entrance_number|entrance_type_stagCollapsed),adjust="fdr"))

cat("\nNon-Tappers")
tapByEntrance.nonTap = subset(tapByEntrance_loops, (tapping_selectivity %in% "non_tapper"));
taprate.nonTap_int = lmer(postEntrance_tappingRate ~ (1 | stimulus_root) + (1 | subject_id) + 
    entrance_type_stagCollapsed*entrance_number + experiment_block, 
    data=tapByEntrance.nonTap,REML=F,na.action=na.omit);
# test significance of entrance_type:entrance_number interaction
taprate.nonTap_noInt = update(taprate.nonTap_int, .~. - entrance_type_stagCollapsed:entrance_number);
print(anova(taprate.nonTap_noInt,taprate.nonTap_int)); # entrance_type:entrance_number not significant
# test main effects
taprate.nonTap_noEntrType = update(taprate.nonTap_noInt, .~. - entrance_type_stagCollapsed);
taprate.nonTap_noEntrNum = update(taprate.nonTap_noInt, .~. - entrance_number);
print(anova(taprate.nonTap_noEntrType,taprate.nonTap_noInt));
print(anova(taprate.nonTap_noEntrNum,taprate.nonTap_noInt)); # no significant effects for non_tappers

# return to original directory
setwd(initial.dir)