# statistical analyses of polystream_move mocap data
# mixed effects models are fit for entrainment ratio
# 10.02.2013 BH
################

# store current directory
initial.dir = getwd()
# change to directory where data reside
data_dir = '/data/polystream_move/tables/motion_analysis'
setwd(data_dir)

# Read in data
data_fname = 'postEntrance_tapping_mocap_dataTbl.txt'
cat("\n\n\nData imported from the following file:",sprintf("%s/%s",data_dir,data_fname),"\n\n\n")
postEntrMocap = read.table(data_fname,header=T,sep='\t',na.strings=c('.'))

# set entrance_number to be a factor
postEntrMocap$entrance_number = factor(postEntrMocap$entrance_number)

# ignore chameleon trials for this set of analyses
postEntrMocap_loops = subset(postEntrMocap, !(stimulus_name %in% "chameleon.mp3"))
# for mocap analyses, ignore subjects for which we have no mocap responses
MCbyEntranceData = subset(postEntrMocap_loops, !(subject_id %in% 
	c("08cem92181","04heg89151","09kir90191","10jzc89111")))
	
######################
# Linear mixed effects models
######################
cat("\n\n\nEvaluating fixed effects for SMR . . .\n")

# mocapEnt.full = lmer(postEntrance_mean_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
		# entrance_type_stagCollapsed*entrance_number*tapping_selectivity + years_training + experiment_block,
		# data=MCbyEntranceData, REML=FALSE, na.action=na.omit);

# # test significance of music training covariate
# cat("\nSignificance test of musical training\n")
# mocapEnt.NoTraining = update(mocapEnt.full, .~. - years_training);
# print(anova(mocapEnt.NoTraining,mocapEnt.full));

cat("\n\nLMM model for HIGH SELECTIVE group\n\n")
postEntr_mcEnt.highSel = subset(MCbyEntranceData, (tapping_selectivity %in% "high_selectivity"));
mocapEnt.highSel.int = lmer(postEntrance_mean_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number + experiment_block,
		data=postEntr_mcEnt.highSel, REML=FALSE, na.action=na.omit);
mocapEnt.highSel.noInt = update(mocapEnt.highSel.int, .~. -entrance_type_stagCollapsed:entrance_number);
print(anova(mocapEnt.highSel.noInt,mocapEnt.highSel.int)); # 2way int not significant
# test main effects
mocapEnt.highSel.noEntrType = update(mocapEnt.highSel.noInt, .~. -entrance_type_stagCollapsed);
mocapEnt.highSel.noEntrNum = update(mocapEnt.highSel.noInt, .~. -entrance_number);
# entrance type
print(anova(mocapEnt.highSel.noEntrType,mocapEnt.highSel.noInt));
# entrance number
print(anova(mocapEnt.highSel.noEntrNum,mocapEnt.highSel.noInt));

cat("\n\nPolynomial contrast across entrance numbers for High Selective group\n\n")
print(lsmeans(mocapEnt.highSel.noInt,poly~entrance_number,adjust="fdr"))

###########################
cat("\n\nLMM model for LOW SELECTIVE group\n\n")
postEntr_mcEnt.lowSel = subset(MCbyEntranceData, (tapping_selectivity %in% "low_selectivity"));
mocapEnt.lowSel.int = lmer(postEntrance_mean_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number + experiment_block,
		data=postEntr_mcEnt.lowSel, REML=FALSE, na.action=na.omit);
mocapEnt.lowSel.noInt = update(mocapEnt.lowSel.int, .~. -entrance_type_stagCollapsed:entrance_number);
print(anova(mocapEnt.lowSel.noInt,mocapEnt.lowSel.int)); # interaction is significant

cat("\n\nContrasts among interaction conditions for Low Selective group\n\n")
print(lsmeans(mocapEnt.lowSel.int,list(pairwise~entrance_type_stagCollapsed|entrance_number,
	poly~entrance_number|entrance_type_stagCollapsed),adjust="fdr"))

#########################
cat("\n\nLMM model for NON-TAPPER group\n\n")
postEntr_mcEnt.nonTap = subset(MCbyEntranceData, (tapping_selectivity %in% "non_tapper"));
mocapEnt.nonTap.int = lmer(postEntrance_mean_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
	entrance_type_stagCollapsed*entrance_number + experiment_block,
	data=postEntr_mcEnt.nonTap, REML=FALSE, na.action=na.omit);
mocapEnt.nonTap.noInt = update(mocapEnt.nonTap.int, .~. -entrance_type_stagCollapsed:entrance_number);
print(anova(mocapEnt.nonTap.noInt,mocapEnt.nonTap.int)); # interaction not significant
# test main effects
mocapEnt.nonTap.noEntrType = update(mocapEnt.nonTap.noInt, .~. -entrance_type_stagCollapsed);
mocapEnt.nonTap.noEntrNum = update(mocapEnt.nonTap.noInt, .~. -entrance_number);
# entrance type
print(anova(mocapEnt.nonTap.noEntrType,mocapEnt.nonTap.noInt)); # not significant
# entrance number
print(anova(mocapEnt.nonTap.noEntrNum,mocapEnt.nonTap.noInt)); # entrance number is significant

cat("\n\nPolynomial contrasts across entrance numbers for Non-Tapper group\n\n")
print(lsmeans(mocapEnt.nonTap.noInt,poly~entrance_number,adjust="fdr"))

# return to original directory
setwd(initial.dir)