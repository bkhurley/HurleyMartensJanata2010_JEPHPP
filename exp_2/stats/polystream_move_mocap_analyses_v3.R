# statistical analyses of polystream_move mocap data
# mixed effects models are fit for mean reson energy
# 10.02.2013 BH - wrote script
# 02.09.2014 BH - changed post hoc contrasts to be done via lsmeans rather than 
#                 multcomp (easier for polynomial contrasts & cleaner code)
################


# store current directory
initial.dir = getwd()
# change to directory where data reside
data_dir = '/data/polystream_move/tables/motion_analysis';
setwd(data_dir);

# Read in data
data_fname = 'postEntrance_tapping_mocap_dataTbl.txt';
cat("\n\n\nData imported from the following file:",sprintf("%s/%s",data_dir,data_fname),"\n\n\n")
postEntrMocap = read.table(data_fname,header=T,sep='\t',na.strings=c('.'));

# set entrance_number to be a factor
postEntrMocap$entrance_number = factor(postEntrMocap$entrance_number);

# ignore chameleon trials for this set of analyses
postEntrMocap_loops = subset(postEntrMocap, !(stimulus_name %in% "chameleon.mp3"));
# for mocap analyses, ignore subjects for which we have no mocap responses
MCbyEntranceData = subset(postEntrMocap_loops, !(subject_id %in% 
	c("08cem92181","04heg89151","09kir90191","10jzc89111")))

######################
# Linear mixed effects models
######################

cat("\n\n\nEvaluating fixed effects for mean_resEnergy . . .\n")

cat("\n\nLMM model for HIGH SELECTIVE group\n\n")
postEntrResNrg.highSel = subset(MCbyEntranceData, (tapping_selectivity %in% "high_selectivity"));
mcEnt.highSel.int = lmer(postEntrance_mean_resEnergy ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number + experiment_block,
		data=postEntrResNrg.highSel, REML=FALSE, na.action=na.omit);
mcEnt.highSel.noInt = update(mcEnt.highSel.int, .~. -entrance_type_stagCollapsed:entrance_number);
print(anova(mcEnt.highSel.noInt,mcEnt.highSel.int)); # 2way int not significant
# test main effects
mcEnt.highSel.noEntrType = update(mcEnt.highSel.noInt, .~. -entrance_type_stagCollapsed);
mcEnt.highSel.noEntrNum = update(mcEnt.highSel.noInt, .~. -entrance_number);
# entrance type
print(anova(mcEnt.highSel.noEntrType,mcEnt.highSel.noInt));
# entrance number
print(anova(mcEnt.highSel.noEntrNum,mcEnt.highSel.noInt));
cat("\n\nPolynomial contrast among entrance numbers for High Selective group\n\n")
print(lsmeans(mcEnt.highSel.noInt,poly~entrance_number,adjust="fdr"))

###########################
cat("\n\nLMM model for LOW SELECTIVE group\n\n")
postEntrResNrg.lowSel = subset(MCbyEntranceData, (tapping_selectivity %in% "low_selectivity"));
mcEnt.lowSel.int = lmer(postEntrance_mean_resEnergy ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number + experiment_block,
		data=postEntrResNrg.lowSel, REML=FALSE, na.action=na.omit);
mcEnt.lowSel.noInt = update(mcEnt.lowSel.int, .~. -entrance_type_stagCollapsed:entrance_number);
print(anova(mcEnt.lowSel.noInt,mcEnt.lowSel.int)); # 2way interaction is significant

cat("\n\nContrasts among interaction conditions for Low Selective group\n\n")
print(lsmeans(mcEnt.lowSel.int,list(pairwise~entrance_type_stagCollapsed|entrance_number,
	poly~entrance_number|entrance_type_stagCollapsed),adjust="fdr"))

#########################
cat("\n\nLMM model for NON-TAPPER group\n\n")
postEntrResNrg.nonTap = subset(MCbyEntranceData, (tapping_selectivity %in% "non_tapper"));
mcEnt.nonTap.int = lmer(postEntrance_mean_resEnergy ~ (1|stimulus_root) + (1|subject_id) + 
	entrance_type_stagCollapsed*entrance_number + experiment_block,
	data=postEntrResNrg.nonTap, REML=FALSE, na.action=na.omit);
mcEnt.nonTap.noInt = update(mcEnt.nonTap.int, .~. -entrance_type_stagCollapsed:entrance_number);
print(anova(mcEnt.nonTap.noInt,mcEnt.nonTap.int)); # not significant
# test for main effects
mcEnt.nonTap.noEntrType = update(mcEnt.nonTap.noInt, .~. -entrance_type_stagCollapsed);
mcEnt.nonTap.noEntrNum = update(mcEnt.nonTap.noInt, .~. -entrance_number);
print(anova(mcEnt.nonTap.noEntrType,mcEnt.nonTap.noInt));
print(anova(mcEnt.nonTap.noEntrNum,mcEnt.nonTap.noInt)); # no significant main effects

# return to original directory
setwd(initial.dir)