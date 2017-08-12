# Assess effects of musical training on the various dependent measures

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
	
cat("\n\nSUBJECT DESCRIPTOR STATS\n")
subSelTrain_data = unique(tapByEntrance_loops[c("subject_id","tapping_selectivity",
	"musically_trained","years_training")])
print(subSelTrain_data)
trained_subs = subSelTrain_data[subSelTrain_data$musically_trained=="trained",]
untrained_subs = subSelTrain_data[subSelTrain_data$musically_trained=="non-trained",]
cat(sprintf("\n%d out of %d subjects were musically trained.\n",
	length(trained_subs[,1]),length(subSelTrain_data[,1])))
cat(sprintf("\nTrained subs yrs trained: M = %1.2f, SD = %1.2f\n",
	mean(trained_subs[,"years_training"]),sd(trained_subs[,"years_training"])))
cat(sprintf("\nUntrained subs yrs trained: M = %1.2f, SD = %1.2f\n",
	mean(untrained_subs[,"years_training"]),sd(untrained_subs[,"years_training"])))
trainingSelectivity_tbl = with(subSelTrain_data,table(musically_trained,tapping_selectivity))
print(trainingSelectivity_tbl)
		
cat("\n################################################\n")

cat("\nEFFECTS OF MUSICIANSHIP ON THE FOLLOWING MEASURES:\n")

# cat("\n\nMOVEMENT SELECTIVITY\n")
# require(nnet)
# select_logtest = multinom(tapping_selectivity ~ musically_trained, data=subSelTrain_data)
# summary(select_logtest)

cat("\n\nTAP RATE\n")
cat("\n\nTest of 3-way interaction (entrance type X entrance period X musically trained)\n")
taprate.full = lmer(postEntrance_tappingRate ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number*musically_trained + experiment_block,
		data=tapByEntrance_loops, REML=FALSE, na.action=na.omit);

taprate.no3wayInt = update(taprate.full, .~. -entrance_type_stagCollapsed:entrance_number:musically_trained)
print(anova(taprate.no3wayInt,taprate.full))

# test significance of 2-way interactions involving musical training
cat("\n\nEntrance type X musically trained\n")
taprate.noEntrTypeMusTrain_int = update(taprate.no3wayInt, .~. -entrance_type_stagCollapsed:musically_trained)
print(anova(taprate.noEntrTypeMusTrain_int,taprate.no3wayInt)) # significant interaction (p<.01)

cat("\n\nEntrance period X musically trained\n")
taprate.noEntrNum_musTrain_int = update(taprate.no3wayInt, .~. -entrance_number:musically_trained)
print(anova(taprate.noEntrNum_musTrain_int,taprate.no3wayInt)) # significant interaction (p<.0001)

# test significance of entrance type X entrance period interaction
cat("\n\nEntrance type X Entrance period\n") 
taprate.noEntrType_EntrNum_int = update(taprate.no3wayInt, .~. -entrance_type_stagCollapsed:entrance_number)
print(anova(taprate.noEntrType_EntrNum_int,taprate.no3wayInt)) # highly significant (p<.0001)

# main effect of musical training
cat("\n\nMusically Trained main effect\n")
taprate.noMusTrain_ints = update(taprate.no3wayInt, .~. -entrance_type_stagCollapsed:musically_trained 
	-entrance_number:musically_trained)
taprate.noMusTrain = update(taprate.noMusTrain_ints, .~. -musically_trained)
print(anova(taprate.noMusTrain,taprate.noMusTrain_ints)) # not significant (p = .094)


cat("\n\nTap Rate Contrasts\n")
print(lsmeans(taprate.no3wayInt,list(
	pairwise~entrance_type_stagCollapsed|musically_trained,
	poly~entrance_number|musically_trained,
	pairwise~entrance_type_stagCollapsed|entrance_number,
	poly~entrance_number|entrance_type_stagCollapsed
	),adjust="fdr"))
	
#############	
cat("\n\nTAP SYNCHRONIZATION\n")
cat("\n3-way interaction (entrance type X entrance period X musically trained) on SMR\n")
tapSMR.full = lmer(postEntrance_tapping_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number*musically_trained + experiment_block,
		data=tapByEntrance_loops, REML=FALSE, na.action=na.omit);
tapSMR.no3wayInt = update(tapSMR.full, .~. -entrance_type_stagCollapsed:entrance_number:musically_trained)
print(anova(tapSMR.no3wayInt,tapSMR.full)) # not significant

# test significance of 2-way interactions involving musical training
cat("\nEntrance type X musically trained\n")
tapSMR.noEntrTypeMusTrain_int = update(tapSMR.no3wayInt, .~. -entrance_type_stagCollapsed:musically_trained)
print(anova(tapSMR.noEntrTypeMusTrain_int,tapSMR.no3wayInt)) # significant (p<.05)

cat("\nEntrance period X musically trained\n")
tapSMR.noEntrNum_musTrain_int = update(tapSMR.no3wayInt, .~. -entrance_number:musically_trained)
print(anova(tapSMR.noEntrNum_musTrain_int,tapSMR.no3wayInt)) # significant (p<.05)

# test significance of entrance type X entrance period interaction
cat("\n\nEntrance type X Entrance period\n") 
tapSMR.noEntrType_EntrNum_int = update(tapSMR.no3wayInt, .~. -entrance_type_stagCollapsed:entrance_number)
print(anova(tapSMR.noEntrType_EntrNum_int,tapSMR.no3wayInt)) # not significant

# main effect of musical training
cat("\n\nMusically Trained main effect\n")
tapSMR.noMusTrain_ints = update(tapSMR.no3wayInt, .~. -entrance_type_stagCollapsed:musically_trained 
	-entrance_number:musically_trained)
tapSMR.noMusTrain = update(tapSMR.noMusTrain_ints, .~. -musically_trained)
print(anova(tapSMR.noMusTrain,tapSMR.noMusTrain_ints)) # significant 

cat("\nInteraction contrasts\n")
print(lsmeans(tapSMR.no3wayInt,list(pairwise~entrance_type_stagCollapsed|musically_trained,
	poly~entrance_number|musically_trained),adjust="fdr"))

#############
cat("\n\nMOCAP ENERGY\n")

cat("\n3-way interaction (entrance type X entrance period X musically trained)\n")
mcEnergy.full = lmer(postEntrance_mean_resEnergy ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number*musically_trained + experiment_block,
		data=tapByEntrance_loops, REML=FALSE, na.action=na.omit);
mcEnergy.no3wayInt = update(mcEnergy.full, .~. -entrance_type_stagCollapsed:entrance_number:musically_trained)
print(anova(mcEnergy.no3wayInt,mcEnergy.full)) # not significant

# test significance of 2-way interactions involving musical training
cat("\n\nEntrance type X musically trained\n")
mcEnergy.noEntrTypeMusTrain_int = update(mcEnergy.no3wayInt, .~. -entrance_type_stagCollapsed:musically_trained)
print(anova(mcEnergy.noEntrTypeMusTrain_int,mcEnergy.no3wayInt)) # not significant (p=.067)

cat("\n\nEntrance period X musically trained\n")
mcEnergy.noEntrNum_musTrain_int = update(mcEnergy.no3wayInt, .~. -entrance_number:musically_trained)
print(anova(mcEnergy.noEntrNum_musTrain_int,mcEnergy.no3wayInt)) # not significant

# test significance of entrance type X entrance period interaction
cat("\n\nEntrance type X Entrance period\n") 
mcEnergy.noEntrType_EntrNum_int = update(mcEnergy.no3wayInt, .~. -entrance_type_stagCollapsed:entrance_number)
print(anova(mcEnergy.noEntrType_EntrNum_int,mcEnergy.no3wayInt)) # not significant (but almost; p = .053)

# main effect of musical training
cat("\n\nMusically Trained main effect\n")
mcEnergy.noMusTrain_ints = update(mcEnergy.no3wayInt, .~. -entrance_type_stagCollapsed:musically_trained 
	-entrance_number:musically_trained)
mcEnergy.noMusTrain = update(mcEnergy.noMusTrain_ints, .~. -musically_trained)
print(anova(mcEnergy.noMusTrain,mcEnergy.noMusTrain_ints)) # not significant


#############
cat("\n\nMOCAP SYNCHRONIZATION\n")

cat("\n\n3-way interaction (entrance type X entrance period X musically trained)\n")
mcSMR.full = lmer(postEntrance_mean_entrainmentRatio ~ (1|stimulus_root) + (1|subject_id) + 
		entrance_type_stagCollapsed*entrance_number*musically_trained + experiment_block,
		data=tapByEntrance_loops, REML=FALSE, na.action=na.omit);
mcSMR.no3wayInt = update(mcSMR.full, .~. -entrance_type_stagCollapsed:entrance_number:musically_trained)
print(anova(mcSMR.no3wayInt,mcSMR.full)) # not significant

# test significance of 2-way interactions involving musical training
cat("\n\nEntrance type X musically trained\n")
mcSMR.noEntrTypeMusTrain_int = update(mcSMR.no3wayInt, .~. -entrance_type_stagCollapsed:musically_trained)
print(anova(mcSMR.noEntrTypeMusTrain_int,mcSMR.no3wayInt)) # not significant

cat("\n\nEntrance period X musically trained\n")
mcSMR.noEntrNum_musTrain_int = update(mcSMR.no3wayInt, .~. -entrance_number:musically_trained)
print(anova(mcSMR.noEntrNum_musTrain_int,mcSMR.no3wayInt)) # significant (p=.02)

# test significance of entrance type X entrance period interaction
cat("\n\nEntrance type X Entrance period\n") 
mcSMR.noEntrType_EntrNum_int = update(mcSMR.no3wayInt, .~. -entrance_type_stagCollapsed:entrance_number)
print(anova(mcSMR.noEntrType_EntrNum_int,mcSMR.no3wayInt)) # significant (p<.001)

# main effect of musical training
cat("\n\nMusically Trained main effect\n")
mcSMR.noMusTrain_ints = update(mcSMR.no3wayInt, .~. -entrance_type_stagCollapsed:musically_trained 
	-entrance_number:musically_trained)
mcSMR.noMusTrain = update(mcSMR.noMusTrain_ints, .~. -musically_trained)
print(anova(mcSMR.noMusTrain,mcSMR.noMusTrain_ints)) # not significant

cat("\n\nInteraction contrasts\n")
print(lsmeans(mcSMR.no3wayInt,list(
	pairwise~entrance_type_stagCollapsed|entrance_number,
	poly~entrance_number|entrance_type_stagCollapsed,
	pairwise~musically_trained|entrance_number,
	poly~entrance_number|musically_trained
	),adjust="fdr"))
	

############
# look at interaction effect between training (non-trained or trained) X effector (hand or head) X stimulus position (beginning or end) on synchronization (SMR)

# only look at baseline and 4th entrance periods
earlyLate_data = subset(MCbyEntranceData, (entrance_number %in% c(1,4)))

# double the lenght of the data frame so we can look can compare SMR for tapping to SMR for head movements
smr_training_data = rbind(earlyLate_data,earlyLate_data)

# create a variable, SMR, containing values of both tapping and head movement SMR
smr_training_data$SMR = c(earlyLate_data$postEntrance_tapping_entrainmentRatio,
	earlyLate_data$postEntrance_mean_entrainmentRatio)

# create categorical factor specifying which movement type the SMR value is associated with
smr_training_data$movement_type = c(
	rep("tapping",length(earlyLate_data$postEntrance_tapping_entrainmentRatio)),
	rep("head_movement",length(earlyLate_data$postEntrance_mean_entrainmentRatio))
	)
smr_training_data$movement_type = factor(smr_training_data$movement_type)
	
cat("\n\n3-way interaction (musically trained X movement type X stim position)\n")
SMR.full = lmer(SMR ~ (1|stimulus_root) + (1|subject_id) + 
		musically_trained*movement_type*entrance_number + experiment_block,
		data=smr_training_data, REML=FALSE, na.action=na.omit);
SMR.no3wayInt = update(SMR.full, .~. -musically_trained:movement_type:entrance_number)
print(anova(SMR.no3wayInt,SMR.full)) # not significant

# 2way interactions
cat("\n\nMusically trained X movement type\n")
SMR.noMusTrain_moveType_int = update(SMR.no3wayInt, .~. -musically_trained:movement_type)
print(anova(SMR.noMusTrain_moveType_int,SMR.no3wayInt)) # significant

cat("\n\nMusically trained X stim position\n")
SMR.noMusTrain_stimPosition_int = update(SMR.no3wayInt, .~. -musically_trained:entrance_number)
print(anova(SMR.noMusTrain_stimPosition_int,SMR.no3wayInt))

cat("\n\nMovement type X stim position\n")
SMR.noMoveType_stimPosition_int = update(SMR.no3wayInt, .~. -movement_type:entrance_number)
print(anova(SMR.noMoveType_stimPosition_int,SMR.no3wayInt))

cat("\n\nContrasts\n")
print(lsmeans(SMR.no3wayInt,list(
	pairwise~movement_type|musically_trained,
	pairwise~musically_trained|movement_type,
	pairwise~musically_trained|entrance_number,
	pairwise~entrance_number|musically_trained
	),adjust="fdr"))
