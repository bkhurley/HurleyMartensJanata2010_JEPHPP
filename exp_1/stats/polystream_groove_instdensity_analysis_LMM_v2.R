# statistical analyses of polystream_groove slider data 
# fit model to instrument density data
# 09.17.2013 BH
################
# Read in & prep data

# store current directory
initial.dir = getwd()

# change to directory where data reside
data_dir = '/data/polystream/groove_rating/tables'
setwd(data_dir);

# Read in mean entrainment ratio data & replace missing data with NA (R's notation for missing data points)
data_fname = 'slider_groove_data.txt';
cat("\n\n\nData imported from the following file:",sprintf("%s/%s",data_dir,data_fname),"\n\n\n")
sliderdata = read.table(data_fname,header=T,sep='\t',na.strings=c('.'));

# ignore staggered trials for this set of analyses (empty string under insturment_density in data table)
sliderdata_dens = subset(sliderdata, !(instrument_density %in% ""));

################
# Fit linear mixed-effects models
################

# RANDOM EFFECTS

# full model
mod.full = lmer(late_mean ~ 
		instrument_density*musically_trained + (1|root_excerpt) + (1|subject_id),
		data = sliderdata_dens,REML=F,na.action = na.omit);
summary(mod.full);

# test the comparitive fit of random effects structure with (a) stimulus root ('root_excerpt') & subject id and (b) subject id only
cat("\nTest of fit goodness with excerpt as an additional random intercept beyond subject_id\n")
mod.noStimRoot = update(mod.full,.~. - (1|root_excerpt));
print(anova(mod.noStimRoot,mod.full)); # fits significantly better (p<.001) with excerpt & subject id as random intercepts

# EVALUATE FIXED EFFECTS
cat("\n\n\nEvaluating fixed effects . . .\n")
print(summary(mod.full))

# test significance of density:musially_trained interaction
cat("\nSignificance test of density:musically_trained interaction")
mod.noInt = update(mod.full, .~. - instrument_density:musically_trained)
print(anova(mod.noInt,mod.full)) # interaction not significant

# test significance of music training factor
cat("\nSignificance test of musical training\n")
sliderMod.training = lmer(late_mean ~ musically_trained + (1|root_excerpt) + (1|subject_id),
		data = sliderdata,REML=F,na.action = na.omit)
summary(sliderMod.training);
sliderMod.noTraining = update(sliderMod.training, .~. - musically_trained)
print(anova(sliderMod.noTraining,sliderMod.training)) # training factor significant (p<.05)

nonTrained_mean = mean(sliderdata$late_mean[sliderdata$musically_trained=='non-trained'])
nonTrained_SE = sd(sliderdata$late_mean[sliderdata$musically_trained=='non-trained'])/sqrt(length(sliderdata$late_mean[sliderdata$musically_trained=='non-trained']))
trained_mean = mean(sliderdata$late_mean[sliderdata$musically_trained=='trained'])
trained_SE = sd(sliderdata$late_mean[sliderdata$musically_trained=='trained'])/sqrt(length(sliderdata$late_mean[sliderdata$musically_trained=='trained']))
cat("\n\nMean & SE groove rating for non-trained and trained participants\n")
cat(sprintf("Non-trained: M = %1.2f, SE = %1.2f\n",nonTrained_mean,nonTrained_SE))
cat(sprintf("Trained: M = %1.2f, SE = %1.2f\n",trained_mean,trained_SE))

cat("\n\nSubject Descriptor Stats\n")
subjTrain_data = unique(sliderdata[c("subject_id","musically_trained","years_music_experience")])
print(subjTrain_data)
trained_subs = subjTrain_data[subjTrain_data$musically_trained=="trained",]
untrained_subs = subjTrain_data[subjTrain_data$musically_trained=="non-trained",]
cat(sprintf("\n%d out of %d subjects were musically trained.",
	length(trained_subs[,1]),length(subjTrain_data[,1])))
cat(sprintf("\nTrained subs yrs trained: M = %1.2f, SD = %1.2f\n",
	mean(trained_subs[,"years_music_experience"]),sd(trained_subs[,"years_music_experience"])))
cat(sprintf("\nUntrained subs yrs trained: M = %1.2f, SD = %1.2f\n",
	mean(untrained_subs[,"years_music_experience"]),sd(untrained_subs[,"years_music_experience"])))

# assess significance of density factor
cat("\nSignificance test of instrument_density factor\n")
mod.noDens = update(mod.noInt, .~. - instrument_density);
print(anova(mod.noDens,mod.noInt)); # instrument density is significant (p<.0001)

# MULTIPLE COMPARISONS
cat("\n\n\nEvaluating contrasts . . .\n")
# test comparisons among all density levels
cat("\nMultiple comparisons among all levels of instrument density\n\n")
print(lsmeans(mod.noInt,pairwise~instrument_density,adjust="fdr"))
print(lsmeans(sliderMod.training,pairwise~musically_trained,adjust="fdr"))

# return to original directory
setwd(initial.dir)