# wrapper script for all scripts that analyze data from polystream_move experiment
# executing this script from the R command line will run each analysis script and print output to the file specified below in sink()
# 09.17.2013 BH
################

# store the current directory
start.dir=getwd()

# load required packages
require(lme4) # for specifying linear mixed-effects models
require(lsmeans) # for carrying out contrasts
require(nnet) # for logistic regression
require(reshape2) # for releveling variables

# load user-defined functions
source("~/svn/private/R/utils/rm_corr_test.R")

# set polystream_groove R directory in svn
setwd('/Users/bkhurley/svn/private/R/groove/polystream_move')

# set output file
sink("/data/polystream_move/analyses/polystream_move_R_stat_results.txt")

# run subjective response analysis
cat("\n\n\n=============================================================================\n")
cat("\nSTATISTICAL ANALYSIS OF POST-STIMULUS SUBJECTIVE RESPONSES . . .\n")
source("polystream_move_poststim_response_analysis_v2.R")

# model tapping rate
cat("\n\n\n=============================================================================\n")
cat("\nSTATISTICAL ANALYSIS OF TAPPING RATE AS A FUNCTION OF ENTRANCE TYPE . . .\n")
source("lmer_polystream_move_timevar_tapping_analyses_stagCollapsed_v2.R")

# model tapping entrainment
cat("\n\n\n=============================================================================\n")
cat("\nSTATSTICAL ANALYSIS OF TAPPING ENTRAINMENT AS A FUNCTION OF ENTRANCE TYPE . . .\n")
source("lmer_polystream_move_timevar_tapentrain_analyses_stagCollapsed_v2.R")

# model mocap energy
cat("\n\n\n=============================================================================\n")
cat("\nSTATSTICAL ANALYSIS OF MOCAP ENERGY AS A FUNCTION OF ENTRANCE TYPE . . .\n")
source("polystream_move_mocap_analyses_v3.R")

# model tapping entrainment
cat("\n\n\n=============================================================================\n")
cat("\nSTATSTICAL ANALYSIS OF MOCAP ENTRAINMENT AS A FUNCTION OF ENTRANCE TYPE . . .\n")
source("polystream_move_mocapEntrain_analysis_v2.R")

# assess effects related to musicianship
cat("\n\n\n=============================================================================\n")
cat("\nSTATSTICAL ANALYSIS OF MOVEMENT MEASURES AS A FUNCTION OF MUSICIANSHIP . . .\n")
source("polystream_move_musTraining_analysis.R")

# close output file
sink()

# change back to original directory
setwd(initial.dir)