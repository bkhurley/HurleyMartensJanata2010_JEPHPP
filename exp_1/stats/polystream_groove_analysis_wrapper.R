# wrapper script for all scripts that analyze data from Experiment 1 in the polystream project.
# executing this script from the R command line will run each analysis script and 
# print output to the file specified below in sink()
#
# 09.17.2013 BH
################

# store the current directory
start.dir<-getwd()

# load required packages
require(lme4) # for specifying linear mixed-effects models
# require(multcomp) 
require(lsmeans) # for carrying out contrasts

# set polystream_groove R directory in svn
setwd('/Users/bkhurley/svn/private/R/groove/polystream_groove')

# load user-defined functions
source("~/svn/private/R/utils/rm_corr_test.R")

# set output file
sink("/data/polystream/groove_rating/analyses/polystream_groove_R_stat_results.txt")

# run entrance type slider-slope analysis
cat("STATISTICAL ANALYSIS OF INTER-ENTRANCE SLIDER-SLOPE RESPONSES AS A FUNCTION OF ENTRANCE TYPE . . .\n")
source("polystream_groove_interentrance_slope_LMM_v3.R")

# run instrument density analysis
cat("\n\n\n=============================================================================\n")
cat("\nSTATISTICAL ANALYSIS OF SLIDER MEANS AS A FUNCTION OF INSTRUMENT DENSITY . . .\n")
source("polystream_groove_instdensity_analysis_LMM_v2.R")

# run subjective response analysis
cat("\n\n\n=============================================================================\n")
cat("\nSTATISTICAL ANALYSIS OF POST-STIMULUS SUBJECTIVE RESPONSES . . .\n")
source("poststim_response_analysis_v2.R")

# close output file
sink()

# change back to original directory
setwd(initial.dir)