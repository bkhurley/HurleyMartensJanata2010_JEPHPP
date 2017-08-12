# HurleyMartensJanata2014_JEPHPP
Code in this repo uses MATLAB to preprocess, clean, and plot and uses R to statistically model the data from my study published in the Journal of Experimental Psychology: Human Perception and Performance with Peter Martens and Petr Janata.

Read the paper [here](https://bkhurley.github.io/assets/HurleyMartensJanata_2014_JEPHPP.pdf).

Note that the MATLAB code is presented for display purposes only and is not reproducable as is. This is because (1) some data are queried from the lab's MySQL server which requires authentication, and (2) some parts of my MATLAB code rely on lab-general functions that live in a private repo.

## Requirements
The R code uses the following packages:
- `lme4`
- `lsmeans`
- `nnet`
- `reshape2`

## Usage
MATLAB and R code wrapped into the following scripts:
- `polystream_groove_analysis_wrapper.R` 
    * Calls functions for statistically modeling Experiment 1 data
- `polystreamMove_analysis_v1.m` 
    * Calls functions for preprocessing, cleaning, plotting tapping and questionnaire data from Experiment 2
- `polystreamMove_mocap_analysis.m` 
    * Calls functions for preprocessing, cleaning, plotting motion capture data from Experiment 2
- `polystream_move_analysis_wrapper.R` 
    * Calls functions for statistically modeling Experiment 2 data
