% analysis job stack for tapping portion of polystream move experiment
%
% Feb 15, 2011  BH

globals = polystreamMove_globals;

%store form defs in a more easily accessed variable
form_defs = globals.ensemble.form_defs;

LOAD_EXISTING = 1;  % Load existing mat-files to circumvent rerunning analyses

fstub = 'polystreamMove_analysis_v1.mat';
matfname = fullfile(globals.paths.matpath,fstub);

if LOAD_EXISTING
    load(matfname);
end

% Create an array of analysis structures
%clear an;
na = 0;



% get experiment info from Ensemble
na=na+1; %1
an{na}.name = 'expinfo';
an{na}.fun  = @ensemble_load_expinfo;
an{na}.params = globals;
an{na}.params.ensemble.experiment_title = globals.experiment_name;
an{na}.params.ensemble.remove_incomplete_sessions = 1;
% NOTE: subject_id filter should only be used for debugging, and should be commented out when not debugging
%an{na}.params.filt.include.all.subject_id = {'02sij87191'};



% parse midi response
na=na+1; %2
an{na}.name = 'preproc_midi';
an{na}.fun  = @polystreamMove_preproc_midi;
an{na}.params = globals;
an{na}.params.use_saved_preproc = 0; % set to 1 if we wish to use .mat file of previous parsing results
an{na}.params.save_parsed_data = 1; 
an{na}.requires = {struct('name','expinfo')};



% Create a data structure that has the responses to post-stim questions,
% as well as the associated groove responses and put all of the data into 
% a data structure that is oriented by question ID 
%
na = na+1; %3
an{na}.name = 'poststim_responses';
an{na}.fun = @attach_polystreamMove_responses;
an{na}.params.ensemble = globals.ensemble;
an{na}.params.mysql = globals.mysql;
an{na}.params.filt.include.all.question_id = [493 1048 1049];
an{na}.params.minNumTrials = 1;
an{na}.params.matpath = globals.paths.matpath;
an{na}.params.qid_labels = {
    493.01,'music_grooved';
    1048.01,'enjoyed_music';
    1049.01,'like_music_to_continue';
    };
an{na}.requires = {
    struct('name','expinfo','vars','response_data')
    struct('name','preproc_midi')
    };



%
% Analyze subject differences and append subject metrics to ongoing response struct (carried
% over from 'poststim_responses'
%

na = na+1; %4
an{na}.name = 'response_w_sub_analysis';
an{na}.fun = @polystreamMove_analyze_sMESS_tapping;
an{na}.requires = {
    struct('name','expinfo','vars','response_data')
    struct('name','poststim_responses')
    };
an{na}.params = globals.ensemble;
an{na}.params.mysql = globals.mysql;
an{na}.params.plot_subjectselectivity = 1;
an{na}.params.plot_training_selectivity_corr = 1;
an{na}.params.plot_training_MESS_corr = 1;
an{na}.params.plot_tappingMESS_corr = 1;
an{na}.params.plot_MESS_selectivity_corr = 1;
an{na}.params.plot_heardGrv_selectGrpMembership = 1;
an{na}.params.selectivity_colormap = 'default';
an{na}.params.filt.include.all.question_id = [921 922 923 924 925];
an{na}.params.minNumTrials = 1;
an{na}.params.minTapCrit = .85; % cutoff for proportion of trials with less than 10 taps
an{na}.params.matpath = globals.paths.matpath;
an{na}.params.sub_musExp.ignore_response = {'N/A'};
an{na}.params.qid_labels = {
    921.01,'move_with_music_playing';
    922.01,'urge_to_move_with_favorite_music';
    923.01,'move_with_favorite_music_regardless_situation';
    924.01,'match_moves_with_rhythm';
    925.01,'move_to_music_playing_in_head';
    };
an{na}.params.report.write2file = 1;
an{na}.params.report.fname = fullfile(globals.paths.tablepath,'sMESS_tapping_stats.txt');
an{na}.params.report_selectivity.write2file = 0;
an{na}.params.report_selectivity.fname = fullfile(globals.paths.tablepath,'subject_selectivity.txt');
an{na}.params.plot_subDist = 0;
an{na}.params.figfname_sMESS_analysis = fullfile(globals.paths.fig_path,'polystreamMove_subject_analysis.ps');
an{na}.params.figname_subDist = fullfile(globals.paths.fig_path,'subject_distribs.ps');



%
% tapping rhythm profile
%

na=na+1; %5
an{na}.name = 'midiResp_rhythm_profiler';
an{na}.fun = @polystreamMove_midiResp_rhythm_profiler;
an{na}.params.process = {'calcRP'}; %'calcRP','loadRP','calcMeans','recalcTappingStats'
%an{na}.params.filt.include.all.subject_id = {};
an{na}.params.savePath = globals.paths.resp_rp_path;
an{na}.params.fname = 'tappingRP.mat';
an{na}.params.saveCalc = 1;
an{na}.params.nmat2aud.Fs = 100;
an{na}.params.nmat2aud.padStop = 0;
% an{na}.params.nmat2aud.ioiThresh = globals.midi.ioiThresh;
an{na}.params.nmat2aud.extractNotes = globals.midi.validMidiNotes;
an{na}.params.nmat2aud.scaleByVelocity = 1;
% an{na}.params.tapping_stats = globals.tapping_stats;
an{na}.params.ipem_proc = globals.respIPEMParams;
an{na}.params.glob.ignore = 'rp.rpPath';
%param for running coherence calc
% an{na}.params.ipem_proc.rp.perform.interpMPP = 1;
% an{na}.params.ipem_proc.rp.perform.calcOnsetInfo = 1;
% an{na}.params.ipem_proc.rp.perform.calcComplexOutput = 1;
% an{na}.params.ipem_proc.rp.perform.calcVonMises = 1;
% an{na}.params.ipem_proc.rp.perform.calcPhaseCoherence = 1;
an{na}.params.storeIpemProc = {'rp'};
an{na}.requires = {struct('name','preproc_midi'),...
    struct('name','expinfo','vars','stimulus_metadata'),...
    struct('name','response_w_sub_analysis')};


%
% analyze tapping RP
%
na = na + 1; %6
an{na}.name = 'analyze_tapping_rp';
an{na}.fun = @polystreamMove_analyze_tapping_rp;
an{na}.params = globals;
an{na}.params.calc_metricNonmetric_resEnergy_ratio = 1;
an{na}.requires = {struct('name','midiResp_rhythm_profiler'),...
    struct('name','expinfo','vars','stimulus_metadata')};


%
% Perform analyses by attribute
%

na = na+1; %7
an{na}.name = 'response_analyses';
an{na}.fun = @polystreamMove_analyze_responses_v2;
an{na}.requires = {
    struct('name','analyze_tapping_rp'), ...
    struct('name','expinfo','vars','stimulus_metadata'), ...
    struct('name','expinfo','vars','stimulus_x_attribute_metadata'), ...
    };
an{na}.params.ensemble = globals.ensemble;
an{na}.params.songatt = globals.songatt;
an{na}.params.matpath = globals.paths.matpath;
an{na}.params.figpath = globals.paths.fig_path;
an{na}.params.analyze_by_genre = 0;
an{na}.params.analyze_by_grvlvl = 0; % set to 1 if wish to analyze using groove level as factor. default is to use genre as factor.
an{na}.params.analyze_response_selectivity = 1;
an{na}.params.report_outData.fname = fullfile(globals.paths.tablepath,'tapping_data_tbl.txt');
an{na}.params.report_outData.write2file = 1;
an{na}.params.report_outData.replace_nan = true;
an{na}.params.report_entrance_analysis.fname = fullfile(globals.paths.tablepath,'postEntrance_tapping_data_tbl.txt');
an{na}.params.report_entrance_analysis.write2file = 1;
an{na}.params.report_entrance_analysis.replace_nan = true;
an{na}.params.report_by_genre.fname = fullfile(globals.paths.tablepath,'entrancetype_by_genre_tbl.txt');
an{na}.params.report_by_genre.write2file = 1;
an{na}.params.report_by_genre.replace_nan = true;
an{na}.params.tap_by_genre_figfname = fullfile(globals.paths.fig_path,'entrancetype_genre.ps');
an{na}.params.report_by_grvlvl.fname = fullfile(globals.paths.tablepath,'entrancetype_by_grvlvl_tbl.txt');
an{na}.params.report_by_grvlvl.write2file = 1;
an{na}.params.report_by_grvlvl.replace_nan = true;
an{na}.params.tap_by_grvlvl_figfname = fullfile(globals.paths.fig_path,'entrancetype_grvlvl.ps');
an{na}.params.selectivity_analysis_root = globals.paths.tablepath;
an{na}.params.report_selectivity_analysis.write2file = 1;
an{na}.params.report_selectivity_analysis.replace_nan = true;
an{na}.params.selectivity_tapping_figname = fullfile(globals.paths.fig_path,'selectivity_tapping.ps');
an{na}.params.report_grvrank.write2file = 1;
an{na}.params.report_grvrank.fname = fullfile(globals.paths.tablepath,'ranked_grooveratings.txt');
an{na}.params.figname_tapByEntrance = fullfile(globals.paths.fig_path,'tapping_by_entrance.ps');
an{na}.params.report_stimMeans.write2file = 1;
an{na}.params.report_stimMeans.fname = fullfile(globals.paths.tablepath,'stim_taprate_means.txt');
an{na}.params.report_stimMeans.replace_nan = 1;



na = na + 1; %8
an{na}.name = 'poststimTap_corr';
an{na}.fun = @poststim_tapping_analysis;
an{na}.requires = {struct('name','subject_analysis')};
an{na}.params = globals.ensemble;
an{na}.params.mysql = globals.mysql;
an{na}.params.matpath = globals.paths.matpath;
an{na}.params.report.write2file = 1;
an{na}.params.report.fname = fullfile(globals.paths.tablepath,'poststimTap_corr.txt');
an{na}.params.plot_corr = 1;
an{na}.params.figfname_poststimTapCorr = fullfile(globals.paths.fig_path,'poststimTapCorr.ps');



%
% trial stats
%
na = na+1; %9
an{na}.name = 'trial_stats';
an{na}.fun = @polystreamMove_trialStats;

an{na}.requires = {
    struct('name','subject_analysis'),...
    struct('name','expinfo','vars','stimulus_metadata')
    };



%
% plots
%
na = na+1; %10
an{na}.name = 'plots';
an{na}.fun = @polystreamMove_plots_v3;
an{na}.requires = {...
    struct('name','response_analyses')...
    };
an{na}.params = globals;
% specify which plots to run
an{na}.params.tapping_dist = 0;
an{na}.params.tcourse_images = 0;
an{na}.params.collapse_stagConds = 1;
an{na}.params.plot_mean_tapRate = 0;
an{na}.params.plot_entrainmentRatio = 1;
an{na}.params.plot_timecourse = 0;
an{na}.params.plot_poststimCorr = 0;
an{na}.params.plot_subjectStats = 0;
an{na}.params.plot_trialStats = 0;
an{na}.params.plot_poststimTapCorrDists = 0;
an{na}.params.figfname_poststim_means = fullfile(globals.paths.fig_path,'polystreamMove_poststim_barMeans');
an{na}.params.figfname_tappingDist = fullfile(globals.paths.fig_path,'polyMovetapping_dist.ps');
an{na}.params.figfname_mean_tapRate = fullfile(globals.paths.fig_path,'polyMoveBarScatter_tapRate.ps');
an{na}.params.figfname_tapping_poststim_scatter = fullfile(globals.paths.fig_path,'polyMove_tappingPostStim_scatter.ps');
an{na}.params.figfname_tappingEntrainment_poststim_scatter = fullfile(globals.paths.fig_path,'polyMove_tapEntrainment_postStim_scatter.ps');
an{na}.params.figfname_tapRate_poststim_scatter = fullfile(globals.paths.fig_path,'polyMove_tapRate_postStim_scatter.ps');
an{na}.params.figfname_tcourseImg = fullfile(globals.paths.fig_path,'polystreamMove_tcourseImg_figs.ps');
an{na}.params.figfname_tcoursemeans = fullfile(globals.paths.fig_path,'polystream_cumsumMeans_figs.ps');
an{na}.params.figfname_trial_stats = fullfile(globals.paths.fig_path,'stim_tapping_distributions.ps');
an{na}.params.figfname_poststimTapCorrDists = fullfile(globals.paths.fig_path,'poststimTapCorrDists.ps');
an{na}.params.figfname_entrainmentRatio = fullfile(globals.paths.fig_path,'polyMove_tappingEntrainmentRatio.ps');
an{na}.params.figfname_stimcond_fig = fullfile(globals.paths.fig_path,'polystreamMove_stimcond_fig');



jobmanParams.ensemble = globals.ensemble;
jobmanParams.run_analyses = [1];

an = ensemble_jobman(an,jobmanParams);

save(matfname,'an');