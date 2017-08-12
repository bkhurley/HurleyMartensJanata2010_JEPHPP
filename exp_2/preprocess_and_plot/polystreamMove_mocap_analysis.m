% analysis list for polystream_move motion capture
%
% 9 Mar 2012 - BH


LOAD_EXISTING = 1;  % Load existing analyses

fstub = 'polystreamMove_mocap_analysis.mat';
globals = polystreamMove_globals;
matfname = fullfile(globals.paths.mocapMatFiles,fstub);

if LOAD_EXISTING
    load(matfname);
end

% Load polystreamMove_preproc_zebris analysis stack
polystreamMove_preproc_zebris;



% get experiment info from Ensemble
na=na+1;
an{na}.name = 'ensemble_load_expinfo';
an{na}.fun  = @ensemble_load_expinfo;
an{na}.params = globals;
an{na}.params.filt.exclude.any.subject_id = [globals.filt.exclude.any.subject_id {'04heg89151','09kir90191','08cem92181'}]; % combining subjects with mocap issues with general to-exclude subjects
an{na}.params.ensemble.experiment_title = globals.experiment_name;
an{na}.params.ensemble.remove_incomplete_sessions = 1;


% remove practice stim data from expinfo
na=na+1;
an{na}.name = 'expinfo_minus_practiceStims';
an{na}.fun = @ensemble_filt_practicedata;
an{na}.params = globals;
an{na}.params.practice_attrib = 'polystream_training';
an{na}.requires = {struct('name','ensemble_load_expinfo')};


% Process rhythm profile for mocap data
na = na + 1;
an{na}.name = 'proc_movement_rp';
an{na}.fun = @polystreamMove_proc_movement_rp;
an{na}.params = globals;
an{na}.params.paths.mocapData = globals.paths.mocapMatFiles;
an{na}.params.paths.metricSurfaces = globals.paths.mocapRP;
an{na}.params.ipem = globals.ipem;
an{na}.params.calcRMS = 0;
an{na}.requires = {struct('name','motion_highpass_filt'),...
    struct('name','expinfo_minus_practiceStims','vars','stimulus_metadata'),...
    struct('name','expinfo_minus_practiceStims','vars','stimulus_x_attribute_metadata')};


% Calculate various metrics
na = na + 1;
an{na}.name = 'analyze_mocap_data';
an{na}.fun = @polystreamMove_analyzeMocapData_v2;
an{na}.params = globals;
an{na}.params.calc_metricNonmetric_resEnergy_ratio = 1;
an{na}.params.phase_analysis.refType = 'stim';
an{na}.params.phase_analysis.targType = 'mocap';
an{na}.params.phase_analysis.plotStimBand = 5;
an{na}.params.phase_analysis.analyze_StimBand = 'max'; % param values can be either an integer corresponding to a 
                                                       % stimulus RP band or the string 'max' corresponding to the 
                                                       % stim band with the highest mean coherence value
an{na}.params.phase_analysis.analyze_axis = 'z';
an{na}.params.phase_analysis.plotMocapAxis = 'z';
an{na}.params.phase_analysis.degreesPerBin = 5;
an{na}.params.phase_analysis.plot_stim_bands = 0;
an{na}.params.phase_analysis.magnitude_plots = 0;
an{na}.params.phase_analysis.plot_stimMPP = 0;
an{na}.params.phase_analysis.relative_phase_plot = 0;
an{na}.params.phase_analysis.rose_plot = 0;
an{na}.params.phase_analysis.coherenceByBand_histogram = 0;
an{na}.params.phase_analysis.hist_figfname = fullfile(globals.paths.movement_figpath,'phaseAnalysis_stimBand_hist.ps');
an{na}.params.phase_analysis.stimBandsPlot_figfname = fullfile(globals.paths.movement_figpath,'phaseAnalysis_stimBandComplex.ps');
an{na}.params.report.fname = fullfile(globals.paths.tablepath,'motion_analysis','tapping_mocap_dataTbl.txt');
an{na}.params.report.write2file = 1;
an{na}.params.report.replace_nan = true;
an{na}.params.report_entrance_analysis.fname = fullfile(globals.paths.tablepath,'motion_analysis','postEntrance_tapping_mocap_dataTbl.txt');
an{na}.params.report_entrance_analysis.write2file = 1;
an{na}.params.report_entrance_analysis.replace_nan = true;
an{na}.params.selectivity_postEntr_table.fname = fullfile(globals.paths.tablepath,'selectivity_postEntr_stim_means.txt');
an{na}.params.selectivity_postEntr_table.write2file = 1;
an{na}.params.selectivity_postEntr_table.replace_nan = 1;
an{na}.params.selectivity_postEntrMeans_table.fname = fullfile(globals.paths.tablepath,'selectivity_postEntr_enterCond_means.txt');
an{na}.params.selectivity_postEntrMeans_table.write2file = 1;
an{na}.params.selectivity_postEntrMeans_table.replace_nan = 1;
an{na}.requires = {struct('name','proc_movement_rp'),...
    struct('name','expinfo_minus_practiceStims','vars','stimulus_metadata')};

%
% Plot Data
%
na = na + 1;
an{na}.name = 'movement_plots';
an{na}.fun = @polystreamMove_movement_plots;
an{na}.params = globals;
an{na}.params.grvLvlCond = 0; % whether or not to divide stims into groove level condition rather than genre (only applies to plot_mean_metricRatio)
an{na}.params.analyzeSelectivity = 1;
an{na}.params.collapse_stagConds = 1;
an{na}.params.plot_metricRatio_tseries = 0;
an{na}.params.plotMocapRP = 0;
an{na}.params.plot_entrainmentRatioFigs = 0;
an{na}.params.plot_mean_tseries.resEnergy = 0;
an{na}.params.plot_mean_tseries.metricEnergy = 0;
an{na}.params.plot_mean_tseries.ratio = 0;
an{na}.params.plot_meanAPS = 0;
an{na}.params.plot_bar_postEntrance_means.resEnergy = 0;
an{na}.params.plot_bar_postEntrance_means.metricEnergy = 0;
an{na}.params.plot_bar_postEntrance_means.metricRatio = 1;
an{na}.params.plot_bar_postEntrance_means.postEntrMPPcorr = 0;
an{na}.params.plot_bar_postEntrance_means.postEntrancePhaseCoh = 0;
an{na}.params.plot_bar_postEntrance_means.barMeans = 0;
an{na}.params.plot_bar_postEntrance_means.tapMocap_corr = 0;
an{na}.params.plot_subject_variability = 0;
an{na}.params.filt.include.all.question_id = [493];
an{na}.params.qid_labels = {493.01,'music_grooved'};
an{na}.params.figfname_metricRatio_tseries = fullfile(globals.paths.movement_figpath,'metricRatio_tseries.ps');
an{na}.params.figfname_meanMetricRatio_tseries = fullfile(globals.paths.movement_figpath,'meanMetricRatio_tseries');
an{na}.params.figfname_meanResEnergy_tseries = fullfile(globals.paths.movement_figpath,'meanResEnergy_tseries.ps');
an{na}.params.figfname_meanAPS = fullfile(globals.paths.movement_figpath,'meanAPS.ps');
an{na}.params.figfname_mocap_tapping_scatter = fullfile(globals.paths.movement_figpath,'mocap_tapping_scatter');
an{na}.params.subRPplot_subid = '03ler90131';
an{na}.params.subRPplot_stimname = 'dance4_strong.mp3';
an{na}.params.plot.savePath = fullfile(globals.paths.movement_figpath,'posterfigs');
an{na}.params.plot.freqLegend = 0;
an{na}.params.plot.legendParams.minHz = 1;
an{na}.params.plot.legendParams.maxHz = 5;
an{na}.params.plot.legendParams.numKeys = 8;
an{na}.params.plot.colorbar = 1;
an{na}.params.plot.color = 0;
an{na}.params.plot.timeTickRes = 4.0;
an{na}.params.plot.numResonTicks = 15;
an{na}.params.plot.perform = {'resonEnergy','stdEnergy'};
an{na}.params.plot.inputSigType = 'impulse';
an{na}.params.plot.separatePages = 0;
an{na}.params.plot.savePlot = 1;
an{na}.params.plot.appendToPlot = 1;
an{na}.requires = {struct('name','analyze_mocap_data'),...
    struct('name','expinfo_minus_practiceStims','vars','stimulus_metadata'),...
    struct('name','expinfo_minus_practiceStims','vars','response_data')};

% Job indices (the order will likely change):
% 1: convert_zebris2mocap
% 2: mocap_matchWithAudio
% 3: despike_zebris_data
% 4: spline_interpolation
% 5: manual_adjustment
% 6: make_motion_movie_v2
% 7: motion_decompose_and_scale
% 8: motion_highpass_filt
% 9: ensemble_load_expinfo
% 10: expinfo_minus_practice_stims
% 11: polystreamMove_proc_movement_rp
% 12: analyze_mocap_data
% 13: movement_plots

params.run_analyses = [12];
an = ensemble_jobman(an,params);

save(matfname,'an');