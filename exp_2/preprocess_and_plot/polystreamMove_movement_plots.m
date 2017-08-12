function outData = polystreamMove_movement_plots(inData,params)

% generates plots for polystream movement data
%
% Nov 19, 2012 BH

clear searchCrit;
searchCrit.name = 'analyze_mocap_data';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
mocapData_st = inData{findIdx};
mdCols = set_var_col_const(mocapData_st.vars);

% remove rows of missing mocap data
missing_mocapRP_idx = cellfun(@(x) all(isnan(x)),mocapData_st.data{mdCols.mocapRP});
nMocapDataVars = length(mocapData_st.data);
for iMCvar = 1:nMocapDataVars
    mocapData_st.data{iMCvar}(missing_mocapRP_idx) = [];
end

mocapRPpaths = mocapData_st.data{mdCols.mocapRP};
mocapRP_metricRatio_tseries = mocapData_st.data{mdCols.metricNonmetric_resEnergy_ratio};
mocapRP_metricEnergy_tseries = mocapData_st.data{mdCols.metric_energy_tseries};
resEnergy_tseries = mocapData_st.data{mdCols.resEnergy_tseries};
mocapSubs = mocapData_st.data{mdCols.subject_id};
subids = unique(mocapSubs);
nsubs = length(subids);
resp_stimNames = mocapData_st.data{mdCols.stim_name};
tapping_selectivity = mocapData_st.data{mdCols.selectivity_level};
groove_level = mocapData_st.data{mdCols.groove_level_medSplit};
nMocapRPs = length(mocapRPpaths);
entranceTimes = mocapData_st.data{mdCols.entrance_times};
% mocap_stim_MPPcorr = mocapData_st.data{mdCols.mocap_stim_MPPcorr};
meanMetricRatio = mocapData_st.data{mdCols.meanMetricRatio};
postEntrance_meanMetricRatio = mocapData_st.data{mdCols.postEntrance_meanMetricRatio};
meanMetricEnergy = mocapData_st.data{mdCols.mean_metric_energy};
postEntrance_meanMetricEnergy = mocapData_st.data{mdCols.postEntrance_meanMetricEnergy};
meanResEnergy = mocapData_st.data{mdCols.meanResEnergy};
postEntrance_meanResEnergy = mocapData_st.data{mdCols.postEntrance_meanResEnergy};
postEntrance_mocapPhaseCoh = mocapData_st.data{mdCols.postEntrance_mocapPhaseCoh};
postEntrance_mppCorrz = mocapData_st.data{mdCols.postEntrance_MPPcorrz};
postEntrance_mppCorr = mocapData_st.data{mdCols.postEntrance_MPPcorr};
data_seq = mocapData_st.data{mdCols.entrance_type};
mppcorr = mocapData_st.data{mdCols.mocap_stim_MPPcorr};
musTrained_groups = unique(mocapData_st.data{mdCols.musically_trained});

clear searchCrit;
searchCrit.type = 'stimulus_metadata';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
stimMetaData = inData{findIdx};
stimMetaDataCols = set_var_col_const(stimMetaData.vars);
stimNames = stimMetaData.data{stimMetaDataCols.name};
stimIDs = stimMetaData.data{stimMetaDataCols.stimulus_id};
nstims = length(stimIDs);
excerptnames = fieldnames(params.songatt.exc);

clear searchCrit;
searchCrit.type = 'response_data';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
respData = inData{findIdx};
repDataCols = set_var_col_const(respData.vars);

% standard plot parameters
label_fontsize = 14;
axes_fontsize = 12;
info_text_fontweight = 'bold';
info_text_fontsize = 14;
title_fontsize = 16;
marker_size = 8;
axis_extra = 0.5;
 
% seq_conditions = fieldnames(params.songatt.seq);
if params.collapse_stagConds
    legend_conds = {'simultaneous','staggered'};
    entrance_type = mocapData_st.data{mdCols.entrance_type_stagCollapsed};
else
    legend_conds = {'simultaneous','staggered strong','staggered weak'};
    entrance_type = mocapData_st.data{mdCols.entrance_type};
end
seq_conditions = unique(entrance_type);
seq_conditions(cellfun(@isempty,seq_conditions)) = [];
% legend_conds = {'simultaneous','staggered strong','staggered weak'};
nconds = size(seq_conditions,1);
genres = fieldnames(params.songatt.genre);
ngenres = length(genres);
groove_levels = unique(groove_level(~cellfun(@isempty,groove_level)));
ngrvLvls = length(groove_levels);

%% plot periodicity surfaces for entrainment ratio methods fig

if params.plot_entrainmentRatioFigs
    params.plot.appendToPlot = 0;
    subid = '11eoj88181';
    stimname = 'dance1_strong';
    subid_mask = ismember(mocapSubs,subid);
    stim_mask = ismember(resp_stimNames,sprintf('%s.mp3',stimname));
    subtrial_mask = subid_mask&stim_mask;
    respRP_fpath = mocapRPpaths{subtrial_mask};
    stimRP_path = fullfile('/data/stimuli/audio/groove/polystream_move',...
        stimname,'rp');
    stimRP_file = dir(fullfile(stimRP_path,'*.mat'));    
    stimRP_fpath = fullfile(stimRP_path,stimRP_file.name);
    metric_reson_mask = mocapData_st.data{mdCols.metric_reson_mask}{subtrial_mask};
    plot_entrainRatio_methodsFigs('respRP',respRP_fpath,'stimRP',stimRP_fpath,...
        'metricResonMask',metric_reson_mask,'params',params);
end

%% plot rhythm profile surfaces and mocap time series

if params.plotMocapRP || params.plot_meanAPS || params.plot_mean_tseries.ratio || ...
        params.plot_mean_tseries.resEnergy
    polystreamMove_plot_mocapAPS_timeseries;
end


%% plot mean energy & entrainment across trial duration & following each entrance

if params.plot_bar_postEntrance_means.metricRatio || params.plot_bar_postEntrance_means.resEnergy ...
        || params.plot_bar_postEntrance_means.metricEnergy || params.plot_bar_postEntrance_means.postEntrMPPcorr ...
        || params.plot_bar_postEntrance_means.postEntrancePhaseCoh
    polystreamMove_plot_mocap_meanData;
end


%% plot subject mocap variability
if params.plot_subject_variability
    polystreamMove_plot_mocap_subj_variability
end

close all
            
outData = [];

return