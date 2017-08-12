function outData = polystreamMove_analyze_tapping_rp(inData,params)

% calculate entrainment ratio of tapping RP/stimulus RP
%
% 13 May 2013 BH wrote script
% 02 Apr 2014 BH added code to select most recent stim RP analysis file if
%                multiple exist

clear searchCrit;
searchCrit.name = 'midiResp_rhythm_profiler';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
tapRP_data_st = inData{findIdx};
rpCols = set_var_col_const(tapRP_data_st.vars);
subject_ids = tapRP_data_st.data{rpCols.subject_id};
unique_subids = unique(subject_ids);
nsubs = length(unique_subids);
tapResp_stimIDs = tapRP_data_st.data{rpCols.stimulus_id};
tapRP_matfiles = tapRP_data_st.data{rpCols.tappingRP};
emptyRPmat_idx = cellfun(@isempty,tapRP_matfiles);
tapRP_matfiles(emptyRPmat_idx) = {''};
nObs = length(tapRP_matfiles);


clear searchCrit;
searchCrit.type = 'stimulus_metadata';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
stimMetaData = inData{findIdx};
stimMetaDataCols = set_var_col_const(stimMetaData.vars);
stimNames = stimMetaData.data{stimMetaDataCols.name}(21:end);
stimIDs = stimMetaData.data{stimMetaDataCols.stimulus_id}(21:end);
nstims = length(stimIDs);



calcRatio = params.calc_metricNonmetric_resEnergy_ratio;

% Initiate parallel compute pool
if exist('matlabpool') && ~matlabpool('size') 
    matlabpool
end

% load mocap RP structs
tapRP_list = cell(size(tapRP_matfiles));
fprintf('Attempting to load %d mocap Rhythm Profile structures . . .\n',length(tapRP_matfiles));

% loading trial data RP files in parallel broken down by subject, as attempting to load 
% all files (2k+ files) in one parfor loop generated an error "Attempt to serialize data which is too large"
for isub = 1:nsubs
    curr_start = 1+(isub-1)*41; % 41 trials for each subject
    curr_end = 41+(isub-1)*41;
    parfor iRP = curr_start:curr_end 
        if ~isnan(tapRP_matfiles{iRP})
            curr_RP = load(tapRP_matfiles{iRP});
            tapRP_list{iRP} = curr_RP.rp;
        else
            tapRP_list{iRP} = '';
        end
    end
end
tapRPcols = set_var_col_const(tapRP_list{end}.vars);

% load table of stimulus part entrance times
entrance_table = csvread('/data2/polystream_move/tables/polymove_entrance_times.csv');
entrance_table = sortrows(entrance_table,1);

% Using parfor to iterate over mocap RP structures in parallel 
fprintf('Extracting energy within metric periodicities from %d rhythm profiles . . .\n',nObs);

% preallocate master loop vars
metricNonmetric_resEnergy_ratio = cell(nObs,1);
tapping_stim_MPPcorr = cell(nObs,1);
meanMetricResEnergyRatio = cell(nObs,1);
postEntrance_meanMetricRatio = cell(nObs,1);
postEntrance_meanMetricEnergy = cell(nObs,1);
resEnergy_tseries = cell(nObs,1);
postEntrance_meanResEnergy = cell(nObs,1);
meanResEnergy = cell(nObs,1);
metric_resEnergy = cell(nObs,1);
meanMetricResEnergy = cell(nObs,1);

stimpeak_counter = 0;
for istim = 1:nstims
    
    % specifiy simultaneous version of each stim, so we can use
    % simultaneous peaks for assessing entrainment ratio
    stimpeak_counter = stimpeak_counter+1;
    curr_peakStim = strrep(stimNames{istim-(stimpeak_counter-1)},'.mp3','');
    curr_entranceVect = entrance_table(istim,2:end);
        
    % remove 0s from entrance column for stims with only 4 entrances
    if ~curr_entranceVect(end)
        curr_entranceVect(end) = []; 
    end
    % load stim RP data
    stimRP_fpath = fullfile('/data/stimuli/audio/groove/polystream_move',...
        curr_peakStim,'rp');
    stimRP_file = dir(fullfile(stimRP_fpath,'*.mat'));
    
    % if multiple stim RP files, select most recent one
    if length(stimRP_file) > 1
        stimRP_date_vect = nan(size(stimRP_file));
        for iStimRP = 1:length(stimRP_file)
            stimRP_date_vect(iStimRP) = stimRP_file(iStimRP).datenum;
        end
        [~,maxDate_idx] = max(stimRP_date_vect);
        stimRP_file = stimRP_file(maxDate_idx);
    end
    stimRP_st = load(fullfile(stimRP_fpath,stimRP_file.name));
    stimRP_st = stimRP_st.rp;
    RP_cols = set_var_col_const(stimRP_st.vars);
    stim_MPP = stimRP_st.data{RP_cols.meanResonatorEnergy};
    
    % index tapping RPs corresponding to current stim
    tapStim_idx = ismember(tapResp_stimIDs,stimIDs(istim));
    currTapRPs = tapRP_list(tapStim_idx)'; 
    % sanity check to make sure all indexed data are associated with expected stim
    curr_stimIDList = tapResp_stimIDs(tapStim_idx);
    curr_stimID = unique(curr_stimIDList);
    
    % calculate various metrics:
    % - isolate reson energy of tapping RPs to the freq bands related to
    %   current stim's metric periodicities, and calculate the ratio of w/in
    %   metric freqs / outside metric freqs (entrainment ratio)
    % - calculate cross correlation betwen tapping MPP & stim MPP
    
    % preallocate subloop vars
    this_tappingStimCorr = cell(size(currTapRPs));
    thisEntranceTime = cell(size(currTapRPs));
    resEnergy_ratio = cell(size(currTapRPs));
    this_meanMetricRatio = cell(size(currTapRPs));
    this_meanMetricEnergy = cell(size(currTapRPs));
    this_postEntranceMeanMetricRatio = cell(size(currTapRPs));
    this_postEntranceMeanMetricEnergy = cell(size(currTapRPs));
    currResEnergy_tseries = cell(size(currTapRPs));
    this_meanResEnergy = cell(size(currTapRPs));
    this_postEntranceMeanResEnergy = cell(size(currTapRPs));
    metricEnergy = cell(size(currTapRPs));
    meanMetricEnergy_tseries = cell(size(currTapRPs));
    
    parfor imcRP = 1:length(currTapRPs) % parfor
        if ~isstruct(currTapRPs{imcRP})
            currResEnergy = NaN;
            currResEnergy_tseries{imcRP} = NaN;
            this_meanResEnergy{imcRP} = NaN;
            resEnergy_ratio{imcRP} = NaN;
            meanMetricEnergy_tseries{imcRP} = NaN;
            this_meanMetricRatio{imcRP} = NaN;
            this_meanMetricEnergy{imcRP} = NaN;
            this_postEntranceMeanMetricRatio{imcRP} = nan(size(curr_entranceVect));
            this_postEntranceMeanMetricEnergy{imcRP} = nan(size(curr_entranceVect));
            this_postEntranceMeanResEnergy{imcRP} = nan(size(curr_entranceVect));
            this_tappingStimCorr{imcRP} = NaN;
        else
            
            currTapRP_st = currTapRPs{imcRP};
            currResEnergy = currTapRP_st.data{RP_cols.resonatorEnergy};
            % Resonator Energy time series
            currResEnergy_tseries{imcRP} = mean(currResEnergy);
            this_meanResEnergy{imcRP} = mean(currResEnergy_tseries{imcRP});
            
            % calculate timeseries of metric:nonmetric energy ratio based on
            % energy found within and energy found outside peak (i.e., metric)
            % resons
            [metricEnergy{imcRP},resEnergy_ratio{imcRP},~] = isolate_peak_resons('stimRP',stimRP_st,...
                'targetRP',currTapRP_st,'calcRatio',calcRatio);
            
            meanMetricEnergy_tseries{imcRP} = mean(metricEnergy{imcRP});
            this_meanMetricEnergy{imcRP} = mean(meanMetricEnergy_tseries{imcRP});
            
            % calculate mean metric:nonmetric energy ratio for each trial
            this_meanMetricRatio{imcRP} = nanmean(resEnergy_ratio{imcRP});
            
            % calculate mean time series values across 3.89 sec prior to
            % the 2nd entrance and following each subsequent entrance 
            % (3.9 sec is shortest inter-entrance interval across all stims)
            for ientr = 1:length(curr_entranceVect)
                currEntranceSample = curr_entranceVect(ientr);
                if ientr == 1
                    currWindow = (curr_entranceVect(2)-101):(curr_entranceVect(2)-1);
                else
                    currWindow = curr_entranceVect(ientr):(curr_entranceVect(ientr)+389);
                end
                this_postEntranceMeanMetricRatio{imcRP}(ientr) = ...
                    mean(resEnergy_ratio{imcRP}(currWindow));
                this_postEntranceMeanMetricEnergy{imcRP}(ientr) = ...
                    mean(meanMetricEnergy_tseries{imcRP}(currWindow));
                this_postEntranceMeanResEnergy{imcRP}(ientr) = ...
                    mean(currResEnergy_tseries{imcRP}(currWindow));
            end
            % mocap by stim MPP xcorr
            tapping_MPP = currTapRP_st.data{RP_cols.meanResonatorEnergy};
            this_tappingStimCorr{imcRP} = xcorr(stim_MPP,tapping_MPP,0,'coeff');
                       
        end
        % vector of entrance times (in samples)
        thisEntranceTime{imcRP} = curr_entranceVect;
        
    end % for imcRP
    
    tapping_stim_MPPcorr(tapStim_idx) = this_tappingStimCorr';
    entranceTimes(tapStim_idx) = thisEntranceTime';
    metricNonmetric_resEnergy_ratio(tapStim_idx) = resEnergy_ratio';
    metric_resEnergy(tapStim_idx) =  meanMetricEnergy_tseries';
    meanMetricResEnergyRatio(tapStim_idx) = this_meanMetricRatio';
    meanMetricResEnergy(tapStim_idx) = this_meanMetricEnergy';
    postEntrance_meanMetricRatio(tapStim_idx) = this_postEntranceMeanMetricRatio;
    postEntrance_meanMetricEnergy(tapStim_idx) = this_postEntranceMeanMetricEnergy;
    resEnergy_tseries(tapStim_idx) = currResEnergy_tseries;
    meanResEnergy(tapStim_idx) = this_meanResEnergy;
    postEntrance_meanResEnergy(tapStim_idx) = this_postEntranceMeanResEnergy;
    
    if stimpeak_counter == 3
        stimpeak_counter = 0;
    end
    
end

% close pool of parallel workers
matlabpool close

%% Output struct

% add RP data to the ongoing data struct
outData = tapRP_data_st;
outData.name = 'analyze_tapping_rp';
outData.vars = [tapRP_data_st.vars,{'tapping_stim_MPPcorr'},{'tapping_entrainment_ratio'},...
    {'tapping_metric_energy_tseries'},{'tapping_meanMetricRatio'},{'tapping_mean_metric_energy'},...
    {'tapping_resEnergy_tseries'},{'tapping_meanResEnergy'},{'tapping_postEntrance_meanMetricRatio'},...
    {'tapping_postEntrance_meanMetricEnergy'},{'tapping_postEntrance_meanResEnergy'}];
outData.data = [tapRP_data_st.data {tapping_stim_MPPcorr} {metricNonmetric_resEnergy_ratio} ...
    {metric_resEnergy} {meanMetricResEnergyRatio} {meanMetricResEnergy} {resEnergy_tseries} ...
    {meanResEnergy} {postEntrance_meanMetricRatio} {postEntrance_meanMetricEnergy} ...
    {postEntrance_meanResEnergy}];

end