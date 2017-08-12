function outData = polystreamMove_analyzeMocapData_v2(inData,params)

% calculate various metrics on mocap data
%
% 11/11/2012 BH
% 03/14/2013 BH - added tapping selectivity level to mocap data structs
% 04/18/2013 BH - recoded to load tapping analysis struct, append mocap
%                 variables to tapping struct, and output a results struct
%                 containing tapping, mocap, and ensemble data (previously
%                 only contained mocap data). This will reduce redundant
%                 coding between analyses and will more easily allow
%                 analyses incorporating all 3 sources of data.
% 05/08/2013 BH - added metric energy variables to analysis output
% 13Jan2013 BH -  added code to handle situations in which multiple MAT
%                 files exist for a stimulus RP or response RP (loads most
%                 recent file)

LOAD_SUB_RP_FILES = 1; % no need to waste time re-loading if still in workspace
CHECK_NONTAPPER_ENTRAINRATIO_VALS = 0;

% load tapping analysis results
tapping_an_fname = fullfile(params.paths.matpath,'tapping_analysis_results_st.mat');
load(tapping_an_fname,'results_st');
tapResCols = set_var_col_const(results_st.vars);
tapResp_subIDs = results_st.data{tapResCols.subject_id};
tapResp_enterCond = results_st.data{tapResCols.entrance_type};
tapResp_enterCondCollapsed = results_st.data{tapResCols.entrance_type_stagCollapsed};
tapResp_stimNames = results_st.data{tapResCols.stim_name};
nsubs = length(unique(tapResp_subIDs));
tapResp_selectivity = results_st.data{tapResCols.selectivity_level};
select_levels = unique(tapResp_selectivity);

clear searchCrit;
searchCrit.name = 'proc_movement_rp';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
mocapData_st = inData{findIdx};
mdCols = set_var_col_const(mocapData_st.vars);
subject_ids = mocapData_st.data{mdCols.subject_id};
unique_subids = unique(subject_ids);
resp_stimNames = mocapData_st.data{mdCols.stim_name};
mocapRP_matfiles = mocapData_st.data{mdCols.mocapRP};

% reorder mocap data to match that of tapping data
mocapRP_matfiles_tapStruct = cell(size(tapResp_subIDs));
for iSub = 1:length(unique_subids)
    curr_sub = unique_subids{iSub};
    mocapRespIdx = ismember(subject_ids,curr_sub);
    tapRespIdx = ismember(tapResp_subIDs,curr_sub);
    mocapRP_matfiles_tapStruct(tapRespIdx) = mocapRP_matfiles(mocapRespIdx);
end
% replace empty cells with ''
emptyRPmat_idx = cellfun(@isempty,mocapRP_matfiles_tapStruct);
mocapRP_matfiles_tapStruct(emptyRPmat_idx) = {''};
nObs = length(mocapRP_matfiles_tapStruct);


clear searchCrit;
searchCrit.type = 'stimulus_metadata';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
stimMetaData = inData{findIdx};
stimMetaDataCols = set_var_col_const(stimMetaData.vars);
stimNames = stimMetaData.data{stimMetaDataCols.name};
stimIDs = stimMetaData.data{stimMetaDataCols.stimulus_id};
nstims = length(stimIDs);

calcRatio = params.calc_metricNonmetric_resEnergy_ratio;

% % load mocap RP structs
% if LOAD_SUB_RP_FILES
%     mocapRP_list = cell(size(mocapRP_matfiles_tapStruct));
%     fprintf('Attempting to load %d mocap Rhythm Profile structures . . .\n',length(mocapRP_matfiles_tapStruct));
%     
%     % loading trial data RP files in parallel broken down by subject, as attempting to load
%     % all files (2k+ files) in one parfor loop generated an error "Attempt to serialize data which is too large"
%     for isub = 1:nsubs
%         curr_start = 1+(isub-1)*41; % 41 trials for each subject
%         curr_end = 41+(isub-1)*41;
%         for iRP = curr_start:curr_end
%             if ~isnan(mocapRP_matfiles_tapStruct{iRP})
%                 mocapRP_list{iRP} = load(mocapRP_matfiles_tapStruct{iRP});
%                 mocapRP_list{iRP} = mocapRP_list{iRP}.rp;
%             else
%                 mocapRP_list{iRP} = '';
%             end
%         end
%     end
% %     mcRPcols = set_var_col_const(mocapRP_list{end}.vars);
% end

% load table of stimulus part entrance times
entrance_table = csvread('/data2/polystream_move/tables/polymove_entrance_times.csv');
entrance_table = sortrows(entrance_table,1);

% Using parfor to iterate over mocap RP structures in parallel 
fprintf('Extracting energy within metric periodicities from %d rhythm profiles . . .\n',nObs);

% preallocate master loop vars
metricNonmetric_resEnergy_ratio = cell(nObs,1);
mc_stim_MPPcorr = cell(nObs,1);
entranceTimes = cell(nObs,1);
meanMetricResEnergyRatio = cell(nObs,1);
postEntrance_meanMetricRatio = cell(nObs,1);
postEntrance_meanMetricEnergy = cell(nObs,1);
postEntrance_PhaseCoh = cell(nObs,1);
resEnergy_tseries = cell(nObs,1);
postEntrance_meanResEnergy = cell(nObs,1);
meanResEnergy = cell(nObs,1);
metric_resEnergy = cell(nObs,1);
meanMetricResEnergy = cell(nObs,1);
stimroot = cell(nObs,1);
metric_reson_mask = cell(nObs,1);
postEntrance_MPPcorr = cell(nObs,1);
postEntrance_MPP_ztransf = cell(nObs,1);

if strcmp(params.phase_analysis.analyze_StimBand,'max')
    maxStimBands = csvread('/data2/polystream_move/tables/polymove_stim_maxCohBands.csv');
    maxStimBands = sortrows(maxStimBands,1);
end

% generate stimname rootname pairing for classification of stimulus root
stim_root = cell(nstims,2);
stimroot_names = fieldnames(params.songatt.exc);
icounter = 0;
for iroot = 1:length(stimroot_names)
    curr_stimroot = stimroot_names{iroot};
    if strcmp('chameleon',curr_stimroot)
        icounter = icounter+1;
        stim_root{icounter,1} = curr_stimroot;
        stim_root(icounter,2) = params.songatt.exc.(curr_stimroot);
    else
        for icond = 1:3
            icounter = icounter+1;
            stim_root{icounter,1} = curr_stimroot;
            stim_root{icounter,2} = params.songatt.exc.(curr_stimroot){icond};
        end
    end
end


% Initiate parallel compute pool
if exist('matlabpool') && ~matlabpool('size') 
    matlabpool
end

stimpeak_counter = 0;
for istim = 1:nstims
    
    % specifiy simultaneous version of each stim, so we can use
    % simultaneous peaks for assessing entrainment ratio
    stimpeak_counter = stimpeak_counter+1;
    curr_stim = strrep(stimNames{istim},'.mp3','');
    curr_peakStim_simCond = strrep(stimNames{istim-(stimpeak_counter-1)},'.mp3','');
%     curr_stimName = strrep(stimNames{istim},'.mp3',''); % remove .mp3 from stim names
    curr_entranceVect = entrance_table(istim,2:end);

%     curr_entranceVect(1) = curr_entranceVect(2)-389;
        
    % remove 0s from entrance column for stims with only 4 entrances
    if ~curr_entranceVect(end)
        curr_entranceVect(end) = []; 
    end
    
    if strcmp(params.phase_analysis.analyze_StimBand,'max')
        curr_StimBand = maxStimBands(istim,2);
    else
        curr_StimBand = params.phase_analysis.analyze_StimBand;
    end
    
    % load stim RP data
    if strcmp(curr_stim(end-2:end),'sim')
        stimRP_simCond_fpath = fullfile('/data/stimuli/audio/groove/polystream_move',...
            curr_peakStim_simCond,'rp');
        stimRP_simCond_file = dir(fullfile(stimRP_simCond_fpath,'*.mat'));
        % if multiple RP files, load most recently modified file
        if length(stimRP_simCond_file) > 1
            dates = datenum(cat(1,stimRP_simCond_file.date));
            [maxdate, date_idx] = max(dates);
            fprintf('\nMultiple RP mat files found for current stimulus. Loading most recent mat file: %s\n',...
                stimRP_simCond_file(date_idx).name)
            stimRP_simCond_file = stimRP_simCond_file(date_idx);
        end
        stimRP_simCond = load(fullfile(stimRP_simCond_fpath,stimRP_simCond_file.name));
        stimRP_simCond = stimRP_simCond.rp;
        stimRP = stimRP_simCond;
        
    else
        stimRP_fpath = fullfile('/data/stimuli/audio/groove/polystream_move',...
            curr_stim,'rp');
        stimRP_file = dir(fullfile(stimRP_fpath,'*.mat'));
        if length(stimRP_file) > 1
            dates = datenum(cat(1,stimRP_file.date));
            [maxdate, date_idx] = max(dates);
            fprintf('\nMultiple RP mat files found for current stimulus. Loading most recent mat file: %s\n',...
                stimRP_file(date_idx).name)
            stimRP_file = stimRP_file(date_idx);
        end
        stimRP = load(fullfile(stimRP_fpath,stimRP_file.name));
        stimRP = stimRP.rp;
    end
    
    % get some stim RP metrics
    RP_cols = set_var_col_const(stimRP.vars);
    stim_resEnergy = stimRP.data{RP_cols.resonatorEnergy};
    stim_MPP = stimRP.data{RP_cols.meanResonatorEnergy};
    
    % index mocap RPs corresponding to current stim
    mocapStim_msk = ismember(tapResp_stimNames,stimNames{istim});
    nStimTrials = sum(mocapStim_msk);
    currMocapRP_fpaths =  mocapRP_matfiles_tapStruct(mocapStim_msk);
    currMocapRPs = cell(nStimTrials,1); 

    % sanity check to make sure all indexed data are associated with expected stim
    curr_stimNameList = tapResp_stimNames(mocapStim_msk);
    curr_stimName = unique(curr_stimNameList);
    
    % find current stim root
    stimroot_mask = ismember(stim_root(:,2),curr_stimName);
    this_stimroot = stim_root{stimroot_mask,1};
    
    % calculate various metrics:
    % - isolate reson energy of mocap RPs to the freq bands related to current stim's metric periodicities
    % - calculate cross correlation betwen mocap MPP & stim MPP
    
    % preallocate subloop vars
    this_mcStimCorr = cell(nStimTrials,1);
    thisEntranceTime = cell(nStimTrials,1);
    resEnergy_ratio = cell(nStimTrials,1);
    metricReson_mask = cell(nStimTrials,1);
    this_meanMetricRatio = cell(nStimTrials,1);
    this_meanMetricEnergy = cell(nStimTrials,1);
    this_postEntranceMeanMetricRatio = cell(nStimTrials,1);
    this_postEntranceMeanMetricEnergy = cell(nStimTrials,1);
    currResEnergy_tseries = cell(nStimTrials,1);
    this_meanResEnergy = cell(nStimTrials,1);
    this_postEntranceMeanResEnergy = cell(nStimTrials,1);
    metricEnergy = cell(nStimTrials,1);
    meanMetricEnergy_tseries = cell(nStimTrials,1);
    this_postEntrance_MPPcorr = cell(nStimTrials,1);
    this_postEntrance_MPP_ztransf = cell(nStimTrials,1);
    
    fprintf('\nCalculating mocap analysis metrics for trials associated with the following stimulus: %s\n',cell2str(curr_stimName))
    
    parfor imcRP = 1:nStimTrials % parfor
        fprintf('.\n')
        if isempty(currMocapRP_fpaths{imcRP}) || any(isnan(currMocapRP_fpaths{imcRP})) % for some reason, fpath indices corresponding to trials rejected due to excessive artifact were populated with NaN
            currMocapRPs{imcRP} = NaN;
            currResEnergy_tseries{imcRP} = NaN;
            this_meanResEnergy{imcRP} = NaN;
            resEnergy_ratio{imcRP} = NaN;
            metricReson_mask{imcRP} = NaN;
            meanMetricEnergy_tseries{imcRP} = NaN;
            this_meanMetricRatio{imcRP} = NaN;
            this_meanMetricEnergy{imcRP} = NaN;
            this_postEntranceMeanMetricRatio{imcRP} = nan(size(curr_entranceVect));
            this_postEntranceMeanMetricEnergy{imcRP} = nan(size(curr_entranceVect));
            this_postEntranceMeanResEnergy{imcRP} = nan(size(curr_entranceVect));
            this_mcStimCorr{imcRP} = NaN;
            this_postEntrance_MPPcorr{imcRP} = nan(size(curr_entranceVect));
            this_postEntrance_MPP_ztransf{imcRP} = nan(size(curr_entranceVect));
        else
            % load current RP mat-file
            currMocapRPs{imcRP} = load(currMocapRP_fpaths{imcRP});
            currMocapRPs{imcRP} = currMocapRPs{imcRP}.rp;
            currResEnergy = currMocapRPs{imcRP}.data{RP_cols.resonatorEnergy};
            % Resonator Energy time series
            currResEnergy_tseries{imcRP} = mean(currResEnergy);
            this_meanResEnergy{imcRP} = mean(currResEnergy_tseries{imcRP});
            
            % calculate timeseries of metric:nonmetric energy ratio based on
            % energy found within and energy found outside peak (i.e., metric)
            % resons
            [metricEnergy{imcRP},resEnergy_ratio{imcRP},metricReson_mask{imcRP}] = isolate_peak_resons('stimRP',stimRP_simCond,...
                'targetRP',currMocapRPs{imcRP},'calcRatio',calcRatio);
            
            meanMetricEnergy_tseries{imcRP} = mean(metricEnergy{imcRP});
            this_meanMetricEnergy{imcRP} = mean(meanMetricEnergy_tseries{imcRP});
            
            % calculate mean metric:nonmetric energy ratio for each trial
            this_meanMetricRatio{imcRP} = mean(resEnergy_ratio{imcRP});
            
            % calculate mean time series values across 1 sec prior to
            % the 2nd entrance and following each subsequent entrance 
            % (3.9 sec is shortest inter-entrance interval across all stims)
            for ientr = 1:length(curr_entranceVect)
                currEntranceSample = curr_entranceVect(ientr);
                if ientr == 1
                    currWindow = (curr_entranceVect(2)-101):(curr_entranceVect(2)-1);
                else
                    currWindow = curr_entranceVect(ientr):(curr_entranceVect(ientr)+389);
                end
                % calculate post-entrance MPP for stimulus & mocap
                this_postEntrance_respMPP = ...
                    mean(currResEnergy(:,currWindow),2);
                this_postEntrance_stimMPP = ...
                    mean(stim_resEnergy(:,currWindow),2);
                this_postEntrance_MPPcorr{imcRP}(ientr) = xcorr(this_postEntrance_stimMPP,this_postEntrance_respMPP,0,'coeff');
                this_postEntrance_MPP_ztransf{imcRP}(ientr) = (1/2)*log((1+this_postEntrance_MPPcorr{imcRP}(ientr))/(1-this_postEntrance_MPPcorr{imcRP}(ientr)));
                this_postEntranceMeanMetricRatio{imcRP}(ientr) = ...
                    mean(resEnergy_ratio{imcRP}(currWindow));
                this_postEntranceMeanMetricEnergy{imcRP}(ientr) = ...
                    mean(meanMetricEnergy_tseries{imcRP}(currWindow));
                this_postEntranceMeanResEnergy{imcRP}(ientr) = ...
                    mean(currResEnergy_tseries{imcRP}(currWindow));
            end
            % mocap by stim MPP xcorr
            mocap_MPP = currMocapRPs{imcRP}.data{RP_cols.meanResonatorEnergy};
            this_mcStimCorr{imcRP} = xcorr(stim_MPP,mocap_MPP,0,'coeff');
            
                       
        end
        % vector of entrance times (in samples)
        thisEntranceTime{imcRP} = curr_entranceVect;
        
    end % parfor imcRP
    metric_reson_mask(mocapStim_msk) = metricReson_mask;
    mc_stim_MPPcorr(mocapStim_msk) = this_mcStimCorr';
    entranceTimes(mocapStim_msk) = thisEntranceTime';
    metricNonmetric_resEnergy_ratio(mocapStim_msk) = resEnergy_ratio';
    metric_resEnergy(mocapStim_msk) =  meanMetricEnergy_tseries';
    meanMetricResEnergyRatio(mocapStim_msk) = this_meanMetricRatio';
    meanMetricResEnergy(mocapStim_msk) = this_meanMetricEnergy';
    postEntrance_meanMetricRatio(mocapStim_msk) = this_postEntranceMeanMetricRatio;
    postEntrance_meanMetricEnergy(mocapStim_msk) = this_postEntranceMeanMetricEnergy;
    resEnergy_tseries(mocapStim_msk) = currResEnergy_tseries;
    meanResEnergy(mocapStim_msk) = this_meanResEnergy;
    postEntrance_meanResEnergy(mocapStim_msk) = this_postEntranceMeanResEnergy;
    stimroot(mocapStim_msk) = {this_stimroot};
    postEntrance_MPPcorr(mocapStim_msk) = this_postEntrance_MPPcorr;
    postEntrance_MPP_ztransf(mocapStim_msk) = this_postEntrance_MPP_ztransf;
    
    if stimpeak_counter == 3
        stimpeak_counter = 0;
    end
    
    postEntrance_PhaseCoh(mocapStim_msk) = polystreamMove_phase_analysis(...
        'stim_name',curr_stimName,'reference_RP',stimRP,'target_RP',currMocapRPs,...
        'entrance_times',curr_entranceVect,'stim_band',curr_StimBand,'params',params.phase_analysis);
%     postEntrance_meanPhaseCoh(mocapStim_msk) = this_mocapPhaseCoh;
end

if CHECK_NONTAPPER_ENTRAINRATIO_VALS
   log_nontapper_stimmeans
end
    

% create cell matrix for post-entrance metric:nonmetric energy ratio &
% write data to txt file (this is a bit of a hack, but gets the job done
% for now)
counter = 0;
for i_mcrp = 1:nObs
    curr_postEntranceRatio = postEntrance_meanMetricRatio{i_mcrp};
    curr_postEntranceResEnergy = postEntrance_meanResEnergy{i_mcrp};
    curr_postEntranceMetricResEnergy = postEntrance_meanMetricEnergy{i_mcrp};
       
    for i_entr = 1:length(curr_postEntranceRatio)
        counter = counter+1;
        entrance_analysis_subids{counter} = tapResp_subIDs{i_mcrp};
        entrance_analysis_stimnames{counter} = tapResp_stimNames{i_mcrp};
        entrance_analysis_tapSelectivity{counter} = results_st.data{tapResCols.selectivity_level}{i_mcrp};
        entrance_analysis_entrType{counter} = results_st.data{tapResCols.entrance_type}{i_mcrp};
        entrance_analysis_entrTypeStagCollapsed{counter} = results_st.data{tapResCols.entrance_type_stagCollapsed}{i_mcrp};
        entrance_analysis_grvLvl{counter} = results_st.data{tapResCols.groove_level_medSplit}{i_mcrp};
        entrance_analysis_genre{counter} = results_st.data{tapResCols.genre}{i_mcrp};
        entrance_analysis_numEntr{counter} = i_entr;
        entrance_analysis_postEntrRatio{counter} = curr_postEntranceRatio(i_entr);
        entrance_analysis_postEntrMetricEnergy{counter} = curr_postEntranceMetricResEnergy(i_entr);
        entrance_analysis_postEntrEnergy{counter} = curr_postEntranceResEnergy(i_entr);
        entrance_analysis_postEntrTapAmt{counter} = results_st.data{tapResCols.post_entrance_taps}{i_mcrp}(i_entr);
        entrance_analysis_postEntrTapRate{counter} = results_st.data{tapResCols.postEntrance_tappingRate}{i_mcrp}(i_entr);
        entrance_analysis_postEntrTapRatio{counter} = results_st.data{tapResCols.tapping_postEntrance_meanMetricRatio}{i_mcrp}(i_entr);
        entrance_analysis_sMESS_movement{counter} = results_st.data{tapResCols.sMESS_movement_subscore}(i_mcrp);
        entrance_analysis_percGroove{counter} = results_st.data{tapResCols.music_grooved}(i_mcrp);
        entrance_analysis_enjoyedMusic{counter} = results_st.data{tapResCols.enjoyed_music}(i_mcrp);
        entrance_analysis_wantMusToCont{counter} = results_st.data{tapResCols.like_music_to_continue}(i_mcrp);
        entrance_analysis_years_training(counter) = results_st.data{tapResCols.years_training}(i_mcrp);
        entrance_analysis_musTrained{counter} = results_st.data{tapResCols.musically_trained}{i_mcrp};
        entrance_analysis_selfConsc_beginning(counter) = results_st.data{tapResCols.self_conscious_beginning}(i_mcrp);
        entrance_analysis_selfCosnc_nearEnd(counter) = results_st.data{tapResCols.self_conscious_near_end}(i_mcrp);
        entrance_analysis_experimentblock(counter) = results_st.data{tapResCols.experiment_block}(i_mcrp);
        entrance_analysis_stimroot{counter} = stimroot{i_mcrp};
        entrance_analysis_MPPcorr(counter) = postEntrance_MPPcorr{i_mcrp}(i_entr);
        entrance_analysis_MPP_ztransf(counter) = postEntrance_MPP_ztransf{i_mcrp}(i_entr);
        entrance_analysis_phaseCoh(counter) = postEntrance_PhaseCoh{i_mcrp}(i_entr);
        entrance_analysis_evrHrdGrv{counter} = results_st.data{tapResCols.ever_heard_groove}{i_mcrp};
    end
end

out_entrData_mtx{1} = entrance_analysis_subids';
out_entrData_mtx{2} = entrance_analysis_stimnames';
out_entrData_mtx{3} = entrance_analysis_stimroot';
out_entrData_mtx{4} = entrance_analysis_entrType';
out_entrData_mtx{5} = entrance_analysis_entrTypeStagCollapsed';
out_entrData_mtx{6} = entrance_analysis_genre';
out_entrData_mtx{7} = entrance_analysis_grvLvl';
out_entrData_mtx{8} = entrance_analysis_tapSelectivity';
out_entrData_mtx{9} = cell2mat(entrance_analysis_sMESS_movement');
out_entrData_mtx{10} = cell2mat(entrance_analysis_percGroove');
out_entrData_mtx{11} = cell2mat(entrance_analysis_enjoyedMusic');
out_entrData_mtx{12} = cell2mat(entrance_analysis_wantMusToCont');
out_entrData_mtx{13} = cell2mat(entrance_analysis_numEntr');
out_entrData_mtx{14} = cell2mat(entrance_analysis_postEntrTapAmt');
out_entrData_mtx{15} = cell2mat(entrance_analysis_postEntrTapRate');
out_entrData_mtx{16} = cell2mat(entrance_analysis_postEntrTapRatio');
out_entrData_mtx{17} = cell2mat(entrance_analysis_postEntrRatio');
out_entrData_mtx{18} = cell2mat(entrance_analysis_postEntrMetricEnergy');
out_entrData_mtx{19} = cell2mat(entrance_analysis_postEntrEnergy');
out_entrData_mtx{20} = entrance_analysis_phaseCoh';
out_entrData_mtx{21} = entrance_analysis_MPPcorr';
out_entrData_mtx{22} = entrance_analysis_MPP_ztransf';
out_entrData_mtx{23} = entrance_analysis_years_training';
out_entrData_mtx{24} = entrance_analysis_musTrained';
out_entrData_mtx{25} = entrance_analysis_selfConsc_beginning';
out_entrData_mtx{26} = entrance_analysis_selfCosnc_nearEnd';
out_entrData_mtx{27} = entrance_analysis_experimentblock';
out_entrData_mtx{28} = entrance_analysis_evrHrdGrv';

data_table = ensemble_init_data_struct;
data_table.vars = {'data','column_labels','column_formats'};
data_table.data{1} = out_entrData_mtx;
data_table.data{2} = {'subject_id','stimulus_name','stimulus_root','entrance_type',...
    'entrance_type_stagCollapsed','genre','groove_level','tapping_selectivity',...
    'sMESS_movement_subscore','music_grooved','enjoyed_music','like_music_to_continue',...
    'entrance_number','postEntrance_taps','postEntrance_tappingRate',...
    'postEntrance_tapping_entrainmentRatio','postEntrance_mean_entrainmentRatio',...
    'postEntrance_mean_metricEnergy','postEntrance_mean_resEnergy',...
    'postEntrance_mocapPhaseCoh','postEntrance_MPPcorr','postEntrance_MPP_ztransf',...
    'years_training','musically_trained','selfconscious_beginning',...
    'selfconsciousness_nearEnd','experiment_block','ever_heard_groove'};
data_table.data{3} = {'%s','%s','%s','%s','%s','%s','%s','%s','%1.2f','%d','%d',...
    '%d','%d','%1.2f','%1.2f','%1.2f','%1.2f','%1.2f','%1.2f','%1.2f','%1.2f',...
    '%1.2f','%d','%s','%d','%d','%d','%s'};
ensemble_display_table(data_table,params.report_entrance_analysis);


% create cell matrix for data output
out_data_mtx{1} = tapResp_subIDs;
out_data_mtx{2} = tapResp_stimNames;
out_data_mtx{3} = stimroot;
out_data_mtx{4} = results_st.data{tapResCols.entrance_type};
out_data_mtx{5} = results_st.data{tapResCols.entrance_type_stagCollapsed};
out_data_mtx{6} = results_st.data{tapResCols.genre};
out_data_mtx{7} = results_st.data{tapResCols.groove_level_medSplit};
out_data_mtx{8} = results_st.data{tapResCols.selectivity_level};
out_data_mtx{9} = results_st.data{tapResCols.sMESS_movement_subscore};
out_data_mtx{10} = results_st.data{tapResCols.music_grooved};
out_data_mtx{11} = results_st.data{tapResCols.enjoyed_music};
out_data_mtx{12} = results_st.data{tapResCols.like_music_to_continue};
out_data_mtx{13} = results_st.data{tapResCols.total_tapping};
out_data_mtx{14} = results_st.data{tapResCols.tapping_rate};
out_data_mtx{15} = cell2mat(results_st.data{tapResCols.tapping_meanMetricRatio});
out_data_mtx{16} = cell2mat(meanMetricResEnergyRatio);
out_data_mtx{17} = cell2mat(meanMetricResEnergy);
out_data_mtx{18} = cell2mat(meanResEnergy);
out_data_mtx{19} = cell2mat(mc_stim_MPPcorr);
out_data_mtx{20} = results_st.data{tapResCols.years_training};
out_data_mtx{21} = results_st.data{tapResCols.musically_trained};
out_data_mtx{22} = results_st.data{tapResCols.self_conscious_beginning};
out_data_mtx{23} = results_st.data{tapResCols.self_conscious_near_end};
out_data_mtx{24} = results_st.data{tapResCols.experiment_block};
out_data_mtx{25} = results_st.data{tapResCols.ever_heard_groove};



% write data to file for stat analysis 
data_table = ensemble_init_data_struct;
data_table.vars = {'data','column_labels','column_formats'};
data_table.data{1} = out_data_mtx;
data_table.data{2} = {'subject_id','stimulus_name','stimulus_root','entrance_type',...
    'entrance_type_stagCollapsed','genre','groove_level','tapping_selectivity',...
    'sMESS_movement_subscore','music_grooved','enjoyed_music','like_music_to_continue',...
    'total_tapping','tapping_rate','mean_tapping_entrainmentRatio','mean_entrainmentRatio','mean_metric_energy',...
    'mean_resEnergy','mocapStim_MPP_xcorr','years_training','musically_trained'...
    'selfconscious_beginning','selfconscious_nearEnd','experiment_block','ever_heard_groove'};
data_table.data{3} = {'%s','%s','%s','%s','%s','%s','%s','%s','%1.2f','%d','%d','%d',...
    '%1.2f','%1.2f','%1.2f','%1.2f','%1.2f','%1.2f','%1.2f','%d','%s','%d','%d','%d','%s'};
ensemble_display_table(data_table,params.report);

% close pool of parallel workers
matlabpool close

% Output struct

% add data to the input mocap data struct
outData = results_st;
outData.name = 'analyze_mocap_data';
outData.vars = [results_st.vars,{'mocapRP','mocap_stim_MPPcorr','metricNonmetric_resEnergy_ratio',...
    'metric_energy_tseries','meanMetricRatio','mean_metric_energy','resEnergy_tseries',...
    'meanResEnergy','postEntrance_meanMetricRatio','postEntrance_meanMetricEnergy','postEntrance_meanResEnergy',...
    'postEntrance_mocapPhaseCoh','stimulus_root','postEntrance_MPPcorr','postEntrance_MPPcorrz','metric_reson_mask'}];
outData.data{39} = mocapRP_matfiles_tapStruct;
outData.data{40} = mc_stim_MPPcorr;
outData.data{41} = metricNonmetric_resEnergy_ratio;
outData.data{42} = metric_resEnergy;
outData.data{43} = meanMetricResEnergyRatio;
outData.data{44} = meanMetricResEnergy;
outData.data{45} = resEnergy_tseries;
outData.data{46} = meanResEnergy;
outData.data{47} = postEntrance_meanMetricRatio;
outData.data{48} = postEntrance_meanMetricEnergy;
outData.data{49} = postEntrance_meanResEnergy;
outData.data{50} = postEntrance_PhaseCoh;
outData.data{51} = stimroot;
outData.data{52} = postEntrance_MPPcorr;
outData.data{53} = postEntrance_MPP_ztransf;
outData.data{54} = metric_reson_mask;

end