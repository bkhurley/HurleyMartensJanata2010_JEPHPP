function results_st = polystreamMove_analyze_responses_v2(data_st,params)

% Organizes subject responses from polystream_move experiment and outputs
% tables for statistical analyses
%
% Modeled after polystream_analyze_responses
%

% 04 Mar 2012 BH - initial script
% 20 Aug 2013 BH - added additionaly variable 'entrance_type_stagCollapsed'
%                  in which staggered-strong & staggered-weak conditions
%                  are collapsed to a common staggered condition to enable
%                  easy switching between collapsed & uncollapsed staggered
%                  conditions in subsequent analysis scripts

results_st = ensemble_init_data_struct;

%
% Get the mean ratings by stimid
%
clear crit
crit.name = 'analyze_tapping_rp';
idx = ensemble_find_analysis_struct(data_st,crit);
resp_st = data_st{idx};
respcols = set_var_col_const(resp_st.vars);
resp_stimIDs = resp_st.data{respcols.stimulus_id};
resp_subIDs = resp_st.data{respcols.subject_id};
stimulus_ids = unique(resp_stimIDs);
nstims = length(stimulus_ids);
total_tapping = resp_st.data{respcols.total_tapping};
tapping_meanEntrainmentRatio = cell2mat(resp_st.data{respcols.tapping_meanMetricRatio});
resamp_midiResp = resp_st.data{respcols.resamp_midiResp};
nObs = length(resp_st.data{respcols.subject_id}); % num observations across subjects & conditions

% Get the stimulus structure
crit.name = 'stimulus_metadata';
idx = ensemble_find_analysis_struct(data_st,crit);
stim_st = data_st{idx};
stimcols = set_var_col_const(stim_st.vars);
stimMeta_stimIDs = stim_st.data{stimcols.stimulus_id};
stim_names = stim_st.data{stimcols.name};

% initiate output struct
results_st = resp_st;

% Create matrix that masks for locations within resp_st.data that
% correspond to unique values in stimids
[stim_mask_mtx, stimids] = make_mask_mtx(resp_st.data{respcols.stimulus_id});
nstims = length(stimids);

% load table of stimulus part entrance times
entrance_table = csvread('/data2/polystream_move/tables/polymove_entrance_times.csv');
entrance_table = sortrows(entrance_table,1);

stim_attrib_vars = {'music_grooved','enjoyed_music','like_music_to_continue','total_tapping',...
    'tapping_meanMetricRatio'};
nvars = length(stim_attrib_vars);
stimatt_cols = set_var_col_const(stim_attrib_vars);

% initialize stim x attribute matrix
stim_stats = {'nanmean','nanstd','nanmin','nanmax'};
for istat = 1:length(stim_stats)
    curr_stat = stim_stats{istat};
    stim_att_mtx.(curr_stat) = zeros(nstims,nvars);
end

for ivar = 1:nvars
    curr_var = stim_attrib_vars{ivar};
    if iscell(resp_st.data{respcols.(curr_var)})
        resp_st.data{respcols.(curr_var)} = cell2mat(resp_st.data{respcols.(curr_var)});
    end
    for istim = 1:nstims
        for istat = 1:length(stim_stats)
            curr_stat = stim_stats{istat};
            fh = str2func(curr_stat);
            stim_att_mtx.(curr_stat)(istim,ivar) = ...
                fh(resp_st.data{respcols.(curr_var)}(stim_mask_mtx(:,istim))); % apply curr_stat to all responses of curr_var of curr_stim
        end % for istat=
    end % for istim=
end % for ivar=

rank_stims_by_perceived_groove

% pair up stimulus-entrance type associations
enterType_stimName{1} = cell(nstims-1,1);
enterType_stimName{2} = cell(nstims-1,1);
entrCond_fields = fieldnames(params.songatt.seq);
for iEntrCond = 1:length(entrCond_fields)
    curr_entrCond = entrCond_fields{iEntrCond};
    start_idx = iEntrCond*20-19;
    stop_idx = iEntrCond*20;
    enterType_stimName{1}(start_idx:stop_idx) = {curr_entrCond};
    enterType_stimName{2}(start_idx:stop_idx) = params.songatt.seq.(curr_entrCond);
end

% pair up stimulus-genre associations
genre_stimName{1} = cell(nstims-1,1);
genre_stimName{2} = cell(nstims-1,1);
genre_fields = fieldnames(params.songatt.genre);
for igenre = 1:length(genre_fields)
    curr_genre = genre_fields{igenre};
    startIdx = igenre*15-14;
    stopIdx = igenre*15;
    genre_stimName{1}(startIdx:stopIdx,1) = {curr_genre};
    genre_stimName{2}(startIdx:stopIdx,1) = params.songatt.genre.(curr_genre);
end

% Initiate parallel compute pool
if exist('matlabpool') && ~matlabpool('size') 
    matlabpool
end

entrance_type = cell(size(resp_stimIDs));
grvlvl_highLow = cell(size(resp_stimIDs));
% grvlvl_continuous = nan(size(resp_stimIDs));
genre = cell(size(resp_stimIDs));
resp_stimNames = cell(size(resp_stimIDs));
postEntranceTotalTapping = cell(size(resp_stimIDs));
tapping_rate = nan(size(resp_stimIDs));
postEntrance_tapping_rate = cell(size(resp_stimIDs));
entrance_times = cell(size(resp_stimIDs));

for istim = 1:nstims
    curr_stimID = stimids(istim);
    % get entrance time info
    curr_entranceVect = entrance_table(istim,2:end);    
    % remove 0s from entrance column for stims with only 4 entrances
    if ~curr_entranceVect(end)
        curr_entranceVect(end) = [];
    end
    
    curr_respIdx = stim_mask_mtx(:,istim);
    curr_stimMetaIdx = ismember(stimMeta_stimIDs,curr_stimID);
    curr_stimName = stim_names{curr_stimMetaIdx};
    resp_stimNames(curr_respIdx) = {curr_stimName};
    
    fprintf('Working on tapping responses associated with stimulus %s\n',curr_stimName);
    
    % address chameleon trials in which entrance typ, genre, and groove
    % level does not apply
    if ismember({curr_stimName},{'chameleon.mp3'})
        entrance_type(curr_respIdx) = {''};
        genre(curr_respIdx) = {''};
        grvlvl_highLow(curr_respIdx) = {''};
    else
        enterType_idx = ismember(enterType_stimName{2},{curr_stimName});
        entrance_type(curr_respIdx) = enterType_stimName{1}(enterType_idx);
        genre_idx = ismember(genre_stimName{2},curr_stimName);
        genre(curr_respIdx) = genre_stimName{1}(genre_idx);
        if sum(ismember(low_groove_stims,curr_stimName))
            grvlvl_highLow(curr_respIdx) = {'low_groove'};
        elseif sum(ismember(high_groove_stims,curr_stimName))
            grvlvl_highLow(curr_respIdx) = {'high_groove'};
        else
            error('no groove level information for %s',curr_stimName)
        end
    end
    
    % find all observations for this stimulus
    curr_totalTapData = total_tapping(curr_respIdx);
    curr_tap_tseries = resamp_midiResp(curr_respIdx);
    this_tappingRate = nan(size(curr_totalTapData));
    this_postEntranceTappingRate = cell(size(curr_totalTapData));
    this_postEntranceTotalTapping = cell(size(curr_totalTapData));
    thisEntranceTime = cell(size(curr_totalTapData));
    
    parfor iObs = 1:length(curr_totalTapData)
        
        % calculate mean time series values across 3.89 sec prior to
        % the 2nd entrance and following each subsequent entrance
        % (3.9 sec is shortest inter-entrance interval across all stims)
        this_postEntranceTotalTapping{iObs} = nan(size(curr_entranceVect));
        for ientr = 1:length(curr_entranceVect)
%             currEntranceSample = curr_entranceVect(ientr);
            if ientr == 1
                currWindow = (curr_entranceVect(2)-101):(curr_entranceVect(2)-1);
            else
                currWindow = curr_entranceVect(ientr):(curr_entranceVect(ientr)+389);
            end
            this_postEntranceTotalTapping{iObs}(ientr) = ...
                sum(curr_tap_tseries{iObs}(currWindow));
            this_postEntranceTappingRate{iObs}(ientr) = ...
                this_postEntranceTotalTapping{iObs}(ientr)/(length(currWindow)/100);
        end
        % vector of entrance times (in samples)
        thisEntranceTime{iObs} = curr_entranceVect;
        this_tappingRate(iObs) = curr_totalTapData(iObs)/(length(curr_tap_tseries{iObs})/100);
    end
    
    postEntranceTotalTapping(curr_respIdx) = this_postEntranceTotalTapping;
    entrance_times(curr_respIdx) = thisEntranceTime;
    tapping_rate(curr_respIdx) = this_tappingRate;
    postEntrance_tapping_rate(curr_respIdx) = this_postEntranceTappingRate;
end

% create a 2nd entrance type variable (entrance_type_stagCollapsed) in which
% staggered-strong and staggered-weak conditions have are collapsed to a
% common staggered ('stag') condition
stag_str_msk = ismember(entrance_type,'stagStrong');
stag_wk_msk = ismember(entrance_type,'stagWeak');
entrance_type_stagCollapsed = entrance_type;
entrance_type_stagCollapsed(stag_str_msk) = {'stag'};
entrance_type_stagCollapsed(stag_wk_msk) = {'stag'};

% create cell matrix for post-entrance responses & write data to txt file 
counter = 0;
for iobs = 1:nObs
    curr_postEntranceTotalTap = postEntranceTotalTapping{iobs};
       
    for i_entr = 1:length(curr_postEntranceTotalTap)
        counter = counter+1;
        entrance_analysis_subids{counter} = resp_subIDs{iobs};
        entrance_analysis_stimNames{counter} = resp_stimNames{iobs};
        entrance_analysis_musicGrooved(counter) = resp_st.data{respcols.music_grooved}(iobs);
        entrance_analysis_enjoyedMusic(counter) = resp_st.data{respcols.enjoyed_music}(iobs);
        entrance_analysis_likeMusContinue(counter) = resp_st.data{respcols.like_music_to_continue}(iobs);
        entrance_analysis_selectivityIndex(counter) = resp_st.data{respcols.selectivity_index}(iobs);
        entrance_analysis_selectivityLvl{counter} = resp_st.data{respcols.selectivity_level}(iobs);
        entrance_analysis_sMESSmvmnt(counter) = resp_st.data{respcols.sMESS_movement_subscore}(iobs);
        entrance_analysis_entrType{counter} = entrance_type{iobs};
        entrance_analysis_grvLvl{counter} = grvlvl_highLow{iobs};
        entrance_analysis_genre{counter} = genre{iobs};
        entrance_analysis_numEntr(counter) = i_entr;
        entrance_analysis_postEnterTap(counter) = postEntranceTotalTapping{iobs}(i_entr);
        entrance_analysis_postEnterTappingRate(counter) = postEntrance_tapping_rate{iobs}(i_entr);
        entrance_analysis_postEnterTappingEntrainment(counter) = resp_st.data{respcols.tapping_postEntrance_meanMetricRatio}{iobs}(i_entr);
        entrance_analysis_years_training(counter) = resp_st.data{respcols.years_training}(iobs);
        entrance_analysis_selfConsc_beginning(counter) = resp_st.data{respcols.self_conscious_beginning}(iobs);
        entrance_analysis_selfCosnc_nearEnd(counter) = resp_st.data{respcols.self_conscious_near_end}(iobs);
        entrance_analysis_experimentblock(counter) = resp_st.data{respcols.experiment_block}(iobs);       
        entrance_analysis_entrTypeStagColl{counter} = entrance_type_stagCollapsed{iobs};
    end
end

out_entrData_mtx = {...
    entrance_analysis_subids' ...
    entrance_analysis_stimNames' ...
    entrance_analysis_entrType' ...
    entrance_analysis_entrTypeStagColl' ...
    entrance_analysis_genre' ...
    entrance_analysis_grvLvl' ...
    entrance_analysis_selectivityIndex' ...
    entrance_analysis_selectivityLvl' ...
    entrance_analysis_numEntr' ...
    entrance_analysis_musicGrooved' ...
    entrance_analysis_enjoyedMusic' ...
    entrance_analysis_likeMusContinue' ...
    entrance_analysis_sMESSmvmnt' ...
    entrance_analysis_postEnterTappingRate' ...
    entrance_analysis_postEnterTap' ...
    entrance_analysis_postEnterTappingEntrainment' ...
    entrance_analysis_years_training' ...
    entrance_analysis_selfConsc_beginning' ...
    entrance_analysis_selfCosnc_nearEnd' ...
    entrance_analysis_experimentblock' ...
    };

data_table = ensemble_init_data_struct;
data_table.vars = {'data','column_labels','column_formats'};
data_table.data{1} = out_entrData_mtx;
data_table.data{2} = {'subject_id','stimulus_name','entrance_type','entrance_type_stagCollapsed',...
    'genre','groove_level','selectivity_index','selectivity_level','entrance_number',...
    'music_grooved','enjoyed_music','like_music_to_continue','sMESS_movement_subscore',...
    'postEntrance_tappingRate','postEntrance_totalTapping','postEntrance_tapEntrainmentRatio',...
    'years_training','selfconscious_beginning','selfconscious_nearEnd','experiment_block'};
data_table.data{3} = {'%s','%s','%s','%s','%s','%s','%1.2f','%s','%d','%d',...
    '%d','%d','%1.2f','%1.2f','%1.2f','%1.2f','%d','%d','%d','%d'};
ensemble_display_table(data_table,params.report_entrance_analysis);

results_st.name = 'response_analyses';
results_st.vars = [results_st.vars {'stim_name'} {'groove_level_medSplit'} {'genre'}...
    {'entrance_type'} {'entrance_type_stagCollapsed'} {'entrance_times'} {'tapping_rate'}...
    {'postEntrance_tappingRate'} {'post_entrance_taps'}];
results_st.data = [results_st.data {resp_stimNames} {grvlvl_highLow} {genre} {entrance_type}...
    {entrance_type_stagCollapsed} {entrance_times} {tapping_rate} {postEntrance_tapping_rate}...
    {postEntranceTotalTapping}];

resCols = set_var_col_const(results_st.vars);

% create cell matrix for data output
out_data_mtx = {...
    results_st.data{resCols.subject_id} ...
    results_st.data{resCols.stim_name} ...
    results_st.data{resCols.entrance_type} ...
    results_st.data{resCols.entrance_type_stagCollapsed} ...
    results_st.data{resCols.genre} ...
    results_st.data{resCols.groove_level_medSplit} ...
    results_st.data{resCols.selectivity_level} ...
    results_st.data{resCols.music_grooved} ...
    results_st.data{resCols.enjoyed_music} ...
    results_st.data{resCols.like_music_to_continue} ...
    results_st.data{resCols.sMESS_movement_subscore} ...
    results_st.data{resCols.total_tapping} ...
    results_st.data{resCols.tapping_rate} ...
    cell2mat(results_st.data{resCols.tapping_meanMetricRatio}) ...
    results_st.data{resCols.years_training} ...
    results_st.data{resCols.self_conscious_beginning} ...
    results_st.data{resCols.self_conscious_near_end} ...
    results_st.data{resCols.experiment_block} ...
};

% write data to file for stat analysis 
data_table = ensemble_init_data_struct;
data_table.vars = {'data','column_labels','column_formats'};
data_table.data{1} = out_data_mtx;
data_table.data{2} = {'subject_id','stimulus_name','entrance_type','entrance_type_stagCollapsed',...
    'genre','groove_level','tapping_selectivity','music_grooved','enjoyed_music','like_music_to_continue',...
    'sMESS_movement_subscore','total_tapping','tapping_rate','mean_tappingEntrainmentRatio',...
    'years_training','selfconscious_beginning','selfconscious_end','experiment_block'};
data_table.data{3} = {'%s','%s','%s','%s','%s','%s','%s','%d','%d','%d','%1.2f','%1.2f','%1.2f','%1.2f',...
    '%d','%d','%d','%d'};
ensemble_display_table(data_table,params.report_outData);

% save results to file so they can be loaded by mocap analysis
matfname = fullfile(params.matpath,'tapping_analysis_results_st.mat');
fprintf('Saving results to file: %s\n',matfname);
save(matfname,'results_st');

% close pool of parallel workers
matlabpool close

return