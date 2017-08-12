function result_st = polystreamMove_analyze_sMESS_tapping(inData, params)

% - Analyzes sMESS movement subscores and tapping distributions for each subjects. 
% - Dynamically removes subjects that produce less than 10 taps on 70% or more trials.
% - Plots results and ouputs revised response struct with omitted subjects
%
% 4/24/2012 Brian Hurley

clear anSearchCrit;
anSearchCrit.name = 'poststim_responses';
anListIdx = ensemble_find_analysis_struct(inData,anSearchCrit);
tapResp_st = inData{anListIdx};
tapRespCols = set_var_col_const(tapResp_st.vars);

clear anSearchCrit;
anSearchCrit.name = 'response_data';
anListIdx = ensemble_find_analysis_struct(inData,anSearchCrit);
resp_st = inData{anListIdx};
resp_cols = set_var_col_const(resp_st.vars);

% midi response vars
resp_vars = {'tapping_mean','tapping_sd','total_tapping','sMESS_score'};
nvars = length(resp_vars);

result_st = ensemble_init_data_struct;
result_st = tapResp_st;
result_st.name = 'response_w_sub_analysis';
resCols = set_var_col_const(result_st.vars);

% Generate a data structure that is oriented by composite question IDs
qid_st = ensemble_data_by_question(resp_st,params);

%% calculate tap stats

% Get a list of subjects
curr_parms.extract_var = 'subject_id';
subid_vect = tapResp_st.data{tapRespCols.subject_id};
subids = unique(subid_vect);
nsub = length(subids);

% initialize tap stat vects
sub_tap_trialVect = cell(nsub,1);
sub_tap_subid = cell(nsub,1);
sub_tap_mean = zeros(nsub,1);
sub_tap_sd = zeros(nsub,1);
sub_tap_total = zeros(nsub,1);
nontapperIdx = 0;
nontapperSubs = {};
nontapperSubVectIdx = [];
resp_select_lvl = cell(length(result_st.data{resCols.subject_id}),1);
resp_select_val = zeros(length(result_st.data{resCols.subject_id}),1);

for isub = 1:nsub
    curr_sub = subids{isub};
    sub_idx = strmatch(curr_sub,subid_vect);
    sub_tap_trialVect{isub,1} = tapResp_st.data{tapRespCols.total_tapping}(sub_idx);
    sub_tap_subid{isub} = curr_sub;
    sub_tap_mean(isub,1) = mean(sub_tap_trialVect{isub});
    sub_tap_sd(isub,1) = std(sub_tap_trialVect{isub});
    sub_tap_total(isub,1) = sum(sub_tap_trialVect{isub}); 
    nZeroTapTrials = length(find(sub_tap_trialVect{isub} < 10));
    
    % find subs that do not meet criteria for minimum amount of tapping
    if nZeroTapTrials >= params.minTapCrit*41
        resp_select_lvl(sub_idx) = {'non_tapper'};
        nontapperIdx = nontapperIdx+1;
        nontapperSubs{nontapperIdx} = curr_sub; 
        nontapperSubVectIdx(nontapperIdx) = isub;
    end
end % for isub = 


% exclude nontapping subs from subsequent analyses
subids_nontapRmvd = subids;
subids_nontapRmvd(nontapperSubVectIdx) = [];
nontapRemoved_mean = sub_tap_mean;
nontapRemoved_mean(nontapperSubVectIdx) = [];
nontapRemoved_sd = sub_tap_sd;
nontapRemoved_sd(nontapperSubVectIdx) = [];
nontappers_mean = sub_tap_mean(nontapperSubVectIdx);
nontappers_sd = sub_tap_sd(nontapperSubVectIdx);


% fprintf('Excluding the following subjects from tapping analyses due to a lack of tapping:\n');
% for iexclude =1:length(nontapperSubs)
%     fprintf('%s\n.\n.\n.\n',nontapperSubs{iexclude});
% end

% for iNonTapSub = 1:length(nontapperSubs)
%     nontap_sub_idx = strmatch(nontapperSubs{iNonTapSub}, tapResp_st.data{tapRespCols.subject_id});
%     for idata = 1:length(tapResp_st.data)
%         result_st.data{idata}(nontap_sub_idx) = [];
%     end
% end



% calculate selectivity index for each subject by dividing st dev / mean
selectivity = sub_tap_sd./sub_tap_mean;
selectivity(isnan(selectivity)) = 0;
selectivity_nontappers = selectivity(nontapperSubVectIdx);
means_nontappers = sub_tap_mean(nontapperSubVectIdx);
sd_nontappers = sub_tap_sd(nontapperSubVectIdx);
selectivity_nontapRemvd = selectivity;
selectivity_nontapRemvd(nontapperSubVectIdx) = [];

% rank subjects' selectivity
[sorted_selectivity,sorted_select_idx] = sort(selectivity_nontapRemvd,'descend');
sorted_selectSubs = subids_nontapRmvd(sorted_select_idx);
sorted_means = nontapRemoved_mean(sorted_select_idx);
sorted_sd = nontapRemoved_sd(sorted_select_idx);
selectivity_lvl = cell(size(selectivity));
% now that tappers have been ranked on selectivity, append nontappers to
% bottom of the list
sorted_selectSubs = [sorted_selectSubs; nontapperSubs'];
sorted_selectivity = [sorted_selectivity; selectivity_nontappers];
sorted_means = [sorted_means; means_nontappers];
sorted_sd = [sorted_sd; sd_nontappers];

for iselectlvl = 1:length(selectivity_lvl)
    if iselectlvl > length(sorted_select_idx)
        sel_lvl = 'non_tapper';
    elseif iselectlvl <=24
        sel_lvl = 'high_selectivity';
    else
        sel_lvl = 'low_selectivity';
    end
    selectivity_lvl{iselectlvl} = sel_lvl;
end

for iSub = 1:length(sorted_selectSubs)
    currSub = sorted_selectSubs{iSub};
    subidx = ismember(result_st.data{resCols.subject_id},currSub);
    resp_select_val(subidx) = sorted_selectivity(iSub);
    resp_select_lvl(subidx) = selectivity_lvl(iSub);
end

selectivity_mtx = {sorted_selectSubs selectivity_lvl sorted_selectivity};

% save selectivity info to .mat file for use in other analysis scripts
select_fname = fullfile(params.matpath, 'selectivity_mtx');
fprintf('Saving selectivity_mtx to %s\n', select_fname);
save_vars = 'selectivity_mtx';
save(select_fname, save_vars)

% Print selectivity info to file
selectivity_table = ensemble_init_data_struct;
selectivity_table.vars = {'data','column_labels','column_formats'};
selectivity_table.data{1} = selectivity_mtx;
selectivity_table.data{2} = {'subject_id' 'selectivity_level' 'selectivity_score'};
selectivity_table.data{3} = {'%s','%s','%1.2f'};
ensemble_display_table(selectivity_table,params.report_selectivity);

%% Calculate & sort sMESS movement subscores

fprintf('\nPulling sMESS data for %d subjects:\n', nsub)

curr_parms.extract_var = 'response_enum';
ensemble_resps = ensemble_extract_matrix(qid_st, curr_parms);
ensemble_resps = enum2data(ensemble_resps);  % convert data to scale values
sMESS_scores = sum(ensemble_resps,2);

curr_parms.extract_var = 'subject_id';
sMESS_subid =  ensemble_extract_matrix(qid_st, curr_parms);
sMESS_subid(:,2:5) = [];
% nsub_sMESS = length(sMESS_subid);



% sort sMESS movement subscores into same order as selectivity stats
sMESS_orderedScores = nan(size(sorted_selectSubs));
res_sMESS = nan(length(result_st.data{resCols.subject_id}),1);
for isub = 1:nsub
    curr_sMESSsub = sorted_selectSubs{isub};
    tapStatSub_idx = strmatch(curr_sMESSsub,sMESS_subid);
    sMESS_orderedScores(isub) = sMESS_scores(tapStatSub_idx);
    resSub_idx = ismember(result_st.data{resCols.subject_id},curr_sMESSsub);
    if ~isempty(resSub_idx) % some subs were removed due to nontapping
        res_sMESS(resSub_idx) = sMESS_scores(isub);
    end        
end
%% get subject stats

% get subjects' musical experience
[musExp_subIDs yrsExperience instruments] = ensemble_extract_music_experience(resp_st,params.sub_musExp);
% some subs had less than 1 yr experience. Remove them from these arrays.
zeroExpSubs = find(yrsExperience==0,1);
if ~isempty(zeroExpSubs)
   musExp_subIDs(zeroExpSubs) = [];
   yrsExperience(zeroExpSubs) = [];
   instruments(zeroExpSubs) = [];
end
nsubs_experience = length(musExp_subIDs);
range_experience = [min(yrsExperience) max(yrsExperience)];
mean_experience = mean(yrsExperience);
std_experience = std(yrsExperience);

% task-related self consciousness ratings
selfconscious_beginning_msk = ismember(resp_st.data{resp_cols.question_id},606);
selfconscious_beginning = cell(1,2);
selfconscious_beginning{1} = resp_st.data{resp_cols.subject_id}(selfconscious_beginning_msk);
selfconscious_beginning{2} = enum2data(resp_st.data{resp_cols.response_enum}(selfconscious_beginning_msk));
selfconscious_end_msk = ismember(resp_st.data{resp_cols.question_id},607);
selfconscious_end = cell(1,2);
selfconscious_end{1} = resp_st.data{resp_cols.subject_id}(selfconscious_end_msk);
selfconscious_end{2} = enum2data(resp_st.data{resp_cols.response_enum}(selfconscious_end_msk));
% need to reverse score self-consciousness ratings, as they were recorded on a reversed scale (1 =
% "always" self-conscious, 5 = "never" self-conscious)
selfconscious_beginning{2} = 5 + 1 - selfconscious_beginning{2};
selfconscious_end{2} = 5 + 1 - selfconscious_end{2};

% normal hearing ratings
normal_hearing_mask = ismember(resp_st.data{resp_cols.question_id},1);
normal_hearing = cell(1,2);
normal_hearing{1} = resp_st.data{resp_cols.subject_id}(normal_hearing_mask);
normal_hearing{2} = resp_st.data{resp_cols.response_enum}(normal_hearing_mask);
normal_hearing{2}(ismember(normal_hearing{2},2)) = 0;

% responses to question of whether subject has heard the term "groove"
% applied to music
everHeardGrv_msk = ismember(resp_st.data{resp_cols.question_id},14);
everHeardGrv = cell(1,2);
everHeardGrv{1} = resp_st.data{resp_cols.subject_id}(everHeardGrv_msk);
everHeardGrv{2} = enum2data(resp_st.data{resp_cols.response_enum}(everHeardGrv_msk));
rmv_grv_def_msk = ~isnan(everHeardGrv{2}); % mask out groove definition text responses
everHeardGrv{1} = everHeardGrv{1}(rmv_grv_def_msk);
everHeardGrv{2} = everHeardGrv{2}(rmv_grv_def_msk);
everHeardGrv{2}(ismember(everHeardGrv{2},2)) = 0; % convert "no" responses from 2 to 0


substat_fname = '/data/polystream_move/tables/subject_stats.txt';
fid = fopen(substat_fname,'w');
fprintf('\nWriting subject stats to the following file: %s\n\n',substat_fname);
fprintf(fid,'%d out of %d subjects reported instrument or voice experience\nyears experience range:%d-%d\nmean:%1.2f years\nSD:%1.2f years\n',...
    nsubs_experience,nsub,range_experience(1),range_experience(2),mean_experience,std_experience);
fprintf(fid,'%d out of %d subjects reported to have normal hearing\n',sum(normal_hearing{2}),length(normal_hearing{2}));
fprintf(fid,'%d out of %d subjects reported to have heard the term groove applied to music\n',sum(everHeardGrv{2}),length(everHeardGrv{2}));

% assign subject stats to result_st output
% plays_inst = nan(size(subid_vect));
yrs_training = nan(size(subid_vect));
selfconscious_firsthalf = nan(size(subid_vect));
selfconscious_lasthalf = nan(size(subid_vect));
expt_block = nan(size(subid_vect));
ever_heard_grv = cell(size(subid_vect));
mus_trained = cell(size(subid_vect));
subStats_yrsTraining = nan(size(sorted_selectSubs));
subStats_selfconsc_firsthalf = nan(size(sorted_selectSubs));
subStats_selfconsc_lasthalf = nan(size(sorted_selectSubs));
subStats_everHeardGrv = cell(size(sorted_selectSubs));
subStats_musTrained = cell(size(sorted_selectSubs));
block_structure = [ones(20,1);2*ones(21,1)];

for isubj = 1:nsub
    curr_subj = sorted_selectSubs{isubj};
    curr_sub_resultMsk = ismember(subid_vect,curr_subj);
    expt_block(curr_sub_resultMsk) = block_structure;
    curr_sub_yrs_trainingMask = ismember(musExp_subIDs,curr_subj);
    curr_sub_everHeardGrv = everHeardGrv{2}(ismember(everHeardGrv{1},curr_subj));
    if curr_sub_everHeardGrv == 1
        curr_sub_everHeardGrv = 'yes';
    else
        curr_sub_everHeardGrv = 'no';
    end
    
    if ~sum(curr_sub_yrs_trainingMask)
        yrs_training(curr_sub_resultMsk) = 0;
        subStats_yrsTraining(isubj) = 0;
    else
        yrs_training(curr_sub_resultMsk) = yrsExperience(curr_sub_yrs_trainingMask);
        subStats_yrsTraining(isubj) = yrsExperience(curr_sub_yrs_trainingMask);
    end
    
    % 0 to 1 yrs training = 'non-trained'; 2+ yrs training = 'trained'
    if subStats_yrsTraining(isubj) <= 1
        mus_trained(curr_sub_resultMsk) = {'non-trained'};
        subStats_musTrained{isubj} = 'non-trained';
    else
        mus_trained(curr_sub_resultMsk) = {'trained'};
        subStats_musTrained{isubj} = 'trained';
    end
    currsub_selfconsc_frsthlf_msk = ismember(selfconscious_beginning{1},curr_subj);
    selfconscious_firsthalf(curr_sub_resultMsk) = selfconscious_beginning{2}(currsub_selfconsc_frsthlf_msk);
    subStats_selfconsc_firsthalf(isubj) = selfconscious_beginning{2}(currsub_selfconsc_frsthlf_msk);
    currsub_selfconsc_lsthlf_msk = ismember(selfconscious_end{1},curr_subj);
    selfconscious_lasthalf(curr_sub_resultMsk) = selfconscious_end{2}(currsub_selfconsc_lsthlf_msk);
    subStats_selfconsc_lasthalf(isubj) = selfconscious_end{2}(currsub_selfconsc_lsthlf_msk);
    
    ever_heard_grv(curr_sub_resultMsk) = {curr_sub_everHeardGrv};
    subStats_everHeardGrv{isubj} = curr_sub_everHeardGrv;
end

%% print & output results

% concatenate response data into a single matrix
resp_mtx = {sorted_selectSubs selectivity_lvl sorted_selectivity ...
    subStats_yrsTraining subStats_musTrained sMESS_orderedScores subStats_selfconsc_firsthalf ...
    subStats_selfconsc_lasthalf subStats_everHeardGrv};
resp_table.vars = {'data','column_labels','column_formats'};
resp_table.data{1} = resp_mtx;
resp_table.data{2} = {'subject_id' 'selectivity_lvl' 'selectivity_score' ...
    'years_training' 'musically_trained' 'sMESS_movement_subscore' 'selfconscious_firsthalf' ...
    'selfconscious_lasthalf' 'ever_heard_groove'};
resp_table.data{3} = {'%s','%s','%1.2f','%d','%s','%1.2f','%d','%d','%s'};

% print the table to file
if isfield(params,'report')
    ensemble_display_table(resp_table,params.report);
end

% output response struct with nontapping subjects removed & sub analysis calcs added
result_st.vars = [result_st.vars,{'selectivity_index'},{'selectivity_level'},...
    {'sMESS_movement_subscore'},{'years_training'},{'musically_trained'},...
    {'self_conscious_beginning'},{'self_conscious_near_end'},{'experiment_block'},...
    {'ever_heard_groove'}];
result_st.data = [result_st.data, {resp_select_val},{resp_select_lvl},{res_sMESS},...
    {yrs_training},{mus_trained},{selfconscious_firsthalf},{selfconscious_lasthalf},...
    {expt_block},{ever_heard_grv}];


% plot data
plot_subject_data;

return