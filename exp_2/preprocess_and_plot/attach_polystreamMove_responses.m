function result_st = attach_polystreamMove_responses(inData, params)

% Creates a struct containing preprocessed MIDI tapping responses and responses to Ensemble 
% post-stim questions across all subjects, and organizes data by question ID.
%
% Modeled after attach_polystream_groove_responses.m
%
% 2/2012 - BH

clear anSearchCrit;
anSearchCrit.name = 'session_x_trial_info';
anListIdx = ensemble_find_analysis_struct(inData,anSearchCrit);
sessinfo = inData{anListIdx};

clear anSearchCrit;
anSearchCrit.name = 'response_data';
anListIdx = ensemble_find_analysis_struct(inData,anSearchCrit);
resp_st = inData{anListIdx};
% resp_st.vars{end} = 'compqid';
% resp_st.data{end} = params.compqid;

% midi response vars
resp_vars = {'resamp_midiResp','resamp_midiSec','cumsum_midiResp','total_tapping'};
nvars = length(resp_vars);

result_st = ensemble_init_data_struct;
result_st.name = 'poststim_responses';

% Generate a data structure that is oriented by composite question IDs
qid_st = ensemble_data_by_question(resp_st,params);

% Get a list of subjects
curr_parms.extract_var = 'subject_id';
subid_vect = ensemble_extract_matrix(qid_st, curr_parms);
subid_vect(:,2:3) = [];  % We don't need the second column
subids = unique(subid_vect);
nsub = length(subids);
    
% Get a list of training stim IDs
practice_attrib = 'polystream_training';
[practice_id,] = mysql_get_stim_by_attribute(...
  'attrib_name', practice_attrib, ...
  'mysql', params.mysql);

% Get a vector of stimulus IDs
curr_parms.extract_var = 'stimulus_id';
stimid_vect = ensemble_extract_matrix(qid_st, curr_parms);
stimid_vect(:,2:3) = [];  % We don't need the second column
practice_mask = ismember(stimid_vect, practice_id{1});
stimid_vect(practice_mask) = []; % remove practice stim IDs
stimids = unique(stimid_vect);

% remove rows of subid_vect associated w/ practice trials
subid_vect(practice_mask) = [];

% count up how many times each stimulus occurred 
count = hist(stimid_vect, stimids);

% Ignore those stims for which we have too few responses
bad_stimids = stimids(count<params.minNumTrials);
stimids(count<params.minNumTrials) = [];
badstim_vect = ismember(stimid_vect, bad_stimids);

nstims = length(stimids);

% Extract the poststim data
curr_parms.extract_var = 'response_enum';
ensemble_resps = ensemble_extract_matrix(qid_st, curr_parms);
ensemble_resps = enum2data(ensemble_resps);  % convert data to scale values
ensemble_resps(practice_mask,:) = []; % remove rows associated with practice trials

% Eliminate entries from all vectors that correspond to bad stims
subid_vect(badstim_vect) = [];
stimid_vect(badstim_vect) = [];
ensemble_resps(badstim_vect,:) = [];
badstim_vect(badstim_vect) = [];

nobs = length(stimid_vect);

%
% Initialize the output variables
%
groove_data = cell(1,nvars);
for ivar = 1:nvars
    curr_var = resp_vars{ivar};
    switch curr_var
        case 'total_tapping'
            groove_data{ivar} = ones(nobs,1)*NaN;
        otherwise
            groove_data{ivar} = cell(nobs,1);
    end % switch curr_var
end

exptname = params.ensemble.expname;

fprintf('Pulling data for %d subjects:\n', nsub)
 
sessinfoCols = set_var_col_const(sessinfo.vars);
sessIDs = sessinfo.data{sessinfoCols.session_id};
resp_subids = sessinfo.data{sessinfoCols.subject_id};

for iSess = 1:length(sessIDs)
    curr_sub = subids{iSess};
    fprintf('%s\n',curr_sub);
    sub_idx = strmatch(curr_sub,resp_subids); % index subject in response struct
    submask = ismember(subid_vect,curr_sub);
    trialinfoCols = set_var_col_const(sessinfo.data{sessinfoCols.trial_info}{sub_idx}.vars);
    resp_stim_id = sessinfo.data{sessinfoCols.trial_info}{sub_idx}.data{trialinfoCols.stimulus_id};
    midiResp_st = sessinfo.data{sessinfoCols.trial_info}{sub_idx}.data{trialinfoCols.midi_resp};
    
    % Loop over the stimulus IDs and find the corresponding response data
    for istim = 1:nstims
        curr_stimid = stimids(istim);
        stimmask = stimid_vect == curr_stimid;
        
        if sum(submask & stimmask) > 1
            msgstr = sprintf('Too many instances of stimid (%d) for subject (%s)', curr_stimid, curr_sub);
            warning(msgstr);
            continue
        end
        
        % Find this stimulus in the list of responses
        resp_idx = find(resp_stim_id == curr_stimid);
        if isempty(resp_idx) 
            continue
        end
        
        midiRespCols = set_var_col_const(midiResp_st{resp_idx}.vars);
        
        % Extract the response variables & assign the data
        for ivar = 1:nvars
        
            curr_var = resp_vars{ivar};
            groove_data_idx = strmatch(curr_var, resp_vars, 'exact');
            
            switch curr_var
                case 'total_tapping'
                    curr_data = midiResp_st{resp_idx}.data{midiRespCols.(curr_var)};
                otherwise
                    curr_data = {midiResp_st{resp_idx}.data{midiRespCols.(curr_var)}};
            end
            
            % Assign the data
            groove_data{groove_data_idx}(submask&stimmask) = curr_data;
        
        end % for ivar
    end % for istim
end % for iSess

% Finalize output data struct
new_labels = params.qid_labels(:,2)';
result_st.vars = [{'subject_id','stimulus_id'}, new_labels, resp_vars];
result_cols = set_var_col_const(result_st.vars);
for ivar = 1:length(result_st.vars)
  curr_var = result_st.vars{ivar};
  switch curr_var
    case 'subject_id'
      result_st.data{ivar} = subid_vect;
    case 'stimulus_id'
      result_st.data{ivar} = stimid_vect;
    case new_labels
      label_idx = strmatch(curr_var,new_labels);
      idx = find([qid_st.meta.question.compqid] == ....
        params.qid_labels{label_idx,1}); % arghh!!!!
      result_st.data{ivar} = ensemble_resps(:,idx);
      result_st.meta.question(label_idx) = qid_st.meta.question(idx);
    case resp_vars
      groove_data_idx = strmatch(curr_var, resp_vars,'exact');
      result_st.data{ivar} = groove_data{groove_data_idx};
  end % switch
end % for ivar

return