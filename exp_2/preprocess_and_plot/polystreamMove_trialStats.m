function stim_hist = polystreamMove_trialStats(inData,params)

% generates output reporting amount that each subject tapped for each
% stimulus. This is then easily plotted as a histogram (as of now, it is
% plotted in polystreamMove_plots_v2.m)
%
% 05/3/12 BH

clear crit
crit.name = 'poststim_responses';
idx = ensemble_find_analysis_struct(inData,crit);
resp_st = inData{idx};
respcols = set_var_col_const(resp_st.vars);

clear crit
crit.name = 'stimulus_metadata';
idx = ensemble_find_analysis_struct(inData,crit);
stimdata_st = inData{idx};
stimdata_cols = set_var_col_const(stimdata_st.vars);


stimids = resp_st.data{respcols.stimulus_id};
tapping = resp_st.data{respcols.total_tapping};
stim_list = unique(stimids);
nstims = length(stim_list);

stim_hist = ensemble_init_data_struct;
stim_hist.name = 'trial_stats';
stim_hist.vars = {'stimulus_id' 'stimulus_name' 'data'};

% initialize data cells
stim_hist.data{1} = cell(nstims,1);
stim_hist.data{2} = cell(nstims,1);
stim_hist.data{3} = cell(nstims,1);

for istim = 1:nstims
    
    % get location of current stim in response data
    curr_stimid = stim_list(istim);
    stimtap_idx = stimids == curr_stimid;
    % get location of current stim in stimulus_metadata (in order to know
    % stim name)
    stimname_idx = stimdata_st.data{stimdata_cols.stimulus_id} == curr_stimid;
    curr_stimName = stimdata_st.data{stimdata_cols.name}(stimname_idx);
    
    stim_hist.data{1}{istim} = curr_stimid;
    stim_hist.data{2}{istim} = curr_stimName;
    stim_hist.data{3}{istim} = tapping(stimtap_idx); 

end

return
