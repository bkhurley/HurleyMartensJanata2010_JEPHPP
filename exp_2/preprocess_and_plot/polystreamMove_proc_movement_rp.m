function outData = polystreamMove_proc_movement_rp(inData,params)

% Processes motion capture data through Rhythm Profiler model, attaches
% session & trial info
%
% Nov 7, 2012 BH

clear searchCrit;
searchCrit.name = 'motion_highPass_filt';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
mocapData_st = inData{findIdx};
mocapCols = set_var_col_const(mocapData_st.vars);

clear searchCrit;
searchCrit.type = 'stimulus_metadata';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
stimMetaData = inData{findIdx};
stimMetaDataCols = set_var_col_const(stimMetaData.vars);
stimNames = stimMetaData.data{stimMetaDataCols.name};
stimIDs = stimMetaData.data{stimMetaDataCols.stimulus_id};
nstims = length(stimIDs);

% process mocap RP
mocapResp_st = mocap_procMetricSurface_v2(mocapData_st,params);
mocapRespCols = set_var_col_const(mocapResp_st.vars);
mocapStimList = mocapResp_st.data{mocapRespCols.stim_id};
mocapSubList = mocapResp_st.data{mocapRespCols.subject_id};
mocapRP = mocapResp_st.data{mocapRespCols.mocapRP};

%% Attach stim & attribute info to mocap data

% Get stim names
stimnameList = cell(size(mocapStimList));
for istim = 1:nstims
    curr_stimMetaIdx = mocapStimList == stimIDs(istim);
    [stimnameList{curr_stimMetaIdx}] = deal(stimNames{istim});
end

% Get genre
genres = fieldnames(params.songatt.genre);
ngenres = length(genres);
genreList = cell(size(stimnameList));
for igenre = 1:ngenres
    curr_genre = genres{igenre};
    [genreList{ismember(stimnameList,params.songatt.genre.(curr_genre))}] = deal(curr_genre);
end

% Get entrance condition
sequence_levels = fieldnames(params.songatt.seq);
seqList = cell(size(stimnameList));
nseq = length(sequence_levels);
for iseq = 1:nseq
    curr_seq = sequence_levels{iseq};
    [seqList{ismember(stimnameList,params.songatt.seq.(curr_seq))}] = deal(curr_seq);
end

% Deal missing values
[mocapRP{cellfun(@isempty,mocapRP)}] = deal(NaN);

% Deal with chameleon
chameleon_idx = ismember(stimnameList,'chameleon.mp3');
[genreList{chameleon_idx}] = deal('chameleon');
[seqList{chameleon_idx}] = deal(NaN);

%% Assign output data
outData = ensemble_init_data_struct;
outData.vars = {'subject_id' 'stim_id' 'stim_name' 'genre' 'seq' 'mocapRP'};
outData.data{1} = mocapSubList;
outData.data{2} = mocapStimList;
outData.data{3} = stimnameList;
outData.data{4} = genreList;
outData.data{5} = seqList;
outData.data{6} = mocapRP;
