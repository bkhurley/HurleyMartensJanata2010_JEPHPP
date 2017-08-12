function outData = polystreamMove_midiResp_rhythm_profiler(inData,params)
%
% Calculates summary measurements of tapping from the rhythm model
%
% INPUT(S):
%   session_x_trial structure, returned from ensemble_session_x_trial_info.
%   stimulus_metadata, returned from ensemble_load_expinfo
%   response_w_sub_analysis, to conform output to polystream analysis struct
%
% PARAM(S):
%  params.process:       'calcMeans','calcRP','loadRP'
%  params.savePath:      path to save individual RPs and summary measurement mat file
%  params.ipem_proc:     ipem_proc_series parameters for rp calculation on tapping data
%  params.storeIpemProc: cell array of strings designating which result of ipem_proc_series to
%                        save for each trial (e.g. {'rp'})
%  params.filt:          ensemble_filter parameters for filtering out
%                        sessions/subs to process
%  params.fname:         Filename to save summary measurement struct (saved to params.savePath)
%
%
% DESCRIPTION:
%  Loops through each session and trial in parallel loops, converts the midi
%  responses from each trial to an aud struct. Each aud struct
%  is then processed with jlmt_proc_series. The RPs for each
%  trial are saved to disk in separate files at params.savePath.
%  
%
%
%  03/29/2013 BH - Adapted code from midiResp_rhythm_profiler (written by
%                  Stephan Tomic)

storeIpemProc = params.storeIpemProc;

if(~isfield(params,'process'))
    params.process = {'calcRP'};
end

if(length(intersect({'calcRP','loadRP'},params.process)) == 0)
    params.process = union(params.process,{'calcRP'});
end

clear searchCrit;
searchCrit.name = 'stimulus_metadata';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
stimulusMetadata = inData{findIdx};
smCols = set_var_col_const(stimulusMetadata.vars);

clear searchCrit;
searchCrit.name = 'session_x_trial_info';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
sessTrials = inData{findIdx};
sessTrialsCols = set_var_col_const(sessTrials.vars);

clear searchCrit;
searchCrit.name = 'response_w_sub_analysis';
findIdx = ensemble_find_analysis_struct(inData,searchCrit);
resp_st = inData{findIdx};
respCols = set_var_col_const(resp_st.vars);
stimid_list = resp_st.data{respCols.stimulus_id};
subid_list = resp_st.data{respCols.subject_id};
subids = unique(subid_list);
tapRP = cell(length(subid_list),1);

outData = resp_st;
outData.name = 'midiResp_rhythm_profiler';
outData.vars = [outData.vars 'tappingRP'];
outDataCols = set_var_col_const(outData.vars);

%filter all trial info with filter settings in params.
if(isfield(params,'filt'))
    clear tmpFilt
    tmpFilt = params.filt;
    sessTrials = ensemble_filter(sessTrials,tmpFilt);
end

savepath = params.savePath;
jlmt_proc = params.process;
sessIDs = sessTrials.data{sessTrialsCols.session_id};

% Fire up parallel compute pool if possible
if exist('matlabpool') && ~matlabpool('size') 
    matlabpool
end

for iSess = 1:length(sessIDs)
    
    thisSubjectID = sessTrials.data{sessTrialsCols.subject_id}{iSess};
    % move on if current subject is non-tapper
    if isempty(strmatch(thisSubjectID,subids))
        continue
    end
    subMask = ismember(subid_list,thisSubjectID);
    thisSessionID = sessTrials.data{sessTrialsCols.session_id}(iSess);
    
    fprintf('Processing Subject %s\n',thisSubjectID);
    
    clear tmpFilt;
    tmpFilt.include.all.session_id = thisSessionID;
    thisSessionStruct = ensemble_filter(sessTrials,tmpFilt);
    
    thisSessTrials = thisSessionStruct.data{sessTrialsCols.trial_info}{1};
    trialInfoCols = set_var_col_const(thisSessTrials.vars);
    
    nTrials = length(thisSessTrials.data{trialInfoCols.stimulus_id});
    currSub_tappingRP = cell(nTrials,1);
    currSub_stimIDs = cell(nTrials,1);
    clear tmpFilt
    
    parfor iTrial = 1:nTrials % parfor
        
        filename = sprintf('%s_%d_trial_%d.mat',thisSubjectID,thisSessionID,iTrial);

        thisStimulusID = thisSessTrials.data{trialInfoCols.stimulus_id}(iTrial);
        thisStartTime = thisSessTrials.data{trialInfoCols.start_time_sec}(iTrial);
        thisStopTime = thisStartTime + thisSessTrials.data{trialInfoCols.dur_sec}(iTrial);
        
        thisTrialStruct = filt_ensemble_trial(thisSessTrials,thisStimulusID,thisStartTime);
        
        thisTrialStruct.vars = {thisTrialStruct.vars{:} 'resp_aud' storeIpemProc{:}};
        thisTrialStructCols = set_var_col_const(thisTrialStruct.vars);

        
        if(ismember('calcRP',jlmt_proc))
            disp(sprintf('JLMT Proc for session %d/%d, trial %d/%d\n',iSess,length(sessIDs),iTrial,nTrials));
            
            thisMidiResp = thisTrialStruct.data{thisTrialStructCols.midi_resp}{1};
            midiRespCols = set_var_col_const(thisMidiResp.vars);
            
            thisNmat =thisMidiResp.data{midiRespCols.nmat};
            
            if(iscell(thisNmat))
                thisNmat = thisNmat{1};
            end
            
            % trials with only 1 midi event will error out in RP model
            % (i.e., no periodicities to analyze)
            % this will exclude such trials from further RP-related analyses
            if size(thisNmat,1) <= 1
                currSub_stimIDs{iTrial,1} = thisStimulusID;
                currSub_tappingRP{iTrial,2} = [''];
                continue
            end
            
            nmat2aud_params = populate_nmat2aud_params(thisStartTime,thisStopTime,params);
            respAud = nmat2aud(thisNmat,nmat2aud_params);
            audCols = set_var_col_const(respAud.vars);
            respAud.data{audCols.filename} = filename;
            respAud.data{audCols.path} = savepath;
            thisTrialStruct.data{thisTrialStructCols.resp_aud}{1} = respAud;

            jlmtOut = jlmt_proc_series(respAud,params.ipem_proc);
            jlmtOutCols = set_var_col_const(jlmtOut.vars);
            thisRP = jlmtOut.data{jlmtOutCols.rp}{1};
            
            currSub_stimIDs{iTrial,1} = thisStimulusID;
            currSub_tappingRP{iTrial,2} = thisRP;

        end
        
    end
   
    for itrial = 1:nTrials
        curr_stim = currSub_stimIDs{itrial};
        stimidMask = ismember(stimid_list,curr_stim);
        sub_stim_mask = subMask&stimidMask;
        currRP = currSub_tappingRP{itrial};
        tapRP{sub_stim_mask} = currRP;
    end
    
end

% terminate parallel worker pool
matlabpool close

outData.data{outDataCols.tappingRP} = tapRP;

end

%% Subfunctions

function nmat2aud_params = populate_nmat2aud_params(startTime,stopTime,params) 
nmat2aud_params = params.nmat2aud;
nmat2aud_params.startTime = startTime;
nmat2aud_params.stopTime  = stopTime;
end

function thisTrialStruct = filt_ensemble_trial(thisSessTrials,stimulusID,startTime)
tmpFilt.include.all.stimulus_id = stimulusID;
tmpFilt.include.all.start_time_sec = startTime;
thisTrialStruct = ensemble_filter(thisSessTrials,tmpFilt);
end