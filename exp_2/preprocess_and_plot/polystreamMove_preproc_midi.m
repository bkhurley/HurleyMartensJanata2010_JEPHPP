function sessinfo = polystreamMove_preproc_midi(indata,params)

% polystreamMove_preprocess_mididata
%
% parses MIDI tapping from polystream_move experiment
% parses audio, matches parsed audio with trial stims from Ensemble
% resamples MIDI data to sample rate specified in params.sampleFreq
%

if params.use_saved_preproc
    load (fullfile(params.paths.matpath,'sessinfo_polystream_move'));

else    
    params.verbose = 2;
    
    % get respdata and sessdata    
    expinfo = indata;
    expcols = set_var_col_const(expinfo.vars);  % get indices into data array
    respdata = expinfo.data{expcols.response_data}; % extract response data
    sessdata = expinfo.data{expcols.session_info}; % extract session data
    
    % parse audio and MIDI, associate with trials from Ensemble session
    sessinfo = ensemble_session_trial_info_v2({respdata, sessdata},params);
end

% Save the parsed data (sessinfo) so future changes to code below this
% point won't require re-running the long parsing routine
sessfname = fullfile(params.paths.matpath, sprintf('sessinfo_%s', params.experiment_name));
fprintf('Saving sessinfo to %s\n', sessfname);
save_vars = 'sessinfo';
save(sessfname, save_vars)

% resample midi responses onto sample-time grid
sampling_rate = params.resamp.sampling_rate;

% Load midi toolbox default parameters
miditoolbox_params;

fprintf('Resampling MIDI responses at %d Hz\n', sampling_rate);

sessCols = set_var_col_const(sessinfo.vars);

for isub = 1:length(sessinfo.data{sessCols.trial_info})
    trial = sessinfo.data{sessCols.trial_info}{isub};
    trialCols = set_var_col_const(trial.vars);
    midiResp = trial.data{trialCols.midi_resp};
    
    for itrial = 1:length(midiResp)
        thisMidiResp = sessinfo.data{sessCols.trial_info}{isub}.data{trialCols.midi_resp}{itrial};
        midirespCols = set_var_col_const(thisMidiResp.vars);
        
        % Resample stimulus vector onto uniform scale
        startTimeSec = trial.data{trialCols.start_time_sec}(itrial);
        durSec = trial.data{trialCols.dur_sec}(itrial);
        stopTime = startTimeSec + durSec;
        resamp_tscale = startTimeSec:1/sampling_rate:stopTime;
        this_nmat = thisMidiResp.data{midirespCols.nmat};
        resamp_midiResp = histc(this_nmat(:,ONSET_SEC_COL),resamp_tscale);
        resamp_midiSec = resamp_tscale;
        resampMatSize = size(resamp_midiResp);
        % for some reason, one subject's resp matrix dimension was different from others.
        % This will conform matrices to common dimension
        if resampMatSize(1) > 1
            resamp_midiResp = resamp_midiResp';
            resamp_midiSec = resamp_midiSec';
        end
        
        % append resampled tapping to sessinfo struct
        sessinfo.data{sessCols.trial_info}{isub}.data{trialCols.midi_resp}{itrial}.vars{end+1} = 'resamp_midiResp';
        sessinfo.data{sessCols.trial_info}{isub}.data{trialCols.midi_resp}{itrial}.data{end+1} = resamp_midiResp;
        
        sessinfo.data{sessCols.trial_info}{isub}.data{trialCols.midi_resp}{itrial}.vars{end+1} = 'resamp_midiSec';
        sessinfo.data{sessCols.trial_info}{isub}.data{trialCols.midi_resp}{itrial}.data{end+1} = resamp_midiSec;
        
        sessinfo.data{sessCols.trial_info}{isub}.data{trialCols.midi_resp}{itrial}.vars{end+1} = 'cumsum_midiResp';
        sessinfo.data{sessCols.trial_info}{isub}.data{trialCols.midi_resp}{itrial}.data{end+1} = cumsum(resamp_midiResp);
        
        sessinfo.data{sessCols.trial_info}{isub}.data{trialCols.midi_resp}{itrial}.vars{end+1} = 'total_tapping';
        sessinfo.data{sessCols.trial_info}{isub}.data{trialCols.midi_resp}{itrial}.data{end+1} = sum(resamp_midiResp);
    end %itrial
    
end %isub