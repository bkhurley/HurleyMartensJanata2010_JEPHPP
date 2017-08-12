%
% Based on preproc_zebris_v2 by Stefan Tomic
%
% This script is being developed to preprocess mocap data specifically from
% the polystream_move experiment, as the original preproc_zebris_v2 had
% various nonfunctioning/incomplete parts. The main difference is 
% mocap_matchWithAudio_v2. This script should eventually be
% generalized as standard script for preprocessing any motion data. 
%
%
% Call this script from your Ensemble analysis script and
% pass parameters to it via your global parameters. The script only
% adds analysis structs to your analysis stream for preprocessing
% motion capture data.
%
% Your analysis script may look something like:
%       globals = polystreamMove_globals;
%       polystreamMove_preproc_zebris;
%       na = na + 1;
%       an{na}.name = 'my_project_func'
%       an{na}.fun = @my_project_func;
%       params.run_analyses = [1 2 3 8];
%       an = ensemble_jobman(an,params);
%       etc.
%
%       my_project_func.m would be a function specific to your
%       project that you would like call after the preprocessing
%       functions.


%if you want to log to a file, just set the path in your globals.logFilepath
try
  logFileID = fopen(globals.preProcLogFile,'w+');
catch
  logFileID = 1;
end


na = 0;

na = na + 1;
% Params:
% zebrisFileList: A cell array of ascii file names to read in.
% sampleFreq: The sample frequency of all input data in Hz.
% srcPath: path from which to read files specified in zebrisFileList
% destPath: path to save mat files
% ignoreSeg.shorterThanSec: ignore segments that are shorter than
%                           this duration in seconds.
% ignoreSeg.indicesInFile: an X-by-2 cell array, with filestubs in
%                          the first column, and vectors in the second column. The second
%                          column specifies the segment indices to ignore for that filestub
%                          (which generally corresponds to a subject or session).
%
an{na}.name = 'convert_zebris2mocap';
an{na}.fun = @convert_zebris2mocap;
an{na}.params.recordedDataFiles = globals.recordedDataFiles;
an{na}.params.markerLabels = globals.zebris.markerLabels;
an{na}.params.sampleFreq = globals.zebris.sampleFreq;
an{na}.params.srcPath = globals.paths.zebrisAscii;
an{na}.params.destPath = globals.paths.mocapMatFiles;
an{na}.params.ignoreSeg.shorterThanSec = globals.zebris.ignoreSeg.shorterThanSec;
an{na}.params.logFileID = logFileID;


na = na + 1;
an{na}.name = 'mocap_matchWithAudio';
an{na}.fun = @mocap_matchWithAudio;
an{na}.requires = {struct('name','convert_zebris2mocap')};
an{na}.params.recordedDataFiles = globals.recordedDataFiles;
an{na}.params.markerSyncPath = globals.paths.recordedMarkerSync;
an{na}.params.musicPath = globals.paths.recordedMusic;
an{na}.params.ignoreSeg.shorterThanSec = globals.zebris.ignoreSeg.shorterThanSec;
an{na}.params.use_saved_parsing = 1;
an{na}.params.preproc_mat_file = fullfile(globals.paths.matpath,'sessinfo_polystream_move.mat');
an{na}.params.freq = 100;
%need an ignore idxs parameter for marker sync sig

na = na + 1;
an{na}.name = 'despike_zebris_data';
an{na}.fun = @despike_data_v2;
an{na}.requires = {struct('name','mocap_matchWithAudio')};
an{na}.params.logFileID = logFileID;
an{na}.params.use_consensus = 1;
an{na}.params.remove_stragglers = 1;
an{na}.params.check_nan_lead_edge = 1;
an{na}.params.amp_thresh = 400;  % max deviation (in mm)
an{na}.params.vel_thresh = 5;  % velocity threshold in mm/sample
an{na}.params.acc_thresh = 10;
an{na}.params.dev_from_start_val = 50;
an{na}.params.term_val_thresh = 10;
an{na}.params.search_win_samps = 25;
an{na}.params.backward_nan_search = 13;
an{na}.params.min_nan_gap = 3;
an{na}.params.nan_lead_thresh = 20;
an{na}.params.dev_idx_method = 'velocity_acceleration'; % 'velocity_acceleration', 'velocity_only'
an{na}.params.stdev_thresh = 4; % Chebyshev's inequality may be relevant.


na = na + 1;
% Params:
% cut_of_nan: Markers that have more than this percentage of NaNs in
%             a time series will just be removed altogether (too much missing data)
an{na}.name = 'spline_interpolation';
an{na}.fun = @spline_interpolation_v2;
an{na}.requires = {struct('name','despike_zebris_data')};
an{na}.params.cut_off_nan = globals.spline_interp.cut_off_nan;
an{na}.params.signal_loss_info_fname = globals.spline_interp.signal_loss_info_fname;
an{na}.params.logFileID = logFileID;
an{na}.params.mcSigs_figpath = globals.paths.zebrisFigs;
an{na}.params.plot_sigs = 0;


% na = na + 1;
% an{na}.name = 'manual_adjustment';
% an{na}.fun = @polystreamMove_manual_adjust_zebris_data;
% an{na}.requires = {struct('name','spline_interpolation')};
% an{na}.params.reportStats = globals.manual_adjustment.reportStats; 
% an{na}.params.ensemble.conn_id = globals.ensemble.conn_id;

na = na + 1;
% Params:
% reportStats: whether or not statistics reporting percentage and
%              number of deletions will be reported.Possible values are 'yes'
%              (automatically report), 'no' (do not report), or 'ask' (user will
%              be prompted if they want a report after each data set).
an{na}.name = 'manual_adjustment';
an{na}.fun = @manual_adjust_zebris_data;
an{na}.requires = {struct('name','spline_interpolation')};
an{na}.params.reportStats = globals.manual_adjustment.reportStats; 
an{na}.params.ensemble.conn_id = 2;


na = na + 1;
an{na}.name = 'make_motion_movie_v2';
an{na}.fun = @make_motion_movie_v2;
an{na}.requires = {struct('name','convert_zebris2mocap')};
an{na}.params.mapar = globals.make_movie.mapar;
an{na}.params.connList = globals.zebris.connList;
an{na}.params.outputPath = globals.paths.moviePath;


na = na + 1;
an{na}.name = 'motion_decompose_and_scale';
an{na}.fun = @motion_decompose_and_scale_v2;
an{na}.requires = {struct('name','spline_interpolation')};
an{na}.params.connList = globals.zebris.connList;
an{na}.params.perform = 'decompose';

na = na + 1;
an{na}.name = 'motion_highpass_filt';
an{na}.fun = @motion_highpass_filt;
an{na}.requires = {struct('name','motion_decompose_and_scale')};
an{na}.params.paths.highPassFilt = globals.paths.highPassFilt;
an{na}.params.hp_cutoff = .5;
an{na}.params.filt_data = 1; % set to 0 if one wishes to bypass filtering
an{na}.params.ipem = globals.ipem;
an{na}.params.connList = globals.zebris.connList;
an{na}.params.conns_to_analyze = [1]; % Useful if only analyzing a subset of marker connections. Multiple connections should be treated as a vector.

if(logFileID ~= 1)
  fclose(logFileID);
end

return



