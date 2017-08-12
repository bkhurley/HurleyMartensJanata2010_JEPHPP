function globals = polystreamMove_globals

% returns global params for polystream_move

% Dec 2011      BH  adapted from polystream_globals_v2 in order to conform to 
%                   polystream_move 
%
% Jan 12, 2011  BH  added Chameleon to stim attributes
%
% Feb 10, 2011  BH  added preprocessing params

globals.version = mfilename;

%% expt specific params
globals.filename = 'polystreamMove_globals';
globals.experiment_name = 'polystream_move';

globals.attrib_names = {'polyMove_simultaneous','polyMove_staggeredStrong','polyMove_staggeredWeak',...
    'polyMove_dance1','polyMove_dance2','polyMove_dance3','polyMove_dance4','polyMove_dance5',...
    'polyMove_folk1','polyMove_folk2','polyMove_folk3','polyMove_folk4','polyMove_folk5',...
    'polyMove_funk1','polyMove_funk2','polyMove_funk3','polyMove_funk4','polyMove_funk5',...
    'polyMove_rock1','polyMove_rock2','polyMove_rock3','polyMove_rock4','polyMove_rock5',...
    'polyMove_chameleon',...
    'polyMove_dance','polyMove_folk','polyMove_funk','polyMove_rock'};

globals.resptbl_name = sprintf('response_%s', globals.experiment_name);

%% database connectivity

% mysql params
[null,user] = unix('whoami');
if ~ismember(deblank(user),{'bkhurley','petr'})                       
    globals.mysql.login_type = 'subject';
else
    globals.mysql.login_type = 'researcher';
end

% log into Ensemble
globals.mysql = mysql_login(globals.mysql);
globals.mysql.conn_id = mysql_make_conn(globals.mysql);
globals.mysql.expname = globals.experiment_name;

% Ensemble params
globals.ensemble = globals.mysql;
globals.ensemble.experiment_title = globals.experiment_name;
globals.ensemble.training_stim_id = 6329:6326;

% List of form name definitions
[form_id_const, form_name_id_const_map, form_name_list] = make_form_name_defs(globals);
globals.ensemble.form_defs.form_id_const = form_id_const;
globals.ensemble.form_defs.form_name_id_const_map = form_name_id_const_map;
globals.ensemble.form_defs.form_name_list = form_name_list;


% question id maps
groove_id_map = {'MUSIC_GROOVED',493.01;...
    'ENJOYED_MUSIC',1048.01;...
    'LIKE_MUSIC_TO_CONTINUE',1049.01;...
    };


globals.groove_question_ids = id_map2struct(groove_id_map);

%% function handles
globals.fun.weight_fun = {'polystreamMove_stimselect_library',...
    'polystreamMove_weight_fun'};
globals.fun.pre_send = {'polystreamMove_stimselect_library',...
    'polystreamMove_presend'};

%% stim selection parameters
globals.selection.type = 'weighted';
globals.selection.replacement = 0;                               % turns weighting of stim to 0 once played (i.e., no replacement)

globals.return_stim_during_init = 1;                             % return stimuli during first call to function

%% weighting parameters

% sets weighting parameter for each attribute
for k=1:(length(globals.attrib_names)-4) % genre attribs not used on stim weighting
  globals.weights.(globals.attrib_names{k}) = 0;
end

%% File paths
globals.paths.project_root = fullfile('/data2/', globals.experiment_name);
globals.paths.analysis_path = fullfile(globals.paths.project_root, 'analyses');
globals.paths.log_path = fullfile(globals.paths.project_root, 'logfiles');
globals.paths.fig_path = fullfile(globals.paths.project_root, 'figures');
globals.paths.movement_figpath = fullfile(globals.paths.fig_path,'movement_figures');
globals.paths.stim_recording_path = fullfile(globals.paths.project_root, 'dp','audio');
globals.paths.resp_midi_path = fullfile(globals.paths.project_root, 'dp','midi');
globals.paths.recordedMarkerSync = fullfile(globals.paths.project_root,'dp','marker_sync');
globals.paths.recordedMusic = globals.paths.stim_recording_path;
globals.paths.matpath = fullfile(globals.paths.project_root,'matfiles');
globals.paths.tablepath = fullfile(globals.paths.project_root,'tables');
globals.paths.resp_rp_path = fullfile(globals.paths.project_root,'resp_rp');
globals.paths.stim_rp_plotPath = fullfile(globals.paths.project_root,'stim_rp');
globals.paths.zebrisFigs = fullfile(globals.paths.project_root,'zebris','figures');
globals.paths.zebrisAscii = fullfile(globals.paths.project_root,'zebris','ascii');
globals.paths.mocapMatFiles = fullfile(globals.paths.project_root,'zebris','mat');
globals.paths.moviePath = fullfile(globals.paths.project_root,'animations');
globals.paths.mocapRP = fullfile(globals.paths.project_root,'move_rp');
globals.paths.highPassFilt = fullfile(globals.paths.project_root,'highpass_filt');

%% debug
globals.debug.mat_path = fullfile('/home','bkhurley','matlab','logs','polyMove_logs');
globals.debug.log_path = globals.debug.mat_path;
globals.debug.on = 1;

%% zebris preproc
% params required for preproc_zebris_v2
%globals.preProcLogFile = fullfile(globals.paths.log_path,'preproc.log');
% excluding 08cem92181 (movement not properly recorded)
globals.recordedDataFiles = {...
    '03ler90131','18256','18256.txt','stim_rec_03ler90131_18256.L.wav','03ler90131_polystream_move_18256_zebris.wav',[],[];...
    '06das92111','18264','18264.txt','stim_rec_06das92111_18264.L.wav','06das92111_polystream_move_18264_zebris.wav',[],[];...
    '07rzs93041','18280','18280.txt','stim_rec_07rzs93041_18280.L.wav','07rzs93041_polystream_move_18280_zebris.wav',[],[];...
    '05blr90291','18293','18293.txt','stim_rec_05blr90291_18293.L.wav','05blr90291_polystream_move_18293_zebris.wav',[],[];...
    '03mha83131','18316','18316.txt','stim_rec_03mha83131_18316.L.wav','03mha83131_polystream_move_18316_zebris.wav',[],[];...
    '11azl90131','18318','18318.txt','stim_rec_11azl90131_18318.L.wav','11azl90131_polystream_move_18318_zebris.wav',[],[];...
    '03fsj92061','18321','18321.txt','stim_rec_03fsj92061_18321.L.wav','03fsj92061_polystream_move_18321_zebris.wav',[],[];...
    '12vvd92141','18366','18366.txt','stim_rec_12vvd92141_18366.L.wav','12vvd92141_polystream_move_18366_zebris.wav',[],[];...
    '09doa90281','18399','18399.txt','stim_rec_09doa90281_18399.L.wav','09doa90281_polystream_move_18399_zebris.wav',[],[];...
    '12dov89231','18405','18405.txt','stim_rec_12dov89231_18405.L.wav','12dov89231_polystream_move_18405_zebris.wav',[],[];...
    '11nnd87261','18406','18406.txt','stim_rec_11nnd87261_18406.L.wav','11nnd87261_polystream_move_18406_zebris.wav',[],[];...
    '04moe93211','18407','18407.txt','stim_rec_04moe93211_18407.L.wav','04moe93211_polystream_move_18407_zebris.wav',[],[];...
    '12gas89241','18429','18429.txt','stim_rec_12gas89241_18429.L.wav','12gas89241_polystream_move_18429_zebris.wav',[],[];...
    '06kns93151','18468','18468.txt','stim_rec_06kns93151_18468.L.wav','06kns93151_polystream_move_18468_zebris.wav',[],[];...
    '09mlb93031','18498','18498.txt','stim_rec_09mlb93031_18498.L.wav','09mlb93031_polystream_move_18498_zebris.wav',[],[];...
    '09mnj90041','18521','18521.txt','stim_rec_09mnj90041_18521.L.wav','09mnj90041_polystream_move_18521_zebris.wav',[],[];...
    '09pnj92201','18523','18523.txt','stim_rec_09pnj92201_18523.L.wav','09pnj92201_polystream_move_18523_zebris.wav',[],[];...
    '04nka85251','18524','18524.txt','stim_rec_04nka85251_18524.L.wav','04nka85251_polystream_move_18524_zebris.wav',[],[];...
    
    % Data for subjects below below this point were collected using 3 markers
    % (vs. 2 markers as used for subjects above)
    '04dog91231','18813','18813.txt','stim_rec_04dog91231_18813.L.wav','04dog91231_polystream_move_18813_zebris.wav',[],[];...
    '08wyv86301','18814','18814.txt','stim_rec_08wyv86301_18814.L.wav','08wyv86301_polystream_move_18814_zebris.wav',[],[];...
    '07nam92101','18829','18829.txt','stim_rec_07nam92101_18829.L.wav','07nam92101_polystream_move_18829_zebris.wav',[],[];...
    '11eoj88181','18874','18874.txt','stim_rec_11eoj88181_18874.L.wav','11eoj88181_polystream_move_18874_zebris.wav',[],[];...
    '10frj91211','18884','18884.txt','stim_rec_10frj91211_18884.L.wav','10frj91211_polystream_move_18884_zebris.wav',[],[];...
    '06hgt92031','18886','18886.txt','stim_rec_06hgt92031_18886.L.wav','06hgt92031_polystream_move_18886_zebris.wav',[],[];...
    '07cuy92181','18892','18892.txt','stim_rec_07cuy92181_18892.L.wav','07cuy92181_polystream_move_18892_zebris.wav',[],[];...
    '11lzg89291','18893','18893.txt','stim_rec_11lzg89291_18893.L.wav','11lzg89291_polystream_move_18893_zebris.wav',[],[];...
    '07dnk91191','18917','18917.txt','stim_rec_07dnk91191_18917.L.wav','07dnk91191_polystream_move_18917_zebris.wav',[],[];...
    '01als91251','18919','18919.txt','stim_rec_01als91251_18919.L.wav','01als91251_polystream_move_18919_zebris.wav',[],[];...
    '11mrn85161','18934','18934.txt','stim_rec_11mrn85161_18934.L.wav','11mrn85161_polystream_move_18934_zebris.wav',[],[];...
    '02sij87191','18962','18962.txt','stim_rec_02sij87191_18962.L.wav','02sij87191_polystream_move_18962_zebris.wav',[],[];...
    '02dnm92151','18968','18968.txt','stim_rec_02dnm92151_18968.L.wav','02dnm92151_polystream_move_18968_zebris.wav',[],[];...
    '02lil92281','19004','19004.txt','stim_rec_02lil92281_19004.L.wav','02lil92281_polystream_move_19004_zebris.wav',[],[];...
    '01lua89011','19025','19025.txt','stim_rec_01lua89011_19025.L.wav','01lua89011_polystream_move_19025_zebris.wav',[],[];...
    '08cks89071','19036','19036.txt','stim_rec_08cks89071_19036.L.wav','08cks89071_polystream_move_19036_zebris.wav',[],[];...
    '09nos92301','19038','19038.txt','stim_rec_09nos92301_19038.L.wav','09nos92301_polystream_move_19038_zebris.wav',[],[];...
    '06wgw92241','19051','19051.txt','stim_rec_06wgw92241_19051.L.wav','06wgw92241_polystream_move_19051_zebris.wav',[],[];...
    '11gor92031','19052','19052.txt','stim_rec_11gor92031_19052.L.wav','11gor92031_polystream_move_19052_zebris.wav',[],[];...
    '11yej85252','19103','19103.txt','stim_rec_11yej85252_19103.L.wav','11yej85252_polystream_move_19103_zebris.wav',[],[];...
    '12rsm90151','19136','19136.txt','stim_rec_12rsm90151_19136.L.wav','12rsm90151_polystream_move_19136_zebris.wav',[],[];...
    % Movement only recorded for first 30 or so trials for this subject.
    % Excluding until I get around to adding code to mocap_matchWithAudio to handle shorter-than-expected subject files
    %'09kir90191','19208','19208.txt','stim_rec_09kir90191_19208.L.wav','09kir90191_polystream_move_19208_zebris.wav',[],[];...
    '10fnh93161','19209','19209.txt','stim_rec_10fnh93161_19209.L.wav','10fnh93161_polystream_move_19209_zebris.wav',[],[];...
    '06lra92281','19217','19217.txt','stim_rec_06lra92281_19217.L.wav','06lra92281_polystream_move_19217_zebris.wav',[],[];...
    '06czk92161','19255','19255.txt','stim_rec_06czk92161_19255.L.wav','06czk92161_polystream_move_19255_zebris.wav',[],[];...
    '09xgs84111','19268','19268.txt','stim_rec_09xgs84111_19268.L.wav','09xgs84111_polystream_move_19268_zebris.wav',[],[];...
    '01vao91202','19284','19284.txt','stim_rec_01vao91202_19284.L.wav','01vao91202_polystream_move_19284_zebris.wav',[],[];...
    '09fgl90261','19355','19355.txt','stim_rec_09fgl90261_19355.L.wav','09fgl90261_polystream_move_19355_zebris.wav',[],[];...
    '10jzc89111','19538','19538.txt','stim_rec_10jzc89111_19538.L.wav','10jzc89111_polystream_move_19538_zebris.wav',[],[];...
    '06mys88031','19374','19374.txt','stim_rec_06mys88031_19374.L.wav','06mys88031_polystream_move_19374_zebris.wav',[],[];...
    '01nnc91091','19393','19393.txt','stim_rec_01nnc91091_19393.L.wav','01nnc91091_polystream_move_19393_zebris.wav',[],[];...
    '06wan91101','19395','19395.txt','stim_rec_06wan91101_19395.L.wav','06wan91101_polystream_move_19395_zebris.wav',[],[];...
    '02bob91111','19396','19396.txt','stim_rec_02bob91111_19396.L.wav','02bob91111_polystream_move_19396_zebris.wav',[],[];...
    '04voc89021','19411','19411.txt','stim_rec_04voc89021_19411.L.wav','04voc89021_polystream_move_19411_zebris.wav',[],[];...
    '05lra91011','19413','19413.txt','stim_rec_05lra91011_19413.L.wav','05lra91011_polystream_move_19413_zebris.wav',[],[];...
    '04loa92251','19437','19437.txt','stim_rec_04loa92251_19437.L.wav','04loa92251_polystream_move_19437_zebris.wav',[],[];...
    '11ria90281','19438','19438.txt','stim_rec_11ria90281_19438.L.wav','11ria90281_polystream_move_19438_zebris.wav',[],[];...
    '11trk92171','19449','19449.txt','stim_rec_11trk92171_19449.L.wav','11trk92171_polystream_move_19449_zebris.wav',[],[];...
    '06lia92271','19471','19471.txt','stim_rec_06lia92271_19471.L.wav','06lia92271_polystream_move_19471_zebris.wav',[],[];...
    };           % ASCII files of raw Zebris data
globals.zebris.markerLabels = {'1','head';...
    '2','chest';...
    '3','nose'};
globals.zebris.sampleFreq = 100;                                                   % motion data were sampled @ 100 Hz
globals.zebris.ignoreSeg.shorterThanSec = 15;          
globals.ignoreSeg.indicesInFile = {};                                       % not sure what this should be??
globals.spline_interp.cut_off_nan = 0.25;
globals.spline_interp.signal_loss_info_fname = fullfile(globals.paths.tablepath,'mocap_sigloss_info.txt');
globals.manual_adjustment.reportStats = 'no';

% animation parameter struct
% can set connection matrix between markers here
globals.zebris.connList = {'head','chest';
                           'head','nose'};
globals.make_movie.mapar.showmnum = 1;
globals.make_movie.mapar.showfnum = 1;
globals.make_movie.mapar.colors = 'wbgcr';
globals.make_movie.mapar.msize = 6;
globals.make_movie.mapar.fps = 15;
globals.make_movie.outputPath = globals.paths.moviePath;

%for processing marker sigs through metric surface model
globals.ipem.glob.process = {'rp'};
globals.ipem.glob.force_recalc = {};
globals.ipem.glob.save_calc = {'rp'};
globals.ipem.glob.insigFs = globals.zebris.sampleFreq;
globals.ipem.glob.outputType = 'filepath';
globals.ipem.rp =  rp_paramGroups_v2('param_group','reson_filterQSpacing_periodBasedDecay',...
    'input_type','mat',... % input type must be set to 'mat', as calc_rp will average across dimenions by default for 'aud' input type
    'gain_type','beta_distribution');
% globals.ipem.rp.perform.calcOnsetInfo = 0;
globals.ipem.rp.perform.calcComplexOutput = 0;
% globals.ipem.rp.perform.interpMPP = 1;
% globals.ipem.rp.perform.calcVonMises = 1;
% globals.ipem.rp.perform.calcPhaseCoherence = 1;


%% Audio and MIDI parsing parameters
globals.parse_audio_stims.clip_pad_sec = 0;
globals.parse_audio_stims.thresh_silence_dB = -30;
globals.parse_audio_stims.detect_cues = 1;
globals.parse_audio_stims.throw_out_bad_cues = 1;
globals.parse_audio_stims.save_wvf = 1;
globals.parse_audio_stims.run_xcorr = 1;
globals.parse_audio_stims.xcorr_thresh = .60;
globals.parse_audio_stims.xcorr_fail_delete = 0;
globals.parse_audio_stims.xcorr_channels = 1;
globals.parse_audio_stims.pad_stim_wvf = 5;     % stims for this expt were padded by 5 sec @ beginning. won't pass xcorr unless I pad stim_wvf
globals.parse_audio_stims.mysql_conn_id = globals.mysql.conn_id;
globals.parse_audio_stims.path = globals.paths.stim_recording_path;
globals.parse_audio_stims.expect.min_clip_dur_sec = 5; % there were artifacts in some stim_rec files. This should take care of them
globals.parse_audio_stims.filename = 'stim_rec_<subject_id>_<session_id>.L.wav';

globals.parse_midi_resps.filename = 'midi_resp_<subject_id>_<session_id>.mid';
globals.parse_midi_resps.path = globals.paths.resp_midi_path;

globals.resamp.sampling_rate = 100;
globals.resamp.appendResampMidiResp = 1;

%% analysis

globals.verbose = 2;

% practice trials
globals.num_practice_trials_database_only = 2; % removes practice trials from ensemble session data without removing from audio recording (PTs not recorded in DP)

% filter criteria

globals.filt.exclude.any.subject_id = {'^tmp_.*','^01ttf.*','^01fmt.*',...
    '^03ttf.*','^04ttf.*','02weh90191','11atp79171','02mad87231',...
    '12hht90011','09jap67021','08huc93271','07rir90101','12nrl83311'};

globals.ensemble.filt = globals.filt;

globals.midi.validMidiNotes = [60 74];
globals.midi.leftHandMidiNote = 60;
globals.midi.rightHandMidiNote = 74;
globals.midi.ioiThresh = 0.050;

% tapping stats
globals.tapping_stats.midi = globals.midi;
globals.tapping_stats.normalizeByPeakIOI = 1;
globals.tapping_stats.ioiMinThresh = 0.050;
globals.tapping_stats.ioiHistBinLim = [0 3];
globals.tapping_stats.ioiHistBinSize = 0.05;
globals.tapping_stats.stdev.nBins = 6;
globals.tapping_stats.stdev.range = [0.75 1.25];

% stim RP params
globals.stimIPEMParams.ensemble = globals.ensemble;
globals.stimIPEMParams.paths.stimulus_root = '/nfs/stimuli';
globals.stimIPEMParams.paths.destroot = '/data2/stimuli';
globals.stimIPEMParams.paths.stimulus_ipem_analysis_root = '/data2/stimuli';
globals.stimIPEMParams.glob.process = {'ani','rp'};
globals.stimIPEMParams.glob.force_recalc = {};
globals.stimIPEMParams.glob.save_calc = {'ani','rp'};
globals.stimIPEMParams.glob.outputType = 'filepath';
globals.stimIPEMParams.ani = ani_paramGroups;
globals.stimIPEMParams.rp =  rp_paramGroups_v2('input_type','ani',...
    'gain_type','beta_distribution',...
    'param_group','reson_filterQSpacing_periodBasedDecay');

% Normalizing denv so that denv energy doesn't
% affect correlations between groove ratings and timing measurements
globals.stimIPEMParams.rp.perform.normalize.denv = 0;

% set parameters to calculate complex reson output,
% interpolated MPP, and temporal expectancy calculations
% globals.stimIPEMParams.rp.perform.calcOnsetInfo = 1;
globals.stimIPEMParams.rp.perform.calcComplexOutput = 1;
% globals.stimIPEMParams.rp.perform.interpMPP = 1;
% globals.stimIPEMParams.rp.perform.calcVonMises = 1;
% globals.stimIPEMParams.rp.perform.calcPhaseCoherence = 1;
% globals.stimIPEMParams.rp.perform.calcMPPByFrame = 1;

% response RP params
globals.respIPEMParams.glob.process = {'rp'};
globals.respIPEMParams.glob.save_calc = {'rp'};
globals.respIPEMParams.glob.force_recalc = {};
globals.respIPEMParams.glob.outputType = 'filepath';
globals.respIPEMParams.rp = rp_paramGroups_v2('param_group',...
				   'reson_filterQSpacing_periodBasedDecay',...
				   'input_type','aud','gain_type','beta_distribution');               
% globals.respIPEMParams.rp.perform.calcOnsetInfo = 1;
globals.respIPEMParams.rp.perform.calcComplexOutput = 0;
% globals.respIPEMParams.rp.perform.interpMPP = 1;
% globals.respIPEMParams.rp.perform.calcVonMises = 1;
% globals.respIPEMParams.rp.perform.calcPhaseCoherence = 1;

%% attribute grouping

% get stim names for each attribute group
for k=1:length(globals.attrib_names)                                   
    [stim_id,] = mysql_get_stim_by_attribute('attrib_name', ...
        globals.attrib_names{k}, 'mysql', globals.mysql);                      
    stimName_att{k} = stim_id{2};                            
end

%sequence style
globals.songatt.seq.sim = stimName_att{1};
globals.songatt.seq.stagStrong = stimName_att{2};
globals.songatt.seq.stagWeak = stimName_att{3};

%exemplar
globals.songatt.exc.dance1 = stimName_att{4};
globals.songatt.exc.dance2 = stimName_att{5};
globals.songatt.exc.dance3 = stimName_att{6};
globals.songatt.exc.dance4 = stimName_att{7};
globals.songatt.exc.dance5 = stimName_att{8};

globals.songatt.exc.folk1 = stimName_att{9};
globals.songatt.exc.folk2 = stimName_att{10};
globals.songatt.exc.folk3 = stimName_att{11};
globals.songatt.exc.folk4 = stimName_att{12};
globals.songatt.exc.folk5 = stimName_att{13};

globals.songatt.exc.funk1 = stimName_att{14};
globals.songatt.exc.funk2 = stimName_att{15};
globals.songatt.exc.funk3 = stimName_att{16};
globals.songatt.exc.funk4 = stimName_att{17};
globals.songatt.exc.funk5 = stimName_att{18};

globals.songatt.exc.rock1 = stimName_att{19};
globals.songatt.exc.rock2 = stimName_att{20};
globals.songatt.exc.rock3 = stimName_att{21};
globals.songatt.exc.rock4 = stimName_att{22};
globals.songatt.exc.rock5 = stimName_att{23};

globals.songatt.exc.chameleon = stimName_att{24};

%genre
globals.songatt.genre.dance = stimName_att{25};
globals.songatt.genre.folk = stimName_att{26};
globals.songatt.genre.funk = stimName_att{27};
globals.songatt.genre.rock = stimName_att{28};