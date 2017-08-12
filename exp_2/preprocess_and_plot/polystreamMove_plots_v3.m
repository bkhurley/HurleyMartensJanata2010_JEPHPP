function outData = polystreamMove_plots_v3(inData,params)

% plots for polystream_move experiment
%
% abstracted from polystreamMove_plots. As analyses develop,this should 
% probably be broken into several smaller functions that each plot different 
% aspects of the data 
%
% 04/28/2012 BH
% 04/29/2013 BH - cleaned up code and implemented plotting of tapping rate
%                 (need to code dynamic switching btw rate & total)
% 05/15/2013 BH - added dynamic switching between plotting tap rate and
%                 tapping stimulus-matching ratio (SMR) as DV

resp_st = inData;
resp_cols = set_var_col_const(resp_st.vars);

label_fontsize = 14;
axes_fontsize = 12;
info_text_fontweight = 'bold';
info_text_fontsize = 14;
title_fontsize = 16;
plotbox_aspectratio = [5 3 1];
marker_size = 8;
axis_extra = 0.5;

if params.collapse_stagConds
    legend_conds = {'simultaneous','staggered'};
    entrance_type = resp_st.data{resp_cols.entrance_type_stagCollapsed};
else
    legend_conds = {'simultaneous','staggered strong','staggered weak'};
    entrance_type = resp_st.data{resp_cols.entrance_type};
end

seq_conditions = unique(entrance_type);
seq_conditions(cellfun(@isempty,seq_conditions)) = [];
nconds = size(seq_conditions,1);
genres = fieldnames(params.songatt.genre);
ngenres = length(genres);
excerpts = fieldnames(params.songatt.exc);
nexcerpts = length(excerpts);
respStimNames = resp_st.data{resp_cols.stim_name};
stimNames = unique(respStimNames);
nstims = length(stimNames);
subids = unique(resp_st.data{resp_cols.subject_id});
nsubs = length(subids);
selectivity_level = resp_st.data{resp_cols.selectivity_level};
tapping_rate = resp_st.data{resp_cols.tapping_rate};
mean_entrainmentRatio = cell2mat(resp_st.data{resp_cols.tapping_meanMetricRatio});
postEntrance_tappingRate = resp_st.data{resp_cols.postEntrance_tappingRate};
postEntrance_entrainmentRatio = resp_st.data{resp_cols.tapping_postEntrance_meanMetricRatio};

% load table of entrance times
entrance_table = csvread('/data2/polystream_move/tables/polymove_entrance_times.csv');

atts = {'seq' 'genre' 'selectivity_level'}';
num_atts = length(atts);

levels = cell(1,num_atts);
% nlevels = nan(1,num_atts);

% enterConds = params.attrib_names(1:3);           % seq conditions
n_enterConds = length(seq_conditions);
selectivity_groups = unique(selectivity_level);       % groove level
n_selectivity_levels = length(selectivity_groups);
selectivity_titles = {'High Selectivity' 'Low Selectivity' 'Non-tappers'};
lowSelect_idx = ismember(selectivity_level,'low_selectivity');
highSelect_idx = ismember(selectivity_level,'high_selectivity');
nonTapper_idx = ismember(selectivity_level,'non_tapper');

musTrained_groups = unique(resp_st.data{resp_cols.musically_trained});

if params.plot_mean_tapRate
%     figname_bar_postEntrance = params.figfname_mean_tapRate;
    dep_var = tapping_rate;
    dep_var_name = 'tapping_rate';
    postEntrance_depVar = postEntrance_tappingRate;
    y_axis_label = 'Mean Taps/s';
%     if params.collapse_stagConds
%         figname_bar_postEntrance = fullfile(params.paths.fig_path,'polyMove_tapRate_stagCollapsed.ps');
%     else
%         figname_bar_postEntrance = fullfile(params.paths.fig_path,'polyMove_tapRate.ps');
%     end
elseif params.plot_entrainmentRatio
%     figname_bar_postEntrance = params.figfname_entrainmentRatio;
    dep_var = mean_entrainmentRatio;
    dep_var_name = 'tapping_meanMetricRatio';
    postEntrance_depVar = postEntrance_entrainmentRatio;
    y_axis_label = 'Mean SMR';
%     if params.collapse_stagConds
%         figname_bar_postEntrance = fullfile(params.paths.fig_path,'polyMove_tappingEntrainmentRatio_stagCollapsed.ps');
%     else
%         figname_bar_postEntrance = fullfile(params.paths.fig_path,'polyMove_tappingEntrainmentRatio.ps');
%     end
end

nfig = 0;
%% plot post-stim responses as a function of entrance type
postStim_resps = [resp_st.data{resp_cols.music_grooved} resp_st.data{resp_cols.enjoyed_music} ...
    resp_st.data{resp_cols.like_music_to_continue}];
postStim_label = {'Music Grooved','Enjoyed Music','Wanted to Continue'};
entrance_conditions = unique(resp_st.data{resp_cols.entrance_type_stagCollapsed});
entrance_conditions(~isempty(entrance_conditions)) = [];
nfig = nfig+1;
figure(nfig),clf
panel = 0;
resp_conds = entrance_type;
x_label = 'Entrance Type';
xtick_label = {'Simult.','Stag.'};  
for iresp = 1:length(postStim_label)
    panel = panel+1;
    postStim_bar_mean = nan(length(entrance_conditions),1);
    postStim_bar_std = nan(length(entrance_conditions),1);
    postStim_bar_nObs = nan(length(entrance_conditions),1);
    for icond = 1:length(entrance_conditions)
        curr_cond = entrance_conditions{icond};
        cond_mask = ismember(resp_conds,curr_cond);
        postStim_bar_mean(icond) = nanmean(postStim_resps(cond_mask,iresp));
        postStim_bar_std(icond) = nanstd(postStim_resps(cond_mask,iresp));
        postStim_bar_nObs(icond) = sum(~isnan(postStim_resps(cond_mask,iresp)));
    end
    subplot(2,3,panel), bh = bar(postStim_bar_mean);
    axis square
    set(gca,'ylim',[1 7]);
    colormap(gray);
    set(gca,'xticklabel',xtick_label) % ,'fontsize',axes_fontsize);
    xlabel(x_label,'fontsize',12);
    text(1.5,6.5,postStim_label{iresp},'fontsize',12,'HorizontalAlignment','center')
%     if iresp == 1
        ylabel('Mean Rating','fontsize',12);
%     end
    add_errorbars(bh, postStim_bar_std./sqrt(postStim_bar_nObs),'k');
    
    
end
print(fullfile(params.paths.fig_path,'polyMove_barMeans_postStim.eps'),'-depsc')

%% tapping distribution 
if params.tapping_dist
    hist(resp_st.data{respcols.total_tapping});
    print(params.figfname_tappingDist,'-dpsc','-adobecset')
end

%% Timecourse images by subject

% Plots matrix images of each subject's timecourse for staggered and
% simultaneous versions of each excerpt.

if params.tcourse_images

    fprintf('Printing figures: %s\n', figfname_tcourseImg);
    
    for iexcpt = 1:nexcerpts
        nfig = nfig+1;
        figure(nfig), clf
        set(gcf,'defaultaxesfontsize',axes_fontsize);
        curr_iexcpt = excerpts{iexcpt};
        hold on
        
        for iseq = 1:nseq
            subplot(3,1,(iseq))sfsfsfsf
            curr_iseq = levels{1}{iseq};
            imagesc(resp_st.data{resp_cols.tcourseMtx_st}.(curr_iexcpt).(curr_iseq));
            old_x_scale = get(gca,'xtick');
            new_x_scale = old_x_scale/100;
            set(gca,'xticklabel',new_x_scale);
            xlabel('time (sec)','fontsize',label_fontsize);
            ylabel('subjects','fontsize',label_fontsize);
            title(sprintf('%s %s',curr_iexcpt,curr_iseq));
            colorbar
        end
        
        switch iexcpt
            case 1
                print(params.figfname_tcourseImg,'-dpsc','-adobecset')
            otherwise
                print(params.figfname_tcourseImg,'-dpsc','-adobecset','-append')
        end
    end % iexcpt
end

%% Mean timecourse
% plots means of each timecourse point for each staggered & simultaneous
% conditions of each exemplar.

if params.plot_timecourse

    fprintf('Printing figures: %s\n', params.figfname_tcoursemeans);
    
    for iexcpt = 1:nexcerpts-1
        nfig = nfig+1;
        figure(nfig), clf
        set(gcf,'defaultaxesfontsize',axes_fontsize)
        tcourseData = resp_st.data{resp_cols.tcourseMean_st}.(excerpts{iexcpt});
        entrTicks_strong = entrance_table(iexcpt*3,2:end);
        entrTicks_weak = entrance_table(iexcpt*3+1,2:end);
        switch entrTicks_strong(end)
            case 0
                entrTicks_strong(end) = []; % delete 5th element for 4-entrance stims
                entrTicks_weak(end) = [];
                enterLables = {'e1','e2','e3','e4'};
            otherwise
                enterLables = {'e1','e2','e3','e4','e5'};
        end
        hold on
        
        sim_means = tcourseData{1};
        e1 = (std(sim_means)/sqrt(nsubs))*ones(length(tcourseData{1}),1)';
        simUpperErr = sim_means + e1/2;
        simLowerErr = sim_means - e1/2;
        
        stagStrong_means = tcourseData{2};
        e2 = (std(stagStrong_means)./sqrt(nsubs))*ones(length(tcourseData{2}),1)';
        stagStrUpperErr = stagStrong_means + e2/2;
        stagStrLowerErr = stagStrong_means - e2/2;
        
        stagWk_means = tcourseData{3};
        e3 = (std(stagWk_means)./sqrt(nsubs))*ones(length(tcourseData{3}),1)';
        stagWkUpperErr = stagWk_means + e3/2;
        stagWkLowerErr = stagWk_means - e3/2;
        
        sim_fill = jbfill(1:length(sim_means),simUpperErr,simLowerErr,'b','b',1,1);
        stgStr_fill = jbfill(1:length(stagStrong_means),stagStrUpperErr,stagStrLowerErr,'r','r',1,1);
        stgWk_fill = jbfill(1:length(stagWk_means),stagWkUpperErr,stagWkLowerErr,'g','g',1,1);
        
        
        
        data1 = line((1:length(sim_means)),sim_means,'Color','k');
        data2 = line((1:length(stagStrong_means)),stagStrong_means,'Color','k');
        data3 = line((1:length(stagWk_means)),stagWk_means,'Color','k');
        
        ylabel('mean number of taps','fontsize',label_fontsize);
        set(gca,'XMinorTick','on');
        %    set(gca,'xtick',1:200:length(imgData),'xlim',[1 length(imgData)]);
        yl = get(gca,'Ylim');
        for iline = 1:length(entrTicks_strong)
            line([entrTicks_strong(iline) entrTicks_strong(iline)],yl,'LineWidth',1,'Color','black');
            line([entrTicks_weak(iline) entrTicks_weak(iline)],yl,'LineWidth',1,'LineStyle','--','Color','black');
        end
        old_x_scale = get(gca,'xtick');
        new_x_scale = old_x_scale/100;
        set(gca,'xticklabel',new_x_scale);
        xlabel('time (sec)','fontsize',label_fontsize);
        title(sprintf('%s',excerpts{iexcpt}),'fontsize',16,...
            'fontweight',info_text_fontweight);
        axis square
        legend([sim_fill,stgStr_fill,stgWk_fill],'simultaneous','staggered strong','staggered weak','Location','Best');
        
        switch iexcpt
            case 1
                print(params.figfname_tcoursemeans,'-dpsc','-adobecset')
            otherwise
                print(params.figfname_tcoursemeans,'-dpsc','-adobecset','-append')
        end
        hold off
    end
end

%% Mean tapping rate plots
if params.plot_mean_tapRate || params.plot_entrainmentRatio
    meanData_mean = nan(n_selectivity_levels,n_enterConds);
    meanData_stdDev = nan(n_selectivity_levels,n_enterConds);
    meanData_nObs = nan(n_selectivity_levels,n_enterConds);
    
    chameleon_postEntrance_tappingRate = nan(5,n_selectivity_levels);
    chameleon_postEntrance_stderr = nan(5,n_selectivity_levels);

    for i = 1:n_selectivity_levels
        curr_selLevel = selectivity_groups{i};
        selectivity_mask = ismember(selectivity_level,curr_selLevel);
        
        % piggy back on the first 2 selectivity loop iterations for musicality
        % levels
        if i <= 2
           curr_musTrn_grp = musTrained_groups{i}; 
           musTrain_msk = ismember(resp_st.data{resp_cols.musically_trained},curr_musTrn_grp);
        end
            
        for iseq = 1:n_enterConds
            curr_seq = seq_conditions{iseq};
            seq_mask = ismember(entrance_type,curr_seq);
            seqCond1_mask = selectivity_mask&seq_mask;
            curr_meanData = dep_var(seqCond1_mask);
            meanData_mean(i,iseq) = nanmean(curr_meanData);
            
            meanData_stdDev(i,iseq) = nanstd(curr_meanData);
            meanData_nObs(i,iseq) = sum(~isnan(curr_meanData));
            
            curr_postEntrance_plotData(iseq,:) = nanmean(cell2mat(postEntrance_depVar(seqCond1_mask)));
            cur_nObs_vect = meanData_nObs(i,iseq)*ones(1,size(curr_postEntrance_plotData,2));
            curr_postEntrance_stderr(iseq,:) = nanstd(cell2mat(postEntrance_depVar(seqCond1_mask)))./(sqrt(cur_nObs_vect));
            
            % prepare musicality plot data
            if i <= 2
                entrtypeTraining_msk = seq_mask&musTrain_msk;
                curr_trialwise_data = dep_var(entrtypeTraining_msk); % trialwise data under the current crossed conditions
                this_nObs = sum(~isnan(curr_trialwise_data));
                switch i
                    case 1
                        idx = iseq;
                    case 2
                        idx = iseq + i;
                end
                % rows of data matrices: (1) nontrained_sim (2) nontrained_stag (3) trained_sim (4) trained_stag
                % columns correspond to entrance periods
                postEntrance_musTrainData(idx,:) = ...
                    nanmean(cell2mat(postEntrance_depVar(entrtypeTraining_msk)));
                musTrain_nObs = this_nObs*ones(1,size(postEntrance_musTrainData,2));
                postEntrance_musTrain_SE(idx,:) = ...
                    nanstd(cell2mat(postEntrance_depVar(entrtypeTraining_msk)))./(sqrt(musTrain_nObs));
            end
        end
        
        if nfig == 1 || params.plot_entrainmentRatio
            nfig = nfig+1;
            figure(nfig),clf
            hold on
            set(gcf,'colormap',gray)
        end
        y = curr_postEntrance_plotData';
        maxval = max(max(y));
        if strcmp(curr_selLevel,'non_tapper')
            maxval = maxval+.2;
        end
        e = curr_postEntrance_stderr';
        lspec = {'k-' 'k--' 'k:'};
        markers = {'o','s','d'};
        for iCond = 1:nconds
            curr_cond_data = y(:,iCond);
            curr_cond_err = e(:,iCond);
            eh = errorbar(curr_cond_data,curr_cond_err,lspec{iCond},...
                'marker',markers{iCond},'markerfacecolor','k','linewidth',1.01);
        end
        set(gca,'xtick',1:4,'xticklabel',{'Baseline' '2' '3' '4'},'fontsize',axes_fontsize);
        curr_ylim = get(gca,'Ylim');
%         ylim_range = curr_ylim(2) - curr_ylim(1);
%         set(gca,'ytick',curr_ylim(1):(ylim_range/8):curr_ylim(2));
        ylims{i} = curr_ylim;
%         set(gca,'Ylim',[.95 1.35]);
        xlabel(sprintf('Entrance\n'),'fontsize',label_fontsize);
        ylabel(y_axis_label,'fontsize',label_fontsize);
        
        if params.plot_mean_tapRate
            text(2.5,maxval,sprintf('%s',selectivity_titles{i}),...
                'horizontalalign','center', ...
                'fontsize', info_text_fontsize, ...
                'fontweight', info_text_fontweight);
            if i == 3
                set(gca,'Ylim',[0 4])
                legend(legend_conds,'Location','SouthOutside','fontsize',14),legend 'boxoff';
                set(gcf, 'PaperPosition', [1 .5 6 10]);
                print(fullfile(params.paths.fig_path,sprintf('polyMove_postEntr_%s.eps',dep_var_name)),'-depsc')
            end
            
        elseif params.plot_entrainmentRatio
            set(gca,'Ylim',[.95 1.35]);
            text(2.5,1.34,sprintf('%s',selectivity_titles{i}),...
                'horizontalalign','center', ...
                'fontsize', 16, ...
                'fontweight', info_text_fontweight);
            
            
            print(fullfile(params.paths.fig_path,...
                sprintf('polyMove_postEntr_%s_%s.eps',dep_var_name,curr_selLevel)),'-depsc')
        end
    
        % prepare chameleon data for plotting
        chameleon_mask = ismember(respStimNames,'chameleon.mp3');
        selectlvl_chameleon_mask = selectivity_mask&chameleon_mask;
        chameleon_postEntrance_tappingRate(:,i) = nanmean(cell2mat(postEntrance_tappingRate(selectlvl_chameleon_mask)));
        chameleon_curr_nObs = sum(~isnan(cell2mat(postEntrance_tappingRate(selectlvl_chameleon_mask))));
        chameleon_postEntrance_stderr(:,i) = nanstd(cell2mat(postEntrance_tappingRate(selectlvl_chameleon_mask)))./(sqrt(chameleon_curr_nObs));
        
    end
    
    mean_data = meanData_mean;
    stdev = meanData_stdDev;
    num_obs = meanData_nObs;
    plot_legend = 1;
    plot_colors = gray;
    
    nfig = nfig+1;
    figure(nfig),clf
    bh = bar(mean_data);
    colormap(plot_colors);
    set(gca,'xtick',1:n_selectivity_levels,'xticklabel',selectivity_titles);
    xlabel('Tapping Selectivity','fontsize',label_fontsize);
    ylabel(y_axis_label,'fontsize',label_fontsize);
    if plot_legend
        legend(legend_conds), legend boxoff
    end
    add_errorbars(bh, stdev./sqrt(num_obs),'k');
    
    if params.plot_entrainmentRatio
        ylims_mat = cell2mat(ylims);
        bar_ylim = [min(ylims_mat) max(ylims_mat)];
        set(gca,'Ylim',bar_ylim);
    end
    
    % print to file
    set(gcf, 'PaperPosition', [0.25 2.5 8.0 6.0]);
    print(fullfile(params.paths.fig_path,sprintf('polyMove_%s_barMeans.eps',dep_var_name)),'-depsc')
    
    % chameleon post-entrance plot
    
    nfig = nfig+1;
    figure(nfig),clf
    hold on
    set(gcf,'colormap',gray)
    y = chameleon_postEntrance_tappingRate;
    e = chameleon_postEntrance_stderr;
    lspec = {'k-' 'k--' 'k:'};
    markers = {'o','s','d'};
    for isel = 1:n_selectivity_levels
        curr_cond_data = y(:,isel);
        curr_cond_err = e(:,isel);
        eh = errorbar(curr_cond_data,curr_cond_err,lspec{isel},'marker',markers{isel},...
            'markerfacecolor','k','linewidth',1.01);
    end
    set(gca,'xtick',1:5,'xticklabel',{'pre-entrance2' 'entrance 2' 'entrance 3' 'entrance 4' 'entrance 5'},...
        'fontsize',axes_fontsize);
    legend(selectivity_titles,'Location','NorthWest'), legend 'boxoff';
    xlabel('Entrance','fontsize',label_fontsize);
    ylabel('Mean Taps/s','fontsize',label_fontsize);
    title('Stimulus: Chameleon by Herbie Hancock','fontsize',16,...
        'fontweight',info_text_fontweight);
    print(fullfile(params.paths.fig_path,sprintf('polyMove_postEntr_chameleon_%s.eps',dep_var_name)),'-depsc')
    
    % musicality post-entrance plot
    nfig = nfig+1;
    figure(nfig),clf
    hold on
    set(gcf,'colormap',gray)
    musTrainPlot_data = postEntrance_musTrainData';
    y_range = [min(min(musTrainPlot_data)) max(max(musTrainPlot_data))]; 
    musTrainPlot_SE = postEntrance_musTrain_SE';
    musTrain_lspec = {'-.ko','--k^','-ks','k:d'};
    marker_face_color = {'w','w','k','k'};
    for musGrp = 1:size(musTrainPlot_data,2)
        curr_musGrp_data = musTrainPlot_data(:,musGrp);
        curr_musGrp_SE = musTrainPlot_SE(:,musGrp);
        eb_h = errorbar(curr_musGrp_data,curr_musGrp_SE,musTrain_lspec{musGrp},'linewidth',1.01,...
            'MarkerFaceColor',marker_face_color{musGrp},'MarkerSize',8);
    end
    if strcmp(dep_var_name,'tapping_rate')
        set(gca,'Ylim',[0 y_range(2)*1.10])
    else
        set(gca,'Ylim',[1 1.35])
        h_leg = legend({'Simultaneous, Non-Trained','Staggered, Non-Trained','Simultaneous, Trained','Staggered, Trained'},...
            'Location','SouthEast'); legend 'boxoff'
        set(h_leg,'FontSize',label_fontsize)
    end
    set(gca,'xtick',1:4,'xticklabel',{'Baseline' '2' '3' '4'},'fontsize',axes_fontsize)
    xlabel(sprintf('Entrance\n'),'fontsize',label_fontsize);
    ylabel(y_axis_label,'fontsize',label_fontsize);
    axis square
    
    % print to file
    print(fullfile(params.paths.fig_path,sprintf('polyMove_postEntr_musTrain_%s.eps',dep_var_name)),'-depsc')
end


%% Correlations between post-stim questions

if params.plot_poststimCorr
    
    % Scatter plots of post-stim ratings by tap response (either tap rate or
    % entrainment)
    if params.plot_mean_tapRate
        tapscatter_figfname = params.figfname_tapRate_poststim_scatter;
    elseif params.plot_entrainmentRatio
        tapscatter_figfname = params.figfname_tappingEntrainment_poststim_scatter;
    end
    resp_vars = [{dep_var_name} {'music_grooved'} {'enjoyed_music'} {'like_music_to_continue'}];
    nVars = length(resp_vars);
    entrance_conds = unique(resp_st.data{resp_cols.entrance_type});
    % remove empty cells (which are associated with chameleon, which has no
    % entrance condition)
    entrance_conds(cellfun(@isempty,entrance_conds)) = [];
    nEntrance_conds = length(entrance_conds);
    stimResp_stats.nanmean = nan(nstims,nVars);
    stimResp_stats.nanstd = nan(nstims,nVars);
    stimResp_stats.numObs = nan(nstims,nVars);
    poststim_stats = fieldnames(stimResp_stats);
    
    % get descriptive stats on subjective ratings for each stimulus
    for istim = 1:nstims
        currStim = stimNames{istim};
        stimResp_mask = ismember(resp_st.data{resp_cols.stim_name},currStim);
        for iVar = 1:nVars
            curr_resp = resp_st.data{resp_cols.(resp_vars{iVar})}(stimResp_mask);
            if iscell(curr_resp)
                curr_resp = cell2mat(curr_resp);
            end
            for istat = 1:length(poststim_stats)
                curr_stat = poststim_stats{istat};
                switch curr_stat
                    case 'numObs'
                        stimResp_stats.(curr_stat)(istim,iVar) = sum(~isnan(curr_resp));
                    otherwise
                        fh = str2func(curr_stat);
                        stimResp_stats.(curr_stat)(istim,iVar) = fh(curr_resp);
                end
            end
        end
    end
    
    for ivar = 2:nVars
        nfig = nfig+1;
        figure(nfig),clf;
        set(gcf,'defaultaxesfontsize',axes_fontsize)
        curr_var = resp_vars{ivar};
        switch curr_var
            case 'music_grooved'
                x_label = 'Music Grooved';
                x_tickLabel = {'least groove','','','','','','most groove'};
            case 'enjoyed_music'
                x_label = 'Enjoyed Music';
                x_tickLabel = {'not at all','','','','','','very much'};
            case 'like_music_to_continue'
                x_label = 'Wanted Music to Continue';
                x_tickLabel = {'not at all','','','','','','very much'};
        end
        
        x_data = stimResp_stats.nanmean(:,ivar)';
        y_data = stimResp_stats.nanmean(:,1)'; % first column is entrainment
        
        p = plot(x_data,y_data,'o','markersize', marker_size,'markeredgecolor','k');
        lh = lsline;
        set(lh,'color','k','linewidth',1);
        
        [rho,pval] = corr(x_data',y_data');
        text(1,2.4,sprintf('r=%1.2f, p<%1.4f', rho, max(0.0001,pval)), ...
            'horizontalalign','left', ...
            'fontsize', info_text_fontsize, ...
            'fontweight', info_text_fontweight);
        text(7,.9,sprintf('N=%d stimuli', size(stimResp_stats.numObs,1)), ...
            'horizontalalign','right', ...
            'fontsize', info_text_fontsize, ...
            'fontweight', info_text_fontweight);
        
        
        xlabel(x_label,'fontsize',label_fontsize)
        ylabel(y_axis_label,'fontsize',label_fontsize)
        x_enum_values = resp_st.meta.question(2).enum_values;
        
        axis_extra = 0.5;
        
        set(gca,'xtick', 1:length(x_enum_values),'xlim',[1-axis_extra length(x_enum_values)+axis_extra])
        set(gca,'xticklabel',x_tickLabel);
        
        axis square
        
        switch ivar
            case 2
                print(tapscatter_figfname,'-dpsc','-adobecset');
            otherwise
                print(tapscatter_figfname,'-dpsc','-adobecset','-append');
        end
        
    end
    
    
    
    %
    % Music Grooved vs Enjoyment
    %
    nfig = nfig+1;
    figure(nfig), clf
    set(gcf,'defaultaxesfontsize',axes_fontsize)
    xdata = stim_att_mtx.nanmean(:,2);
    ydata = stim_att_mtx.nanmean(:,1);
    p = plot(xdata,ydata,'o','markersize', marker_size);
    lh = lsline;
    set(lh,'color',[1 0 0]); % make it red
    
    [rho,pval] = corr(xdata,ydata);
    text(1,7,sprintf('r=%1.2f, p<%1.4f', rho, max(0.0001,pval)), ...
        'horizontalalign','left', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    text(7,1,sprintf('N = %d excerpts', length(xdata)), ...
        'horizontalalign','right', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    
    xlabel('Mean Enjoyment','fontsize',label_fontsize)
    ylabel('Mean Groove in Music','fontsize',label_fontsize)
    x_enum_values = resp_st.meta.question(2).enum_values;
    y_enum_values = resp_st.meta.question(1).enum_values;
    
    axis_extra = 0.5;
    
    set(gca,'xtick', 1:length(x_enum_values),'xlim',[1-axis_extra length(x_enum_values)+axis_extra])
    set(gca,'xticklabel',{'not at all','','','','','','very much'});
    
    set(gca,'ytick', 1:length(y_enum_values),'ylim',[1-axis_extra length(y_enum_values)+axis_extra])
    set(gca,'yticklabel',{'least groove','','','','','','most groove'});
    axis square
    
    print(params.figfname,'-dpsc','-adobecset');
    
    %
    % Urge to move vs want to continue
    %
    nfig = nfig+1;
    figure(nfig), clf
    set(gcf,'defaultaxesfontsize',axes_fontsize)
    xdata = stim_att_mtx.nanmean(:,3);
    ydata = stim_att_mtx.nanmean(:,1);
    p = plot(xdata,ydata,'o','markersize', marker_size);
    lh = lsline;
    set(lh,'color',[1 0 0]); % make it red
    
    [rho,pval] = corr(xdata,ydata);
    text(1,7,sprintf('r=%1.2f, p<%1.4f', rho, max(0.0001,pval)), ...
        'horizontalalign','left', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    text(7,1,sprintf('N = %d excerpts', length(xdata)), ...
        'horizontalalign','right', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    
    xlabel('Mean Like to Continue','fontsize',label_fontsize)
    ylabel('Mean Groove in Music','fontsize',label_fontsize)
    x_enum_values = resp_st.meta.question(3).enum_values;
    y_enum_values = resp_st.meta.question(1).enum_values;
    
    axis_extra = 0.5;
    
    set(gca,'xtick', 1:length(x_enum_values),'xlim',[1-axis_extra length(x_enum_values)+axis_extra])
    set(gca,'xticklabel',{'not at all','','','','','','very much'});
    
    set(gca,'ytick', 1:length(y_enum_values),'ylim',[1-axis_extra length(y_enum_values)+axis_extra])
    set(gca,'yticklabel',{'least groove','','','','','','most groove'});
    axis square
    
    print(params.figfname,'-dpsc','-adobecset','-append');
    
    %
    % Enjoyment vs. Want to Continue
    %
    nfig = nfig+1;
    figure(nfig), clf
    set(gcf,'defaultaxesfontsize',axes_fontsize)
    xdata = stim_att_mtx.nanmean(:,2);
    ydata = stim_att_mtx.nanmean(:,3);
    p = plot(xdata,ydata,'o','markersize', marker_size);
    lh = lsline;
    set(lh,'color',[1 0 0]); % make it red
    
    [rho,pval] = corr(xdata,ydata);
    text(1,7,sprintf('r=%1.2f, p<%1.4f', rho, max(0.0001,pval)), ...
        'horizontalalign','left', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    text(7,1,sprintf('N = %d excerpts', length(xdata)), ...
        'horizontalalign','right', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    
    xlabel('Mean Enjoyment','fontsize',label_fontsize)
    ylabel('Mean Like to Continue','fontsize',label_fontsize)
    x_enum_values = resp_st.meta.question(2).enum_values;
    y_enum_values = resp_st.meta.question(3).enum_values;
    
    set(gca,'xtick', 1:length(x_enum_values),'xlim',[1-axis_extra length(x_enum_values)+axis_extra])
    set(gca,'xticklabel',{'not at all','','','','','','very much'});
    
    set(gca,'ytick', 1:length(y_enum_values),'ylim',[1-axis_extra length(y_enum_values)+axis_extra])
    set(gca,'yticklabel',{'not at all','','','','','','very much'});
    axis square
    
    print(params.figfname,'-dpsc','-adobecset','-append');
    
    %% Correlations between post-stim questions & total taps
    
    %
    % Enjoyment vs. Number of taps
    %
    nfig = nfig+1;
    figure(nfig), clf
    set(gcf,'defaultaxesfontsize',axes_fontsize)
    xdata = stim_att_mtx.nanmean(:,2);
    ydata = stim_att_mtx.nanmean(:,4);
    p = plot(xdata,ydata,'o','markersize', marker_size, 'color', ones(1,3)*.03);
    lh = lsline;
    set(lh,'color',[0 0 0],'linewidth',2);
    
    [rho,pval] = corr(xdata,ydata);
    text(1,130,sprintf('r=%1.2f, p<%1.4f', rho, max(0.0001,pval)), ...
        'horizontalalign','left', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    text(7,10,sprintf('N = %d excerpts', length(xdata)), ...
        'horizontalalign','right', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    
    xlabel('Mean Enjoyment','fontsize',label_fontsize)
    ylabel(sprintf('Mean Number of Taps'),'fontsize',label_fontsize)
    x_enum_values = resp_st.meta.question(2).enum_values;
    
    set(gca,'xtick', 1:length(x_enum_values),'xlim',[1-axis_extra length(x_enum_values)+axis_extra])
    set(gca,'xticklabel',{'not at all','','','','','','very much'});
    maxTaps = max(stim_att_mtx.nanmean(:,4));
    set(gca,'ylim',[0 maxTaps+10]);
    axis square
    
    print(params.figfname,'-dpsc','-adobecset','-append');
    
    %
    % Music Grooved vs. Amount of Tapping
    %
    nfig = nfig+1;
    figure(nfig), clf
    set(gcf,'defaultaxesfontsize',axes_fontsize)
    xdata = stim_att_mtx.nanmean(:,1);
    ydata = stim_att_mtx.nanmean(:,4);
    p = plot(xdata,ydata,'o','markersize', marker_size, 'color', ones(1,3)*.03);
    lh = lsline;
    set(lh,'color',[0 0 0],'linewidth',2);
    
    [rho,pval] = corr(xdata,ydata);
    text(1,130,sprintf('r=%1.2f, p<%1.4f', rho, max(0.0001,pval)), ...
        'horizontalalign','left', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    text(7,10,sprintf('N = %d excerpts', length(xdata)), ...
        'horizontalalign','right', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    
    xlabel('Mean Groove in Music','fontsize',label_fontsize)
    ylabel(sprintf('Mean Number of Taps'),'fontsize',label_fontsize)
    x_enum_values = resp_st.meta.question(1).enum_values;
    
    set(gca,'xtick', 1:length(x_enum_values),'xlim',[1-axis_extra length(x_enum_values)+axis_extra])
    set(gca,'xticklabel',{'least groove','','','','','','most groove'});
    maxTaps = max(stim_att_mtx.nanmean(:,4));
    set(gca,'ylim',[0 maxTaps+10]);
    axis square
    
    print(params.figfname,'-dpsc','-adobecset','-append');
    
    %
    % Like to continue vs. Number of Taps
    %
    nfig = nfig+1;
    figure(nfig), clf
    set(gcf,'defaultaxesfontsize',axes_fontsize)
    xdata = stim_att_mtx.nanmean(:,3);
    ydata = stim_att_mtx.nanmean(:,4);
    p = plot(xdata,ydata,'o','markersize', marker_size, 'color', ones(1,3)*.03);
    lh = lsline;
    set(lh,'color',[0 0 0],'linewidth',2);
    
    [rho,pval] = corr(xdata,ydata);
    text(1,135,sprintf('r=%1.2f, p<%1.4f', rho, max(0.0001,pval)), ...
        'horizontalalign','left', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    text(7,10,sprintf('N = %d excerpts', length(xdata)), ...
        'horizontalalign','right', ...
        'fontsize', info_text_fontsize, ...
        'fontweight', info_text_fontweight);
    
    xlabel('Mean Like to Continue','fontsize',label_fontsize)
    ylabel(sprintf('Mean Number of Taps'),'fontsize',label_fontsize)
    x_enum_values = resp_st.meta.question(3).enum_values;
    
    set(gca,'xtick', 1:length(x_enum_values),'xlim',[1-axis_extra length(x_enum_values)+axis_extra])
    set(gca,'xticklabel',{'not at all','','','','','','very much'});
    maxTaps = max(stim_att_mtx.nanmean(:,4));
    set(gca,'ylim',[0 maxTaps+10]);
    axis square
    
    print(params.figfname,'-dpsc','-adobecset','-append');
    
    close all
    
end % if params.plot_poststimCorr

%% tapping distributions by stimulus

%
% Distributions of tapping by stimulus
%

if params.plot_trialStats
    clf
    stimids = trialStat_st.data{trialStat_cols.stimulus_id};
    nstims = length(stimids);
    stimNames = trialStat_st.data{trialStat_cols.stimulus_name};
    tapping = trialStat_st.data{trialStat_cols.data};
    fig_pos = 0;
    
    for istim = 1:nstims
        curr_name = cell2str(stimNames{istim});
        % remove '.mp3' from stimulus name
        curr_name = curr_name(1:end-4);
        
        if strcmp(curr_name,'chameleon')
            clf
            hist(tapping{istim})
            title(curr_name)
            print(params.figfname_trial_stats,'-dpsc','-adobecset','-append')
            xlim([0 400])
            set(gca,'xtick', 0:200:400)
        else
            fig_pos = fig_pos+1; %position of figure in subplot
            subplot(5,3,fig_pos); hist(tapping{istim})
            title(curr_name)
            xlim([0 400])
            set(gca,'xtick', 0:200:400)
            
            if fig_pos == 15               
                fig_pos = 0; %start new position sequence after every 3rd hist
                switch istim
                    case 15
                        print(params.figfname_trial_stats,'-dpsc','-adobecset')
                    otherwise
                        print(params.figfname_trial_stats,'-dpsc','-adobecset','-append')
                end
                % clear current figure for next group of genre plots
                clf
                
            end % if fig_pos
        
        end % if strcmp
        
    end % for istim =
    
end % if params.plot_trialStats

%% Distributions of poststim-tapping correlations for each poststim question

if params.plot_poststimTapCorrDists
    
    nfig = 0;
    poststim_vars = poststimTap_st.vars;
    for ivar = 2:length(poststim_vars)
        switch ivar
            case 2
                curr_title = 'Music Grooved';
            case 3
                curr_title = 'Enjoyed Music';
            case 4
                curr_title = 'Want to Continue';
        end
        
        % distribution of subject poststim-tapping correlations
        nfig = nfig+1;
        figure(nfig), clf;
        set(gcf,'defaultaxesfontsize',axes_fontsize);
        curr_data = poststimTap_st.data{ivar};
        hist(curr_data);
        
        title(curr_title,'fontsize',title_fontsize);
        xlabel('correlation','fontsize',label_fontsize);
        ylabel('#subjects','fontsize',label_fontsize);
        axis square
                
        if nfig == 1
            print(params.figfname_poststimTapCorrDists,'-dpsc','-adobecset');
        else
            print(params.figfname_poststimTapCorrDists,'-dpsc','-adobecset','-append');
        end
        
    end % for ivar
    
end % if params.plot_poststimTapCorrDists

close all
outData = [];