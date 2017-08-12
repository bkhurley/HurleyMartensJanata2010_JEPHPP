function result_st = poststim_tapping_analysis(inData,params)

%
% Generates a table of correlations between tapping and poststim questions for each
% subject. May identify (1) subjects that do not tap regardless of how much they
% like stimuli, and (2) subjects that tap regardless of how much/how little
% they like stimuli.
% 
% 5-8-12 BH
%

tapResp_st = inData;
tapRespCols = set_var_col_const(tapResp_st.vars);

% response variables
resp_vars = {'tap_musicGroove_corr','tap_enjoy_corr','tap_continue_corr','tap_compositeQuest_corr'};
nvars = length(resp_vars);

result_st = ensemble_init_data_struct;
result_st.name = 'poststimTap_corr';

subid_vect = tapResp_st.data{tapRespCols.subject_id};
subids = unique(subid_vect);
nsub = length(subids);

% initialize vects
meanTapVect = zeros(nsub,1);
meanGrooveVect = zeros(nsub,1);
meanEnjoyVect = zeros(nsub,1);
meanContinueVect = zeros(nsub,1);
tapGroove_corr = zeros(nsub,1);
tapEnjoy_corr = zeros(nsub,1);
tapCont_corr = zeros(nsub,1);

% calculate correlations
for isub = 1:nsub
    curr_sub = subids{isub};
    sub_idx = strmatch(curr_sub,subid_vect);
    curr_tapVect = tapResp_st.data{tapRespCols.total_tapping}(sub_idx);
    curr_grooveVect = tapResp_st.data{tapRespCols.music_grooved}(sub_idx);
    curr_enjoyVect = tapResp_st.data{tapRespCols.enjoyed_music}(sub_idx);
    curr_continueVect = tapResp_st.data{tapRespCols.like_music_to_continue}(sub_idx);
    tapGroove_corr(isub,1) = corr(curr_tapVect,curr_grooveVect);
    tapEnjoy_corr(isub,1) = corr(curr_tapVect,curr_enjoyVect);
    tapCont_corr(isub,1) = corr(curr_tapVect,curr_continueVect);
    % correlation between mean tapping & mean ratings by subject
    meanTapVect(isub) = mean(tapResp_st.data{tapRespCols.total_tapping}(sub_idx));
    meanGrooveVect(isub) = mean(tapResp_st.data{tapRespCols.music_grooved}(sub_idx));
    meanEnjoyVect(isub) = mean(tapResp_st.data{tapRespCols.enjoyed_music}(sub_idx));
    meanContinueVect(isub) = mean(tapResp_st.data{tapRespCols.like_music_to_continue}(sub_idx));
    
end

resp_mtx{1} = subids;
resp_mtx{2} = tapGroove_corr;
resp_mtx{3} = tapEnjoy_corr;
resp_mtx{4} = tapCont_corr;

%% plot correlations
nfig = 0;
if params.plot_corr
    
    % set figure parameters
    label_fontsize = 14;
    axes_fontsize = 12;
    info_text_fontweight = 'bold';
    info_text_fontsize = 14;
    marker_size = 8;
    axis_extra = 0.5;
    
    
    for ivar = 1:length(resp_vars)
        nfig = nfig+1;
        figure(nfig), clf
        switch ivar
            case 1
                xvar = 'Mean Groove Rating';
                xdata = meanGrooveVect;
                xticklabel = {'least groove','','','','','','most groove'};
            case 2
                xvar = 'Mean Enjoyment Rating';
                xdata = meanEnjoyVect;
                xticklabel = {'not at all','','','','','','very much'};
            case 3
                xvar = 'Mean Want-to-continue Rating';
                xdata = meanContinueVect;
                xticklabel = {'not at all','','','','','','very much'};
        end
        set(gcf,'defaultaxesfontsize',axes_fontsize)
        ydata = meanTapVect;
        maxTaps = max(ydata);
        p = plot(xdata,ydata,'o','markersize', marker_size);
        lh = lsline;
        set(lh,'color',[1 0 0]);
        [rho,pval] = corr(xdata,ydata);
        text(1,maxTaps+10,sprintf('r=%1.2f, p<%1.4f', rho, max(0.0001,pval)), ...
            'horizontalalign','left', ...
            'fontsize', info_text_fontsize, ...
            'fontweight', info_text_fontweight);
        text(7,-10,sprintf('N = %d subjects', length(xdata)), ...
            'horizontalalign','right', ...
            'fontsize', info_text_fontsize, ...
            'fontweight', info_text_fontweight);
    
        xlabel(xvar,'fontsize',label_fontsize)
        ylabel(sprintf('Mean Number of Taps Per Trial'),'fontsize',label_fontsize)
        x_enum_values = 1:7;
        
        set(gca,'xtick', 1:length(x_enum_values),'xlim',[axis_extra length(x_enum_values)+axis_extra])
        set(gca,'xticklabel',xticklabel);
        set(gca,'ylim',[-20 maxTaps+20+axis_extra]);
        axis square

        if nfig == 1
            print(params.figfname_poststimTapCorr,'-dpsc','-adobecset');
        else
            print(params.figfname_poststimTapCorr,'-dpsc','-adobecset','-append');
        end
        
    end % for ivar
end % if params.plot_corr

%% output results to table

resp_table.vars = {'data','column_labels','column_formats'};
resp_table.data{1} = resp_mtx;
resp_table.data{2} = [{'subject_id'}, resp_vars];
resp_table.data{3} = {'%s','%1.2f','%1.2f','%1.2f'};

% print the table to file
if isfield(params,'report')
    ensemble_display_table(resp_table,params.report);
end

% output struct

result_st.vars = resp_table.data{2};
result_st.data = resp_table.data{1};

return

