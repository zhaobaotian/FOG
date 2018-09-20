
for i = 1:D_HighPass_Notch.nchannels
    for j = 1:D_HighPass_Notch.nchannels
        % Calculate the imaginary part coherence
        if isequal(i,j)
            continue
        end
            cfg            = [];
            cfg.method     = 'coh';
            cfg.complex    = 'imag';
            cfg.channelcmb = {iCOHData.label{i} iCOHData.label{j}};
            fdfourier_iCOH = ft_connectivityanalysis(cfg, freqfourier_iCOH);
            
            cfg                  = [];
            cfg.parameter        = 'cohspctrm';
            cfg.xlim             = [1 37];
            cfg.refchannel       = iCOHData.label{i};
            cfg.showlabels       = 'yes';
            cfg.channel = iCOHData.label{j};
            ft_singleplotER(cfg, fdfourier_iCOH);
            
            x = abs(fdfourier_iCOH.cohspctrm);
            windowSize = 150;
            b = (1/windowSize)*ones(1,windowSize);
            a = 1;
            y = filtfilt(b,a,x);
                
        figure
        set(gcf,'Color',[1 1 1])
        BarHeight    = 0.8;
        BarWith      = [10 23];
        BarPosition  = [3 13];
        BarColor     = [0.85 0.85 0.85;0.94 0.94 0.94];
        bar(4:11,repmat(BarHeight,[8 1]),2,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none')
        hold on
        bar(14:34,repmat(BarHeight,[21 1]),2,'FaceColor',[0.94 0.94 0.94],'EdgeColor','none')
        plot(fdfourier_iCOH.freq,y,'Color',[0.8 0.2 0.247],'LineWidth',4)
        set(gca,'FontSize',14);
        xlabel('Frequency [Hz]', 'FontSize', 18)
        ylabel('Imaginary part of coherence', 'FontSize', 18)
        ylim([0 0.8])
        xlim([0 37])
        xticks(0:10:30)
        title([D_HighPass_Notch.chanlabels{i},'_',D_HighPass_Notch.chanlabels{j}],...
            'Interpreter', 'none', 'FontSize', 20)
        print([D_HighPass_Notch.chanlabels{i},'_',D_HighPass_Notch.chanlabels{j},...
            '_' 'iCOH'],'-dpng','-r300')
        close
    end
end