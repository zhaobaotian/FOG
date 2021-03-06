% Scripts to replicate the signal analysis methods from paper:
% Neumann, Wolf-Julian, et al. "Pallidal and thalamic neural oscillatory 
% patterns in Tourette syndrome." Annals of neurology (2018).

% Baotian Zhao @ Beijing 20180917
clear;
SubjectFolder = 'D:\FOG\Data\Dystonia002';
cd(SubjectFolder)
% Functional Neurosurgery Department of Tian Tan Hospital
%% -----------------DEFINE PARAMETERS HERE!!!----------------- %%
%--------------------------------------------------------------------%
% Define your interested timewindow here and the unseleted time series data
% will be discarded during format conversion
ROITimeWindow           = [0 120000]; % in ms
%--------------------------------------------------------------------%
% Define your interested frequencies for time frequency decomposition,
% please pay attention to the Nyquist Law.
ROIFrequencies          = [1:100];
%--------------------------------------------------------------------%
% Wavelet Cycles for TF transform through SPM wavelet transform
WaveletsCycles          = 10;      % Define wavelets cycles here 
                                   % (ballence between time and frequency resolution)
%--------------------------------------------------------------------%
% Define the sampling rate when you collecting the data.
srate                   = 500;
%--------------------------------------------------------------------%
% Define the channel names, COLUMN-WISE. There should be 6 channels
% representing 6 bipolar referenced traces from 2 depth DBS electrodes
% The first 3 channels are located in the STN and the following 3 channels
% are located in GPi
ChannelLabels           = {'LD_A_1';'LD_A_2';'LD_A_3';'LD_B_1';'LD_B_2';'LD_B_3'};
%--------------------------------------------------------------------%
% Frequency band of interest for power spectra plot
ThetaBand               = [3 12];  % theta band
BetaBand                = [13 35]; % beta band
%--------------------------------------------------------------------%
% Theta bursts detection threshold: 75th percentile of averged power of 
% theta band (3-12Hz)
AvgBandPower            = [3 12];  % theta band
ThetaTFFrequencies      = [1 30];  % frequencies for aligned theta burst TF map
ThresholdPercentage     = 75;      % 75%
BurstLengthThreshold    = 200;     % 200ms minimum
PreburstBaseline        = 500;     % 500ms baseline
BurstEpoch              = 1500;    % 1500ms Window for burst
%--------------------------------------------------------------------%
% Parameters for TF map smooth
FrequencySmooth         = 2;       % smooth along frequency 
TimeAxisSmooth          = 250;     % 250ms
%--------------------------------------------------------------------%
% Highpass filter detrending
HighPassFrequency       = 0.5;     % High pass filter cutoff frequency
%--------------------------------------------------------------------%
% 50Hz line noise filter
Notch50                 = 1;       % 0: NOT do the 50Hz line noise filter
                                   % 1: DO the 50Hz line noise filter
StopBand                = [47 53;97 103]; % Stop Band for linenoise notch filter 

%--------------------------------------------------------------------%
% Visually check indicator
VisualCheck            = 1;        % 1: yes; 0: no
%--------------------------PARAMETERS DONE---------------------------%

%% convert data to SPM format
% Initialize SPM EEG module
spm('defaults', 'EEG');
D_Raw = DBS_Signal_SPM_Convert([],srate,ChannelLabels,ROITimeWindow);

% for visually check the raw data
if VisualCheck
    for i = 1:size(D_Raw,1)
        plotECG(D_Raw.time,D_Raw(i,:,1))
    end
    plotECG(D_Raw.time,D_Raw(:,:,1)',...
        'AutoStackSignals',D_Raw.chanlabels)
end
%% Baseline correction and detrending
% baseline correction using spm
clear S
S.D  = D_Raw;
D_bc = spm_eeg_bc(S);

% for visually check baseline corrected data
if VisualCheck
    for i = 1:size(D_bc,1)
        plotECG(D_bc.time,D_bc(i,:,1))
    end
    plotECG(D_bc.time,D_bc(:,:,1)',...
        'AutoStackSignals',D_bc.chanlabels)
end
% High pass filter as detrending
clear S
S.D      = D_bc;
S.band   = 'high';
S.freq   = HighPassFrequency;
S.order  = 3;
S.prefix = 'HP_';
D_HighPass = spm_eeg_filter(S);

% for visually check highpass data
if VisualCheck
    for i = 1:size(D_HighPass,1)
        plotECG(D_HighPass.time,D_HighPass(i,:,1))
    end
    plotECG(D_HighPass.time,D_HighPass(:,:,1)',...
        'AutoStackSignals',D_HighPass.chanlabels)
    
    % for visually check PSD using pwelch method
    for i = 1:size(D_HighPass,1)
        figure
        pwelch(D_HighPass(i,:,1),500,250,[],srate);
    end
end
%% 50Hz line noise notch filter
if Notch50
    % 50 Hz notch filter
    clear S
    S.D              = D_HighPass;
    S.band           = 'stop';
    S.freq           = StopBand(1,:);
    S.order          = 3;
    S.prefix         = 'Notch50_';
    D_HighPass_Notch = spm_eeg_filter(S);
    % 100 Hz notch filter
    clear S
    S.D              = D_HighPass_Notch;
    S.band           = 'stop';
    S.freq           = StopBand(2,:);
    S.order          = 3;
    S.prefix         = 'Notch100_';
    D_HighPass_Notch = spm_eeg_filter(S);
end
% for visually check after notch filter
if VisualCheck
    for i = 1:size(D_HighPass_Notch,1)
        figure
        pwelch(D_HighPass_Notch(i,:,1),500,250,[],srate);
    end
end
%% Time frequency decomposition using SPM wavelet
clear S
if Notch50
    S.D       = D_HighPass_Notch;
else
    S.D       = D_HighPass;
end
S.channels    = 'all';
S.frequencies = ROIFrequencies;
S.timewin     = [-Inf Inf];
S.method      = 'morlet';
S.settings.ncycles     = 10;
[D_tf,~] = spm_eeg_tf(S);

%% Manipulate the TF power matrix
% TFMatrix = channels * frequencies * samples
TFMatrix = D_tf(:,:,:,1);

% Averaging over time resulted in resting power spectra
% that were normalized to the standard deviation of the 
% 5 - 45 Hz and 55 - 95 Hz band power
mkdir('PCS')
cd('PCS')
PowerSpectrumNormalized = zeros(D_tf.nchannels,40);
for i = 1:D_tf.nchannels
    PowerSpectrum = mean(squeeze(TFMatrix(i,:,:)),2);
    STD = std(PowerSpectrum([5:45 55:95]));

    figure
    set(gcf,'Color',[1 1 1])
    BarHeight    = 15;
    BarWith      = [10 23];
    BarPosition  = [3 13];
    BarColor     = [0.85 0.85 0.85;0.94 0.94 0.94];
    bar(4:11,repmat(BarHeight,[8 1]),2,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none')
    hold on
    bar(14:34,repmat(BarHeight,[21 1]),2,'FaceColor',[0.94 0.94 0.94],'EdgeColor','none')
    
    PowerSpectrumNormalized(i,:) = PowerSpectrum(1:40)/STD;
    plot(PowerSpectrumNormalized(i,:),'Color',[0.8 0.2 0.247],'LineWidth',4)
    set(gca,'FontSize',14);
    xlabel('Frequency [Hz]', 'FontSize', 18)
    ylabel('Relative spectral power [a.u.]', 'FontSize', 18)
    ylim([0 15])
    xlim([0 37])
    xticks(0:10:30)
    title(D_tf.chanlabels{i},'Interpreter', 'none', 'FontSize', 20,'Position',[5 15.1 0])
    print([D_tf.chanlabels{i} '_' 'Normalized_Power_Spectrum'],'-dpng','-r300')
    close
end
%% Find peaks on each channel power spectrum and aligned to the respective peaks
% print the peaks found for visually check
Spectrum_Peaks          = cell(D_tf.nchannels,1);
Spectrum_Peaks_Location = cell(D_tf.nchannels,1);
for i = 1:D_tf.nchannels
    figure
    set(gcf,'DefaultLineLineWidth',4)
    findpeaks(PowerSpectrumNormalized(i,:));
    [Spectrum_Peaks_temp, Spectrum_Peaks_Location_temp] = findpeaks(PowerSpectrumNormalized(i,:));
    set(gca,'FontSize',14);
    xlabel('Frequency [Hz]', 'FontSize', 18)
    ylabel('Relative spectral power [a.u.]', 'FontSize', 18)
    title(D_tf.chanlabels{i},'Interpreter', 'none', 'FontSize', 20)
    print([D_tf.chanlabels{i} '_' 'Peak_Detection'],'-dpng','-r300')
    close
    Spectrum_Peaks{i}          = Spectrum_Peaks_temp;
    Spectrum_Peaks_Location{i} = Spectrum_Peaks_Location_temp;
end

% extract the peak location and plot them as aligned power spectrum
for i = 1:D_tf.nchannels
    Spectrum_Peaks_Loc_temp   = Spectrum_Peaks_Location{i};
    Spectrum_Peaks_Value_temp = Spectrum_Peaks{i};
    % search for highest theta peak
    ThetaPeakLocIndex = find(Spectrum_Peaks_Loc_temp > 3 & Spectrum_Peaks_Loc_temp < 12);
    if length(ThetaPeakLocIndex) > 1
        [~, Idx_temp] = max(Spectrum_Peaks_Value_temp(ThetaPeakLocIndex));
        ThetaPeakLocIndex = ThetaPeakLocIndex(Idx_temp);
    end
    if isempty(ThetaPeakLocIndex)
        fprintf('%s\n',[D_tf.chanlabels{i} ': No Peaks were found at Theta Band on this channel'])
        ThetaPeakLocMatrix(i) = 0;
    elseif ~isempty(ThetaPeakLocIndex)
        fprintf('%s\n',[D_tf.chanlabels{i} ': Theta peak was found at ' num2str(Spectrum_Peaks_Loc_temp(ThetaPeakLocIndex)) 'Hz'])
        ThetaPeakLocMatrix(i) = Spectrum_Peaks_Loc_temp(ThetaPeakLocIndex);
    end
    % search for highest beta peak
    BetaPeakLocIndex = find(Spectrum_Peaks_Loc_temp > 13 & Spectrum_Peaks_Loc_temp < 35);
    if length(BetaPeakLocIndex) > 1
        [~, Idx_temp] = max(Spectrum_Peaks_Value_temp(BetaPeakLocIndex));
        BetaPeakLocIndex = BetaPeakLocIndex(Idx_temp);
    end
    if isempty(BetaPeakLocIndex)
        fprintf('%s\n',[D_tf.chanlabels{i} ': No Peaks were found at Beta Band on this channel'])
        BetaPeakLocMatrix(i) = 0;
    elseif ~isempty(BetaPeakLocIndex)
        fprintf('%s\n',[D_tf.chanlabels{i} ': Beta peak was found at ' num2str(Spectrum_Peaks_Loc_temp(BetaPeakLocIndex)) 'Hz'])
        BetaPeakLocMatrix(i) = Spectrum_Peaks_Loc_temp(BetaPeakLocIndex);
    end
end

% data preparation for writing in an excel file
% header preparation
ExcelFileName      = 'PeaksInfo.xlsx';
ExcelSheet         = 1;
TitlesOfEachColumn = {'Channel_Names', 'Theta_Peak(Hz)', 'Beta_Peak(Hz)'};
ChannelLabels      = D_tf.chanlabels';
% write the excel file
xlswrite(ExcelFileName,TitlesOfEachColumn,ExcelSheet,'A1')
xlswrite(ExcelFileName,ChannelLabels,ExcelSheet,'A2')
xlswrite(ExcelFileName,ThetaPeakLocMatrix',ExcelSheet,'B2')
xlswrite(ExcelFileName,BetaPeakLocMatrix',ExcelSheet,'C2')

% plot the aligned peaks
% for theta peaks
ThetaPeakCounter = 0;
for i = 1:D_tf.nchannels    
    if ThetaPeakLocMatrix(i) > 0
        ThetaPeakCounter = ThetaPeakCounter + 1;
        ThetaPeakMatrix(ThetaPeakCounter,:) = PowerSpectrumNormalized(i,ThetaPeakLocMatrix(i) - 3:ThetaPeakLocMatrix(i) + 8);
    end
end

% for beta peaks
BetaPeakCounter = 0;
for i = 1:D_tf.nchannels    
    if BetaPeakLocMatrix(i) > 0
        BetaPeakCounter = BetaPeakCounter + 1;
        BetaPeakMatrix(BetaPeakCounter,:) = PowerSpectrumNormalized(i,BetaPeakLocMatrix(i) - 5:BetaPeakLocMatrix(i) + 5);
    end
end

%% plot the aligned theta peaks
figure
set(gcf,'Color',[1 1 1]) % set background to white
% plot individual line first
for i = 1:size(ThetaPeakMatrix,1)
    plot(ThetaPeakMatrix(i,:),'Color',[0.8 0.8 0.8],'LineWidth',1)
    hold on
end
MeanThetaPeak = mean(ThetaPeakMatrix);
% STDThetaPeak  = std(ThetaPeakMatrix);
plot(MeanThetaPeak,'Color',[0.8 0.2 0.247],'LineWidth',4)
% plot([MeanThetaPeak + STDThetaPeak;MeanThetaPeak - STDThetaPeak]','Color',[0.8 0.2 0.247],'LineWidth',1)
axis tight
set(gca,'FontSize',12);
xlabel('Frequency [Hz]', 'FontSize', 14)
ylabel('Relative spectral power [a.u.]', 'FontSize', 14)
ylim([0 10])
xticklabels(-2:2:8)
% title(D_tf.chanlabels{i},'Interpreter', 'none', 'FontSize', 20)
text(0.3,0.95,'Theta','Units','normalized','FontSize',30)
set(gcf,'units','pixels','position',[0 0 250 600])
print(['Aligned_Theta_Peaks'],'-dpng','-r300')
close


%% plot the aligned beta peaks
figure
set(gcf,'Color',[1 1 1]) % set background to white
% plot individual line first
for i = 1:size(BetaPeakMatrix,1)
    plot(BetaPeakMatrix(i,:),'Color',[0.8 0.8 0.8],'LineWidth',1)
    hold on
end
MeanBetaPeak = mean(BetaPeakMatrix);
% STDThetaPeak  = std(BetaPeakMatrix);
plot(MeanBetaPeak,'Color',[0.8 0.2 0.247],'LineWidth',4)
% plot([MeanThetaPeak + STDThetaPeak;MeanThetaPeak - STDThetaPeak]','Color',[0.8 0.2 0.247],'LineWidth',1)
axis tight
set(gca,'FontSize',12);
xlabel('Frequency [Hz]', 'FontSize', 14)
ylabel('Relative spectral power [a.u.]', 'FontSize', 14)
ylim([0 5])
xticks(1:5:11)
xticklabels(-5:5:5)
% title(D_tf.chanlabels{i},'Interpreter', 'none', 'FontSize', 20)
text(0.3,0.95,'Beta','Units','normalized','FontSize',30)
set(gcf,'units','pixels','position',[0 0 250 600])
print(['Aligned_Beta_Peaks'],'-dpng','-r300')
close

%% theta bursts detection
% theta band filter for visualization later
clear S
if Notch50
    S.D          = D_HighPass_Notch;
else
    S.D          = D_HighPass;
end
S.band           = 'bandpass';
S.freq           = ThetaBand;
S.order          = 3;
S.prefix         = 'Theta_';
D_Theta = spm_eeg_filter(S);
% for visually check theta band traces
if VisualCheck
    for i = 1:size(D_Theta,1)
        plotECG(D_Theta.time,D_Theta(i,:,1))
    end
    plotECG(D_Theta.time,D_Theta(:,:,1)',...
        'AutoStackSignals',D_Theta.chanlabels)
end
% Calculate the power envelope and 75 percent threshold
ThetaPowerEnvelope = squeeze(mean(TFMatrix(:,ThetaBand(1):ThetaBand(2),:),2));
Threshold75th      = prctile(ThetaPowerEnvelope,ThresholdPercentage,2);

% for visually check the thresholds on each channel
if VisualCheck
    for i = 1:D_tf.nchannels
        plotECG(D_Theta.time,[D_Theta(i,:,1)' (ThetaPowerEnvelope(i,:))' ...
            repmat(Threshold75th(i),[D_Theta.nsamples,1])])
    end
end
% initialize the output
ThetaBurstTimestamps = cell(length(D_tf.nchannels),1); % one cell per channel

for i = 1:D_tf.nchannels
    % temporal variables for the current channel
    ThetaPowerEnvelope_temp = ThetaPowerEnvelope(i,:);
    BurstInd                = find(ThetaPowerEnvelope_temp >= Threshold75th(i));
    ThetaBurstGap           = find(diff(BurstInd)>1);
    ThetaBurstMatrix        = [];
    ThetaBurst_counter      = 0;
    for j = 1:length(ThetaBurstGap)+1
        if j == 1
            onset            = BurstInd(1);
            offset           = BurstInd(ThetaBurstGap(j));
            %             PeakValue        = max(ThetaPowerEnvelope_temp(onset:offset));
            ThetaBurstLength = ((offset - onset)/D_tf.fsample)*1000; % in ms
        elseif j == length(ThetaBurstGap)+1
            onset            = BurstInd(ThetaBurstGap(j-1)+1);
            offset           = BurstInd(end);
            %             PeakValue  = max(ThetaPowerEnvelope_temp(onset:offset));
            ThetaBurstLength = ((offset - onset)/D_tf.fsample)*1000; % in ms
        else
            onset = BurstInd(ThetaBurstGap(j-1)+1);
            offset = BurstInd(ThetaBurstGap(j));
            %             PeakValue = max(ThetaPowerEnvelope_temp(onset:offset));
            ThetaBurstLength = ((offset - onset)/D_tf.fsample)*1000; % in ms
        end
        if ThetaBurstLength > BurstLengthThreshold % exclude burst less than 200ms
            ThetaBurst_counter = ThetaBurst_counter + 1;
            ThetaBurstMatrix(ThetaBurst_counter,:) = ...
                [(onset/D_tf.fsample)*1000 ThetaBurstLength (offset/D_tf.fsample)*1000]; 
        end
    end
    ThetaBurstTimestamps{i} = ThetaBurstMatrix;
end

% for visulization and mannual check
if VisualCheck
    for i = 1:D_Theta.nchannels
        timetamps_temp = ThetaBurstTimestamps{i};
        Signal_temp = D_Theta(i,:,1);
        TimeInterval = zeros(D_Theta.nsamples,1);
        ThetaPowerEnvelope_temp = ThetaPowerEnvelope(i,:)';
        for j = 1:size(timetamps_temp,1)
            TimeInterval(int64((timetamps_temp(j,1)/1000)*500):int64((timetamps_temp(j,3)/1000)*500)) = max(Signal_temp);
        end
        plotECG(D_Theta.time,[Signal_temp' repmat(Threshold75th(i),[D_Theta.nsamples 1]) TimeInterval ThetaPowerEnvelope_temp])
    end
end

% data preparation for writing in an excel file
ExcelFileName      = 'ThetaBurstInfo.xlsx';
ExcelSheet         = 1;
ChannelLabels      = D_Theta.chanlabels;
ChannelLabelTitles = {ChannelLabels{1},[],[],ChannelLabels{2},[],[],...
                      ChannelLabels{3},[],[],ChannelLabels{4},[],[],...
                      ChannelLabels{5},[],[],ChannelLabels{6},[],[]};
ThetaBurstContents = repmat({'Onset(ms)','Duration(ms)','Offset(ms)'},[1 6]);

xlswrite(ExcelFileName,ChannelLabelTitles,ExcelSheet,'A1')
xlswrite(ExcelFileName,ThetaBurstContents,ExcelSheet,'A2')
xlswrite(ExcelFileName,ThetaBurstTimestamps{1},ExcelSheet,'A3')
xlswrite(ExcelFileName,ThetaBurstTimestamps{2},ExcelSheet,'D3')
xlswrite(ExcelFileName,ThetaBurstTimestamps{3},ExcelSheet,'G3')
xlswrite(ExcelFileName,ThetaBurstTimestamps{4},ExcelSheet,'J3')
xlswrite(ExcelFileName,ThetaBurstTimestamps{5},ExcelSheet,'M3')
xlswrite(ExcelFileName,ThetaBurstTimestamps{6},ExcelSheet,'P3')

% Bar plot of the Theta Burst duration distribution on each channel
for i = 1:D_tf.nchannels
    % the second row denotes the durations of theta bursts
    ThetaBurstTimestampsDuration_temp = ThetaBurstTimestamps{i}(:,2);
    nBins = 13;
    ThetaBarColors = repmat(linspace(180,76,13)/255,[3 1]);
    figure
    for j = 1:nBins
        if j < 13
            nThetaBurst = length(find(ThetaBurstTimestampsDuration_temp <= 300 + (j - 1) * 100 & ...
                ThetaBurstTimestampsDuration_temp >= 200 + (j - 1) * 100));
            bar(j,nThetaBurst,1,'FaceColor',ThetaBarColors(:,j))
            hold on
            %         errorbar(j,nThetaBurst,std(ThetaBurstTimestampsDuration_temp(find(ThetaBurstTimestampsDuration_temp < 300 + (j - 1) * 100 & ...
            %             ThetaBurstTimestampsDuration_temp >= 200 + (j - 1) * 100))))
        elseif j == 13
            bar(j,length(find(ThetaBurstTimestampsDuration_temp > 200 + (j - 1) * 100)),1,'FaceColor',ThetaBarColors(:,j))
        end
        ThetaBurstHist(i,j) = nThetaBurst;
    end
    xticks(1:13)
    xticklabels([(3:14)*100])
    ax = gca;
    ax.XTickLabel{13} = '>1400';
    ylim([0 max(ThetaBurstHist(i,:))+5]);
    set(gca,'FontSize',12);
    ylabel('N bursts','FontSize', 14)
    xlabel('Theta Burst length [ms]','FontSize', 14)
    xtickangle(30)
    title(D_tf.chanlabels{i},'Interpreter', 'none', 'FontSize', 16)
    set(gcf,'units','pixels','position',[0 0 700 250])
    print(['Theta_Bursts_Length_Histogram' '_' D_tf.chanlabels{i}],'-dpng','-r300')
    close
end

%% Alligned theta burst for time frequency presentation
% Time frequency representations were smoothed with a full width half 
% maximum gaussian smoothing kernel of 2 Hz and 250 ms length for burst 
% analysis.
% %% smooth before alignment seems incorrect
% TFMatrixSmoothed = zeros(size(TFMatrix));
% for i = 1:size(TFMatrix,1)
%     TF_temp = squeeze(TFMatrix(i,:,:));
%     TFMatrixSmoothed(i,:,:) = imgaussfilt(TF_temp,...
%                               [FrequencySmooth TimeAxisSmooth/1000*D_Theta.fsample]);
% end
% 
% Extract smoothed TF data of theta burst from -500ms to 1500ms
% ThetaBurstTF = cell(D_tf.nchannels,1);
% for i = 1:D_tf.nchannels
%     temporal variables for the current channel
%     ThetaBurstTF_temp         = squeeze(TFMatrixSmoothed(i,ThetaTFFrequencies(1):ThetaTFFrequencies(2),:));
%     ThetaBurstTFMatrix        = zeros(size(ThetaBurstTimestamps{i},1),30,1001);
%     for j = 1:size(ThetaBurstTimestamps{i},1)
%         if ThetaBurstTimestamps{i}(j,1) > PreburstBaseline &...
%                 ThetaBurstTimestamps{i}(j,1) + 1500 < ROITimeWindow(2) - ROITimeWindow(1)
%             ThetaBurstInterval        = int16((ThetaBurstTimestamps{i}(j,1) - PreburstBaseline)/1000*D_tf.fsample):...
%                 int16((ThetaBurstTimestamps{i}(j,1) + BurstEpoch)/1000*D_tf.fsample);
%             ThetaBurstTFMatrix(j,:,:) = ThetaBurstTF_temp(:,ThetaBurstInterval);
%         end
%     end
%     ThetaBurstTF{i} = ThetaBurstTFMatrix;
% end

% Extract TF data of theta burst from -500ms to 1500ms
ThetaBurstTF = cell(D_tf.nchannels,1);
for i = 1:D_tf.nchannels
    % temporal variables for the current channel
    ThetaBurstTF_temp          = squeeze(TFMatrix(i,ThetaTFFrequencies(1):ThetaTFFrequencies(2),:));
    ThetaBurstTFMatrix         = [];
    ThetaBurstTFMatrix_Counter = 0;
    for j = 1:size(ThetaBurstTimestamps{i},1)
        if ThetaBurstTimestamps{i}(j,1) > PreburstBaseline &&...
                ThetaBurstTimestamps{i}(j,1) + 1500 < ROITimeWindow(2) - ROITimeWindow(1)
            ThetaBurstTFMatrix_Counter = ThetaBurstTFMatrix_Counter + 1;
            ThetaBurstInterval        = int16((ThetaBurstTimestamps{i}(j,1) - PreburstBaseline)/1000*D_tf.fsample):...
                int16((ThetaBurstTimestamps{i}(j,1) + BurstEpoch)/1000*D_tf.fsample);
            ThetaBurstTFMatrix(ThetaBurstTFMatrix_Counter,:,:) = ThetaBurstTF_temp(:,ThetaBurstInterval);
        end
    end
    ThetaBurstTF{i} = ThetaBurstTFMatrix;
end

% % Power normalization to the baseline (percent change)
% % single epoch normalizaiton before the average
% ThetaBurstTFNormalization = cell(D_tf.nchannels,1);
% for i = 1:D_tf.nchannels
%     % temporal variables for the current channel
%     ThetaBurstTFNormalization_temp         = ThetaBurstTF{i};
%     ThetaBurstTFNormalizationMatrix        = zeros(size(ThetaBurstTFNormalization_temp));
%     for j = 1:size(ThetaBurstTFNormalization_temp,1)
%         BaselinePower = mean(squeeze(ThetaBurstTFNormalization_temp(j,:,1:int16(PreburstBaseline/1000*500))),2);
%         ThetaBurstTFNormalizationMatrix(j,:,:) = (squeeze(ThetaBurstTFNormalization_temp(j,:,:))-BaselinePower)./BaselinePower;
%         if VisualCheck
%             figure
%             contourf(squeeze(ThetaBurstTFNormalizationMatrix(j,:,:)),'LineStyle','none');
%         end
%     end
%     ThetaBurstTFNormalization{i} = ThetaBurstTFNormalizationMatrix;
% end
% % plot the averaged theta burst TF map
% for i = 1:D_tf.nchannels
%     figure
%     set(gcf,'Color',[1 1 1])
%     % Time frequency representations were smoothed with a full width half
%     % maximum gaussian smoothing kernel of 2 Hz and 250 ms length for burst
%     % analysis.
%     ThetaTFmap = imgaussfilt(squeeze(mean(ThetaBurstTFNormalization{i}))*100,[FrequencySmooth TimeAxisSmooth/1000*D_Theta.fsample]);
%     % No smooth
%     % ThetaTFmap = squeeze(mean(ThetaBurstTFNormalization{3}));
%     contourf(ThetaTFmap,40,'LineStyle','none');
%     set(gca,'FontSize',14);
%     xlabel('Time [s]', 'FontSize', 18)
%     ylabel('Frequency [Hz]', 'FontSize', 18)
%     xticks(1:250:1001)
%     xticklabels(-0.5:0.5:1.5)
%     ThetaColorBar = colorbar;
%     ThetaColorBar.Label.String = 'Relative Spectral Power [%]';
%     grid on
%     set(gca,'GridColor',[0.8 0.8 0.8])
%     title(D_tf.chanlabels{i},'Interpreter', 'none', 'FontSize', 20)
%     print([D_tf.chanlabels{i} '_' 'Theta_Burst_Averaged_TF_1'],'-dpng','-r300')
%     close
% end


% Power normalization to the baseline (percent change)
% Average the raw TF representations before the baseline normalization
ThetaBurstTFNormalization = cell(D_tf.nchannels,1);
for i = 1:D_tf.nchannels
    % temporal variables for the current channel
    ThetaBurstTFNormalization_temp         = ThetaBurstTF{i};
    ThetaBurstTFNormalizationMatrix        = zeros(size(ThetaBurstTFNormalization_temp));
    ThetaBurstTFNormalizationMean = squeeze(mean(ThetaBurstTFNormalization_temp));
    BaselinePower  = mean(ThetaBurstTFNormalizationMean(:,1:int16(PreburstBaseline/1000*500)),2);
    ThetaBurstTFNormalizationMatrix = squeeze((ThetaBurstTFNormalizationMean-BaselinePower)./BaselinePower);
    ThetaBurstTFNormalization{i} = ThetaBurstTFNormalizationMatrix;
end

% plot the averaged theta burst TF map
for i = 1:D_tf.nchannels
    figure
    set(gcf,'Color',[1 1 1])
    % Time frequency representations were smoothed with a full width half
    % maximum gaussian smoothing kernel of 2 Hz and 250 ms length for burst
    % analysis.
    ThetaTFmap = imgaussfilt(ThetaBurstTFNormalization{i}*100,[FrequencySmooth TimeAxisSmooth/1000*D_Theta.fsample]);
    contourf(ThetaTFmap,40,'LineStyle','none');
    set(gca,'FontSize',14);
    xlabel('Time [s]', 'FontSize', 18)
    ylabel('Frequency [Hz]', 'FontSize', 18)
    xticks(1:250:1001)
    xticklabels(-0.5:0.5:1.5)
    ThetaColorBar = colorbar;
    ThetaColorBar.Label.String = 'Relative Spectral Power [Percent Change (%)]';
    grid on
    set(gca,'GridColor',[0.8 0.8 0.8])
    title(D_tf.chanlabels{i},'Interpreter', 'none', 'FontSize', 20)
    print([D_tf.chanlabels{i} '_' 'Theta_Burst_Averaged_TF'],'-dpng','-r300')
    close
end


%% Imaginary part coherence
% % Add fieldtrip path
% [SPMpath,~,~] = fileparts(which('spm'));
% addpath([SPMpath filesep 'external' filesep 'fieldtrip'])


% The coherence values reflect the consistency of the phase difference between the two signals at a given frequency
% DATA preparation for the coherence calculation
% iCOHDataHeader = ft_read_header(spm_select(1,'mat')); % select the .mat file
iCOHData.label   = D_HighPass_Notch.chanlabels';
iCOHData.time    = {D_HighPass_Notch.time};
iCOHData.trial   = {D_HighPass_Notch(:,:,1)};
iCOHData.fsample = D_HighPass_Notch.fsample;

if VisualCheck
    figure
    for i = 1:D_HighPass_Notch.nchannels
        subplot(D_HighPass_Notch.nchannels,1,i);
        plot(iCOHData.time{1},iCOHData.trial{1}(i,:));
        axis tight;
        legend(iCOHData.label(i));
    end
end

% fourier tf transform
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [1 37];
cfg.tapsmofrq  = 1;
cfg.keeptrials = 'yes';
cfg.channel    = {'all'};
freqfourier_iCOH    = ft_freqanalysis(cfg, iCOHData);

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
     
        % smooth the output
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





