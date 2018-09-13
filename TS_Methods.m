% Scripts to replicate the signal analysis methods from paper:
% Neumann, Wolf-Julian, et al. "Pallidal and thalamic neural oscillatory 
% patterns in Tourette syndrome." Annals of neurology (2018).
clear;
SubjectFolder = 'D:\FOG\Data\Dystonia001';
cd(SubjectFolder)
% Functional Neurosurgery Department of Tian Tan Hospital
%% -----------------DEFINE PARAMETERS HERE!!!----------------- %%
%--------------------------------------------------------------------%
% Define your interested timewindow here and the unseleted time series data
% will be discarded during format conversion
ROITimeWindow           = [0 60000]; % in ms
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
ChannelLabels           = {'LD_A_1';'LD_A_2';'LD_A_3';'LD_B_1';'LD_B_2';'LD_B_3'};
%--------------------------------------------------------------------%
% Frequency band of interest for power spectra plot
ThetaBand               = [3 12];  % theta band
BetaBand                = [13 35]; % beta band
%--------------------------------------------------------------------%
% Theta bursts detection threshold: 75th percentile of averged power of 
% theta band (3-12Hz)
AvgBandPower            = [3 12];  % theta band
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
NormMethod              = 1;       % 1: Normalize to the std of the Band Power spectrum
                                   % 2: Normalize to the std of the Band
                                   %    Power maxtrix
%--------------------------------------------------------------------%
% Visually check indicator
VisualCheck            = 0;        % 1: yes; 0: no
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
    switch NormMethod
        case 1
            STD = std(PowerSpectrum([5:45 55:95]));
        case 2
            Power  = squeeze(TFMatrix(i,[5:45 55:95],:));
            STD    = std(reshape(Power,[numel(Power) 1]));
    end
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
for i = 1:D_tf.nchannels
    figure
    findpeaks(PowerSpectrumNormalized(i,:))
end
% Pending, To be contunued...
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
    ThetaPowerEnvelope_temp = (ThetaPowerEnvelope(i,:));
    BurstInd                = find(ThetaPowerEnvelope_temp >= Threshold75th(i));
    ThetaBurstGap           = find(diff(BurstInd)>1);
    ThetaBurstMatrix        = [];
    ThetaBurst_counter      = 0;
    for j = 1:length(ThetaBurstGap)
        if j == 1
            onset            = BurstInd(1);
            offset           = BurstInd(ThetaBurstGap(j));
            %             PeakValue        = max(ThetaPowerEnvelope_temp(onset:offset));
            ThetaBurstLength = ((offset - onset)/D_tf.fsample)*1000; % in ms
        elseif j == length(ThetaBurstGap)
            onset            = BurstInd(ThetaBurstGap(j)+1);
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
            TimeInterval(int16((timetamps_temp(j,1)/1000)*500):int16((timetamps_temp(j,3)/1000)*500)) = max(Signal_temp);
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
%% Alligned theta burst for time frequency presentation
% Time frequency representations were smoothed with a full width half 
% maximum gaussian smoothing kernel of 2 Hz and 250 ms length for burst 
% analysis.
TFMatrixSmoothed = zeros(size(TFMatrix));
for i = 1:size(TFMatrix,1)
    TF_temp = squeeze(TFMatrix(i,:,:));
    TFMatrixSmoothed(i,:,:) = imgaussfilt(TF_temp,...
                              [FrequencySmooth TimeAxisSmooth/1000*D_Theta.fsample]);
end

% Extract TF data of theta burst from -500ms to 1500ms
ThetaBurstTF = cell(length(D_tf.nchannels),1);
for i = 1:D_tf.nchannels
    % temporal variables for the current channel
    ThetaPowerEnvelope_temp = (ThetaPowerEnvelope(i,:));
    BurstInd                = find(ThetaPowerEnvelope_temp >= Threshold75th(i));
    ThetaBurstGap           = find(diff(BurstInd)>1);
    ThetaBurstMatrix        = [];
    ThetaBurst_counter      = 0;
    for j = 1:length(ThetaBurstGap)
        if j == 1
            onset            = BurstInd(1);
            offset           = BurstInd(ThetaBurstGap(j));
            %             PeakValue        = max(ThetaPowerEnvelope_temp(onset:offset));
            ThetaBurstLength = ((offset - onset)/D_tf.fsample)*1000; % in ms
        elseif j == length(ThetaBurstGap)
            onset            = BurstInd(ThetaBurstGap(j)+1);
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


% ROI = [1000:3500];
% figure
% imagesc(log10(squeeze(D_tf(2,:,ROI,1))))
% colormap jet
% axis xy
% 
% figure
% imagesc(log10(squeeze(TFMatrixSmoothed(2,:,ROI))))
% colormap jet
% axis xy

%% 


figure;
imagesc((squeeze(D_TF(2,:,1:1500,1))));

ROI = [1000:3500];
figure
imagesc(log10(squeeze(D_TF(2,:,ROI,1))))
colormap jet
axis xy
figure
imagesc(log10(SmoothedPowerSPD(:,ROI)))
colormap jet
axis xy

% plotECG(D_f(4,:,1))

% pwelch(D_f(4,:,1),50,25,[],500)
powerspd = squeeze(D_TF(2,:,:,1));
SmoothedPowerSPD = imgaussfilt(powerspd,[FrequencySmooth TimeAxisSmooth/1000*500]);

figure
plot(mean(SmoothedPowerSPD,2))


%% Imaginary part coherence
% The input data input should be an array organized as:
%   Repetitions x Channel x Channel (x Frequency) (x Time)
% or
%   Repetitions x Channelcombination (x Frequency) (x Time)
%

[c, v, outcnt] = ft_connectivity_corr(input, varargin);


%% theta (3-12Hz) burst detection
theta_trace = spm_eeg_load();
thetaBand = LowerBandNormalized(8,:);
ThetaPower = thetaBand;
thresholdline = prctile(ThetaPower,75);
plotECG(theta_trace.time,[theta_trace(2,:,1)' (ThetaPower.*200)' repmat(thresholdline*200,[length(ThetaPower),1])])

% figure
% plot(ThetaPower)

%%
plotECG(theta_trace.time,theta_trace(2,:,1))
[envelope_lineup,~] = envelope(theta_trace(2,:,1),300,'peak');
plotECG(theta_trace.time,[theta_trace(2,:,1)' envelope_lineup'])
threshold = prctile(envelope_lineup,75);
thresholdline = repmat(threshold,[length(envelope_lineup),1]);
plotECG(theta_trace.time,[theta_trace(2,:,1)' envelope_lineup' thresholdline])

% A= [1:6]
% std(A)
% A_nornalized = A/ans

D_TF.chanlabels


fname = spm_select();
D_5_40 = spm_eeg_load(fname);
plotECG(D_5_40.time,D_5_40(4,:,1));
[envelope_lineup,~] = envelope(D_5_40(4,:,1),5,'peak');
plotECG(D_5_40.time,[D_5_40(4,:,1)' envelope_line'])
threshold = prctile(envelope_lineup,70);
thresholdline = repmat(threshold,[length(envelope_lineup),1]);
plotECG(D_5_40.time,[D_5_40(4,:,1)' envelope_lineup' thresholdline])







