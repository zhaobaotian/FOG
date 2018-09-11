% Scripts to replicate the signal analysis methods from paper:
% Neumann, Wolf-Julian, et al. "Pallidal and thalamic neural oscillatory 
% patterns in Tourette syndrome." Annals of neurology (2018).
clear
% Functional Neurosurgery Department of Tian Tan Hospital
%% -----------------DEFINE PARAMETERS HERE!!!----------------- %%
%--------------------------------------------------------------------%
% Define your interested timewindow here and the unseleted time series data
% will be discarded during format conversion
ROITimeWindow           = [0 150000]; % in ms
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
ROIFreqBand{1}          = [3 12];  % theta band
ROIFreqBand{2}          = [13 35]; % beta band
%--------------------------------------------------------------------%
% Theta bursts detection threshold: 75th percentile of averged power of 
% theta band (3-12Hz)
AvgBandPower            = [3 12];  % theta band
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
%----------------------PARAMETERS DONE-------------------------------%

%% convert data to SPM format
% Initialize SPM EEG module
spm('defaults', 'EEG');
D_Raw = DBS_Signal_SPM_Convert([],srate,ChannelLabels,ROITimeWindow);

% for visually check the raw data
for i = 1:size(D_Raw,1)
    plotECG(D_Raw.time,D_Raw(i,:,1))
end
plotECG(D_Raw.time,D_Raw(:,:,1)',...
    'AutoStackSignals',D_Raw.chanlabels)

%% Baseline correction and detrending
% baseline correction using spm
clear S
S.D  = D_Raw;
D_bc = spm_eeg_bc(S);

% for visually check baseline corrected data
for i = 1:size(D_bc,1)
    plotECG(D_bc.time,D_bc(i,:,1))
end
plotECG(D_bc.time,D_bc(:,:,1)',...
    'AutoStackSignals',D_bc.chanlabels)
% High pass filter as detrending
clear S
S.D      = D_bc;
S.band   = 'high';
S.freq   = HighPassFrequency;
S.order  = 3;
S.prefix = 'HP_';
D_HighPass = spm_eeg_filter(S);

% for visually check highpass data
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

%% Time frequency decomposition using SPM wavelet
clear S
S.D           = D_HighPass;
S.channels    = 'all';
S.frequencies = ROIFrequencies;
S.timewin     = [-Inf Inf];
S.method      = 'morlet';
S.settings.ncycles     = 10;
[D_tf,~] = spm_eeg_tf(S);

%% Manipulate the TF power matrix
% TFMatrix = channels * frequencies * samples
TFMatrix = D_tf(:,:,:,1);

% Smooth the TFMatrix channel by channel
TFMatrixSmoothed = zeros(size(TFMatrix));
for i = 1:size(TFMatrix,1)
    TF_temp = squeeze(TFMatrix(i,:,:));
    TFMatrixSmoothed(i,:,:) = imgaussfilt(TF_temp,...
                              [FrequencySmooth TimeAxisSmooth/1000*500]);
end

% Normalize the power spectra
Power  = squeeze(TFMatrixSmoothed(1,[5:45 55:95],:));
STD = std(reshape(Power,[numel(Power) 1]));

% 
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


% Averaging over time resulted in resting power spectra
PowerSpectrum = mean(squeeze(TFMatrixSmoothed(1,:,:)),2);
figure
plot(PowerSpectrum(1:40)/STD,'Color',[0.8 0.2 0.247],'LineWidth',3)
xlabel('Frequency [Hz]')
ylabel('Relative spectral power [a.u.]')
%%
D = HighPass_50Hz_Filter(fname);
data  = D_f(1,:,1);
plotECG(D_f.time,data')
figure
PSDWelch = pwelch(D(3,:,1),500,250,[],500);
figure
plot(log10(PSDWelch))
findpeaks(PSDWelch)
findpeaks(PSDWelch,'MinPeakProminence',0.00001,'Annotate','extents')

fname = spm_select();
D_f = spm_eeg_load(fname);

fname = spm_select();
D_TF = spm_eeg_load(fname);

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

%% normalize to std
LowerBandPower = SmoothedPowerSPD(1:45,:);
HigherBandPower = SmoothedPowerSPD(55:95,:);
LowerBandPowerSTD = std(LowerBandPower,0,2);
HigherBandPowerSTD = std(HigherBandPower,0,2);

LowerBandNormalized = LowerBandPower./LowerBandPowerSTD;
HigherBandNormalized = HigherBandPower./HigherBandPowerSTD;
figure
imagesc((LowerBandNormalized(:,5000:10000)))
colormap jet
axis xy

Spectra_data = mean(LowerBandNormalized(:,:),2);
figure
plot(Spectra_data)
xticks(1:45)
xticklabels(1:45)

%% COHERENCE
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







