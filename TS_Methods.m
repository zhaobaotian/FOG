%% convert data to SPM format
D = DBS_Signal_SPM_Convert();
D= spm_eeg_load();
% figure;
plotECG(D.time,D(1,:))
% data  = D.time,D(1,:);
% figure
% pwelch(data,500,250,[],500)
%% filter the signal 4 times
% 1st: highpass 1Hz
% 2nd, 3rd, 4th, 5th are bandstop around 50Hz, 100Hz, 150Hz, 200Hz
% respectively
fname = fullfile(D.path,D.fname);
D = HighPass_50Hz_Filter(fname);

fname = spm_select()
D_f = spm_eeg_load(fname)

fname = spm_select()
D_TF = spm_eeg_load(fname)

figure;
imagesc((squeeze(D_TF(2,:,1:1500,1))))

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
SmoothedPowerSPD = imgaussfilt(powerspd,[2 125]);

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
input(:,:,1) = D(2,:,1)';
input(:,:,2) = D(4,:,1)';
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


fname = spm_select()
D_5_40 = spm_eeg_load(fname)
plotECG(D_5_40.time,D_5_40(4,:,1))
[envelope_lineup,~] = envelope(D_5_40(4,:,1),5,'peak');
plotECG(D_5_40.time,[D_5_40(4,:,1)' envelope_line'])
threshold = prctile(envelope_lineup,70);
thresholdline = repmat(threshold,[length(envelope_lineup),1]);
plotECG(D_5_40.time,[D_5_40(4,:,1)' envelope_lineup' thresholdline])







