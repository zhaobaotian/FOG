%% convert data to SPM format
D = DBS_Signal_SPM_Convert();

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
imagesc((squeeze(D_TF(6,:,1:1500,1))))
imagesc((squeeze(D_TF(6,:,1000:3500,1))))
colormap jet
axis xy
plotECG(D_f(4,:,1))

figure
pwelch(D_f(4,:,1),50,25,[],500)
powerspd = squeeze(D_TF(4,:,:,1));
plot(mean(powerspd,2))

D_TF.chanlabels



fname = spm_select()
D_5_40 = spm_eeg_load(fname)
plotECG(D_5_40.time,D_5_40(4,:,1))
[envelope_lineup,~] = envelope(D_5_40(4,:,1),5,'peak');
plotECG(D_5_40.time,[D_5_40(4,:,1)' envelope_line'])
threshold = prctile(envelope_lineup,70);
thresholdline = repmat(threshold,[length(envelope_lineup),1]);
plotECG(D_5_40.time,[D_5_40(4,:,1)' envelope_lineup' thresholdline])







