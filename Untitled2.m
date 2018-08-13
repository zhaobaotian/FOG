clc
clear

fname = spm_select();
D = spm_eeg_load(fname);
%% plot traces
plotECG(D.time,squeeze(D(:,:,1)),'AutoStackSignals',D.chanlabels)
figure
plot(squeeze(D(1,15,1:500,1)))

figure
hold on
pwelch(D(2,:,1))
legend

y = wgn(114850,1,1);
y = y';
figure
pwelch(y)
%% use spm gui to filter
% NOTE: set the eeg type before you filter
%% visualization for inspection
fname_f = spm_select()
D_f = spm_eeg_load(fname_f)
FilteredSignal = squeeze(D_f(1,:,1))';
plotECG(D_f.time,[FilteredSignal])
RawSignal = squeeze(D(1,:,1))';
% plot(RawSignal)
% hold on
% plot(squeeze(D(1,:,1)),'--');
plotECG(D_f.time,[FilteredSignal RawSignal])

%% Time¨CFrequency analysis and denoising %%
% short time fourier transform of the filtered events for each event

% % from article
% The STFT was computed in a 256 ms segment starting 128 ms
% before the center of HFO and extending 128 ms after it. The
% Fourier transform was computed in a 64 ms Hanning window
% which was shifted sample by sample to create a time¨C
% frequency map.
% for i = 1%:length(length(WinCentersFilteredCut)
    
    Sig = D(4,[501:1500],1);
    HanningWin = hann(64*D.fsample/1000);
    Noverlap = length(HanningWin) - 1;
    frequencies = 1:200;
    fsample = D.fsample;
   figure
   spectrogram(Sig,HanningWin,Noverlap,frequencies,fsample,'yaxis')
   colormap('jet')
    
    
   [s,w,t,ps] = spectrogram(Sig,HanningWin,Noverlap,frequencies,fsample,'yaxis');
   figure
   spectrogram(Sig,HanningWin,Noverlap,frequencies,fsample,'yaxis')
   colormap('jet')

   % use threshold to keep 95% energy
   % Calculate the 95% percentile energy threshold from ps
   NormPs = 10*log10(ps + eps);
   NormPsVector = reshape(NormPs,[],1);
   PercentThreshold = prctile(NormPsVector,5);
   [sThresh,w,t,psThresh] = spectrogram(Sig,HanningWin,Noverlap,frequencies,fsample,'MinThreshold',PercentThreshold,'yaxis');
   
%    % for visualization
%    figure
%    imagesc(10*log10(psThresh + eps),[-40,25])
%    colormap('jet')
   
   % do the threshold manually
%    for i = 1:length(NormPsVector)
%        if NormPs(i) <  PercentThreshold
%           psThresh(i) = 0;
%        else 
%           psThresh(i) = ps(i);
%        end
%    end
%    
%    NormPsThreshold = 10*log10(psThresh + eps);
%    figure
%    imagesc(NormPsThreshold,[-40,25]) % plus eps to avoid log10(0)
%    set(gca,'Ydir','normal')
% %    imagesc(psThresh,[-40,25])
%    colormap jet
% end