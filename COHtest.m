%% Imaginary part coherence
% Add fieldtrip path
[SPMpath,~,~] = fileparts(which('spm'));
addpath([SPMpath filesep 'external' filesep 'fieldtrip'])

% For simulation
rng default

Fs = 1000;
t = 0:1/Fs:1-1/Fs;

x = cos(2*pi*100*t) + sin(2*pi*200*t) + 0.5*randn(size(t));
y = 0.5*cos(2*pi*100*t-pi/4) + 0.35*sin(2*pi*200*t-pi/2) + 0.5*randn(size(t));
figure
subplot(2,1,1); plot(x)
subplot(2,1,2); plot(y)

xtf = fft(x);
xtf = xtf(1:500);

ytf = fft(y);
ytf = ytf(1:500);
xy=xcorr(xtf,ytf);

xy=xcorr(fft(x),fft(y));
xx=abs(fft(x));
xx = xx(1:500)*2/1000;
yy=abs(fft(y));
yy = yy(1:500)*2/1000;
figure
subplot(2,1,1); plot(xx)
subplot(2,1,2); plot(yy)
imag_coherence=xy(1:500)./sqrt(xx.*yy);

figure
plot(abs(imag_coherence))
 

% divide data in N equal length blocks for averaging later on
N = 1;
L  = floor(length(x)/N);
xt = reshape(x(1:L*N), L, N);
yt = reshape(y(1:L*N), L, N);

% transform to frequency domain
Xf = fft(xt,L,1);
Xf = Xf(1:500)*2/1000;
Yf = fft(yt,L,1);
Yf = Yf(1:500)*2/1000;

% estimate expectations by taking the average over N blocks
xy = Xf .* conj(Yf);
xx = Xf .* conj(Xf);
yy = Yf .* conj(Yf);

% combine terms to get final result
result=xy./(xx.*yy);
figure
plot(abs(result))



% simulation data preparation
DataSim.label   = {'A';'B'};
DataSim.fsample = 1000;
DataSim.trial   = {[x;y]};
DataSim.time    = {t};

% tf transform
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [1 500];
cfg.tapsmofrq  = 20;
cfg.keeptrials = 'yes';
cfg.channel    = {'all'};
Simfreqfourier    = ft_freqanalysis(cfg, DataSim);
% figure;
% imagesc(abs(squeeze(Simfreqfourier.fourierspctrm(:,2,:)))')
% axis xy

cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = {'A' 'B'};
fdSim          = ft_connectivityanalysis(cfg, Simfreqfourier);


cfg                  = [];
cfg.parameter        = 'cohspctrm';
cfg.xlim             = [1 500];
cfg.refchannel       = 'A';
% cfg.layout           = 'CTF151_helmet.mat';
cfg.showlabels       = 'yes';
% figure; ft_multiplotER(cfg, fd)
cfg.channel = 'B';
figure; ft_singleplotER(cfg, fdSim);




% example
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel    = {'MEG' 'EMGlft' 'EMGrgt'};
freqfourier    = ft_freqanalysis(cfg, data);

figure;
imagesc(abs(squeeze(freqfourier.fourierspctrm(:,1,:)))')
axis xy


figure;
imagesc(abs(squeeze(freqfourier.fourierspctrm(:,1,:)))')
axis xy


cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = {'MEG' 'EMG'};
fd             = ft_connectivityanalysis(cfg, freq);


cfg                  = [];
cfg.parameter        = 'cohspctrm';
cfg.xlim             = [5 80];
cfg.refchannel       = 'EMGlft';
% cfg.layout           = 'CTF151_helmet.mat';
cfg.showlabels       = 'yes';
% figure; ft_multiplotER(cfg, fd)
cfg.channel = 'MRC21';
figure; ft_singleplotER(cfg, fd);

% built-in function for coherence
x = D_HighPass_Notch(1,:,1);
y = D_HighPass_Notch(4,:,1);
[cxy,f] = mscohere(x,y,hann(500),250,[],500);
figure
plot(f(1:39),sqrt(cxy(1:39)))