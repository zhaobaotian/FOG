channel_pair = [1 5];

cohtimerange = [80 200];

coh_data = chan_data(channel_pair,cohtimerange(1)*fs:cohtimerange(2)*fs);

epochs = zeros([1 size(coh_data)]);

epochs(1,:,:)=coh_data;

data4freq     = transform2ft_freq(epochs,fs);


freqvec       = [2:0.5:35];
min_freq      = min(freqvec);
extndtime     = 2.5/min_freq;
max_time      = length(data4freq.trial{1,1})/fs;
time_interval = 0.005;

cfg              = [];
cfg.method       = 'wavelet';
cfg.width        = 5;
cfg.output       = 'powandcsd';
cfg.pad          = 'nextpow2';
cfg.keeptrials   = 'yes';

cfg.feedback     = 'no';
cfg.polyremoval  = 1;

cfg.foi          = freqvec; %[min_freq:0.3:6 7:1:20];

cfg.toi          = [extndtime:time_interval:max_time-extndtime];

tfr              = ft_freqanalysis(cfg, data4freq);

powspctrm = tfr.powspctrm;

crsspctrm_tmp = squeeze(tfr.crsspctrm);

crsspctrm = zeros([1 size(crsspctrm_tmp)]);

crsspctrm(1,:,:) = crsspctrm_tmp;



freqvector = tfr.freq;

timevector = tfr.time;

cyclenum = 3;

[ave_coh_otime,~] = time_coherence(powspctrm,crsspctrm,cyclenum,freqvector,timevector,time_interval);

coh_mean = mean(ave_coh_otime,2);

coh_ste = std(ave_coh_otime,0,2)./sqrt(size(ave_coh_otime,2));


figure;

boundedline(freqvector,coh_mean,coh_ste);