clear all;

%% define some parameters
% task = 'Before_Drug';
% % define parameters
% fieldepoch = 'start';
% twepoch = [-300 1100];
% bc = 1;
% bcfield = 'start';
% twbc = [-200 0];
% smoothwin = 50;
% twResc = [-200 0];
% twsmooth = [-200 1000];
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

figure
pwelch(D_f(1,:,1))