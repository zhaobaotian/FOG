function D = DBS_Signal_SPM_Convert(fname,srate,ChannelLabels,ROITimeWindow)
% This function is used to convert the EEG LFP signal collected from PD DBS
% patient to SPM M/EEG format.
if nargin<1 || isempty(fname)
    fname = spm_select(inf,'any','Select file to import',{},pwd,'.txt');
end
%% import and organize the txt EEG file
EEG = DBS_LFP_Import(fname,srate,ChannelLabels,ROITimeWindow);

%% convert it to SPM MEEG object for process later
% Initialize SPM
%--------------------------------------------------------------------------
spm('defaults','EEG');

% Some details about the data
%--------------------------------------------------------------------------
Nchannels = length(EEG.chanlabels);
Nsamples  = length(EEG.data);
Ntrials   = 1;
TimeOnset = 0; % in sec
Fsample = EEG.srate;

chlabels = EEG.chanlabels;

% define the output file name
%--------------------------------------------------------------------------
if chlabels{1}(1) == 'R'
fname = strcat('Right_',EEG.filename);
elseif chlabels{1}(1) == 'L'
fname = strcat('Left_',EEG.filename);
end
% create data array 
%--------------------------------------------------------------------------
data = EEG.data';

% create the time axis (should be the same for all trials)
%--------------------------------------------------------------------------
timeaxis = [0:(Nsamples-1)]./Fsample + TimeOnset;

% Create the Fieldtrip raw struct
%--------------------------------------------------------------------------

ftdata = [];

for i = 1:Ntrials
   ftdata.trial{i} = squeeze(data(:, :, i));
   ftdata.time{i} = timeaxis;
end


ftdata.fsample = Fsample;
ftdata.label = chlabels;
ftdata.label = ftdata.label(:);

% Convert the ftdata struct to SPM M\EEG dataset
%--------------------------------------------------------------------------
D = spm_eeg_ft2spm(ftdata, fname);

% Examples of providing additional information in a script
% [] comes instead of an index vector and means that the command
% applies to all channels/all trials.
%--------------------------------------------------------------------------
D = type(D, 'continuous');                        % Sets the dataset type
D = chantype(D, ':', 'EEG');                   % Sets the channel type 
% D = conditions(D, 1:Ntrials, 'Condition 1');  % Sets the condition label

% save
%--------------------------------------------------------------------------
save(D);