function EEG = DBS_LFP_Import(filename,srate)
% Function to import data of intracranial LFP collected from DBS patients
% 
% Inputs:
% - filename: name of data file in edf format (optional)
% - sodata  : name of SODATA .mat (optional)
% - bch     : indexes or names (cell array) of pathological channels to
%             exclude.
% Outputs :
% - EEG     : Organized            
% -------------------------------------------------------------------------

%% Get inputs
if nargin<1 || isempty(filename)
    filename = spm_select(inf,'any','Select file to import',{},pwd,'.txt');
end

if nargin<2 || isempty(srate)
    prompt = 'Please type in the sampling rate [e.g.500]:';
    srate = input(prompt);
end

%% import data
% file start with 7 is on the right side
% file start with 8 is on the left side
% merge the 2 txt files into one and name the channels

% rename channels
% naming convention
% Start each electrode name with RD_, RS_, RG_, LD_, LS, or LG_ followed by the name of the electrode (e.g., RD_RHD). 
% The first letter indicates which hemisphere the electrode lies in and the second letter indicates if the electrode is a depth, 
% subdural strip, or subdural grid electrode

EEG.data = [];
EEG.srate = srate;
[~,FileName,~] = fileparts(filename(1,:));
EEG.filename = FileName;
if str2double(FileName(1)) == 7
    EEG.chanlabels = {'RS_1';'RS_2';'RS_3';'RD_1';'RD_2';'RD_3'};
    data_temp = importdata(filename);
    data_temp = data_temp(:,1:6);
    EEG.data = data_temp;
elseif str2double(FileName(1)) == 8
    EEG.chanlabels = {'LS_1';'LS_2';'LS_3';'LD_1';'LD_2';'LD_3'};
    data_temp = importdata(filename);
    data_temp = data_temp(:,1:6);
    EEG.data = data_temp;
end
















