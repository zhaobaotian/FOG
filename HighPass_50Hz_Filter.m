function fd = HighPass_50Hz_Filter(fname)
% List of open inputs
% Filter: File Name - cfg_files
nrun = 1; % enter the number of runs here
jobfile = {'D:\FOG\HighPass_50Hz_Filter_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = {fname}; % Filter: File Name - cfg_files
end
spm('defaults', 'EEG');
fd = spm_jobman('run', jobs, inputs{:});
