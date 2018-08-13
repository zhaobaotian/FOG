%-----------------------------------------------------------------------
% Job saved on 17-Jun-2018 16:03:41 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.preproc.filter.D = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.preproc.filter.type = 'butterworth';
matlabbatch{1}.spm.meeg.preproc.filter.band = 'high';
matlabbatch{1}.spm.meeg.preproc.filter.freq = 1;
matlabbatch{1}.spm.meeg.preproc.filter.dir = 'twopass';
matlabbatch{1}.spm.meeg.preproc.filter.order = 5;
matlabbatch{1}.spm.meeg.preproc.filter.prefix = 'f';
matlabbatch{2}.spm.meeg.preproc.filter.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{2}.spm.meeg.preproc.filter.type = 'fir';
matlabbatch{2}.spm.meeg.preproc.filter.band = 'stop';
matlabbatch{2}.spm.meeg.preproc.filter.freq = [47 53];
matlabbatch{2}.spm.meeg.preproc.filter.dir = 'twopass';
matlabbatch{2}.spm.meeg.preproc.filter.order = 64;
matlabbatch{2}.spm.meeg.preproc.filter.prefix = 'f';
matlabbatch{3}.spm.meeg.preproc.filter.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{3}.spm.meeg.preproc.filter.type = 'fir';
matlabbatch{3}.spm.meeg.preproc.filter.band = 'stop';
matlabbatch{3}.spm.meeg.preproc.filter.freq = [97 103];
matlabbatch{3}.spm.meeg.preproc.filter.dir = 'twopass';
matlabbatch{3}.spm.meeg.preproc.filter.order = 64;
matlabbatch{3}.spm.meeg.preproc.filter.prefix = 'f';
matlabbatch{4}.spm.meeg.preproc.filter.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{4}.spm.meeg.preproc.filter.type = 'fir';
matlabbatch{4}.spm.meeg.preproc.filter.band = 'stop';
matlabbatch{4}.spm.meeg.preproc.filter.freq = [147 153];
matlabbatch{4}.spm.meeg.preproc.filter.dir = 'twopass';
matlabbatch{4}.spm.meeg.preproc.filter.order = 64;
matlabbatch{4}.spm.meeg.preproc.filter.prefix = 'f';
matlabbatch{5}.spm.meeg.preproc.filter.D(1) = cfg_dep('Filter: Filtered Datafile', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{5}.spm.meeg.preproc.filter.type = 'fir';
matlabbatch{5}.spm.meeg.preproc.filter.band = 'stop';
matlabbatch{5}.spm.meeg.preproc.filter.freq = [197 203];
matlabbatch{5}.spm.meeg.preproc.filter.dir = 'twopass';
matlabbatch{5}.spm.meeg.preproc.filter.order = 64;
matlabbatch{5}.spm.meeg.preproc.filter.prefix = 'f';
