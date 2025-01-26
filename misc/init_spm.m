function init_spm(path)

if nargin == 0
    path = fullfile(proj_path,'spm');
end

if isempty(which('spm'))
    disp('Initialising SPM...')
    restoredefaultpath;
    addpath(path);
    spm('defaults','eeg');
    spm_jobman('initcfg');

    fprintf('SPM environment loaded\n');

else
    disp('SPM already on path!')
end