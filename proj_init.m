function proj_init
 
[proj_path,~,~] = fileparts(mfilename('fullpath'));

addpath(fullfile(proj_path,'misc'));

% initialise SPM if its not already
if isempty(which('spm'))
    init_spm;
end

% initialise torso_tools
if isempty(which('tt_path'))
    addpath(fullfile(proj_path,'torso_tools'));
    tt_add_bem;
end

% initalise FEM tools and iso2mesh
if isempty(which('iso2meshver'))
    addpath(genpath(fullfile(proj_path,'fem_tools')));
end

