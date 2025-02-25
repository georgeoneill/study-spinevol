clearvars
close all
clc

proj_init;

triaxial = 1; % Turn on or off triaxial sensors

files = [];
files.root = proj_path;

if triaxial
    files.results = fullfile(files.root,'results','3axis');
else
    files.results = fullfile(files.root,'results','1axis');
end

if ~exist(files.results,'dir')
    mkdir(files.results)
end


%% Load in geometries and convert to be FieldTrip ready

load(fullfile(files.root,'geometries','all_geometries_v2.mat'));

% boundaries

bnds(1).pos = mesh_wm.p;
bnds(1).tri = mesh_wm.e;
bnds(1).unit = 'm';

bnds(2).pos = mesh_bone.p;
bnds(2).tri = mesh_bone.e;
bnds(2).unit = 'm';

bnds(3).pos = mesh_blood.p;
bnds(3).tri = mesh_blood.e;
bnds(3).unit = 'm';

bnds(4).pos = mesh_lungs.p;
bnds(4).tri = mesh_lungs.e;
bnds(4).unit = 'm';

bnds(5).pos = mesh_torso.p;
bnds(5).tri = mesh_torso.e;
bnds(5).unit = 'm';

% sensors

if triaxial
    coils = coils_3axis;
else
    coils = coils_1axis;
end

ncoils = length(coils.p);
grad = [];
grad.coilpos = coils.p;
grad.coilori = coils.n;
grad.tra = eye(ncoils);
grad.chanunit = repmat({'T'}, ncoils, 1);
grad.unit = 'm';
for ii = 1:ncoils
    grad.label{ii} = sprintf('coil_%03d',ii);
end
grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'm');

% sources
src = [];
src.pos = sources_cent.p;
src.inside = ones(length(src.pos),1);
src.unit = 'm';



%% Forward modelling - Infinite Medium

models = [];
L = [];
exec_time = [];
ti = tic;

cfg                     = [];
cfg.method              = 'infinite';
vol                     = ft_prepare_headmodel(cfg);
vol.type                = 'infinite_currentdipole';
vol.unit                = 'm';


% calculate forwards
cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat(fwd_tmp.leadfield);
exec_time{end+1} = toc(ti);
models{end+1} = 'Inf';


%% Forward modelling - BIG SPHERE

ti = tic;

vol = [];
vol.r = 1.1012;
vol.o = [0.8104 -1.7138 0.1324];  
vol.cond = 1;
vol.type = 'singlesphere';
vol.unit = 'm';

% calculate forwards
cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat(fwd_tmp.leadfield);
exec_time{end+1} = toc(ti);
models{end+1} = 'LS';

%% Forward modelling - Nolte's Corrected Sphere

ti = tic;

id = 5; % Just need the torso mesh;

bnd = [];
bnd.tri         = bnds(id).tri;
bnd.pos         = bnds(id).pos;
bnd.unit        = bnds(id).unit;
bnd             = ft_convert_units(bnd,'m');

cfg                     = [];
cfg.method              = 'singleshell';
vol                     = ft_prepare_headmodel(cfg,bnd);
vol.skin_surface        = [];

% calculate forwards
cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat({fwd_tmp.leadfield{:}});
exec_time{end+1} = toc(ti);
models{end+1} = 'CS';

%% Forward modelling - 1 Tissue BEM

ti = tic;

id = 5; % Just need the torso mesh;

bnd = [];
bnd.tri         = bnds(id).tri;
bnd.pos         = bnds(id).pos;
bnd.unit        = bnds(id).unit;
bnd             = ft_convert_units(bnd,'m');

ci = 0.23;
co = 0;

cfg                     = [];
cfg.method              = 'bem_hbf';
cfg.conductivity        = [ci;co];
vol_bem1c               = ft_prepare_headmodel(cfg,bnd);

cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol_bem1c;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat({fwd_tmp.leadfield{:}});
exec_time{end+1} = toc(ti);
models{end+1} = '1c LC';

%% 3 Tissue BEM

ti = tic;

clear bnd

id = [3 4 5]; % Blood, Lungs, Torso

for ii = 1:numel(id)

    tmp = [];
    tmp.tri         = bnds(id(ii)).tri;
    tmp.pos         = bnds(id(ii)).pos;
    tmp.unit        = bnds(id(ii)).unit;
    tmp             = ft_convert_units(tmp,'m');

    bnd(ii) = tmp;

end

ci = [.62 .05 .23];
co = [.23 .23  0 ];

cfg                     = [];
cfg.method              = 'bem_hbf';
cfg.conductivity        = [ci;co];
vol_bem3c               = ft_prepare_headmodel(cfg,bnd);

cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol_bem3c;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat({fwd_tmp.leadfield{:}});
exec_time{end+1} = toc(ti);
models{end+1} = '3c LC';

%% 4 Tissue BEM

ti = tic;

cratio = 40;
clear bnd
id = [1 3 4 5]; % Spinal Cord, Blood, Lungs, Torso

for ii = 1:numel(id)

    tmp = [];
    tmp.tri         = bnds(id(ii)).tri;
    tmp.pos         = bnds(id(ii)).pos;
    tmp.unit        = bnds(id(ii)).unit;
    tmp             = ft_convert_units(tmp,'m');

    bnd(ii) = tmp;

end

ci = [.33 .62 .05 .23];
co = [.23 .23 .23 .0];

cfg                     = [];
cfg.method              = 'bem_hbf';
cfg.conductivity        = [ci;co];
vol_bem4c               = ft_prepare_headmodel(cfg,bnd);

cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol_bem4c;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat({fwd_tmp.leadfield{:}});
exec_time{end+1} = toc(ti);
models{end+1} = '4c LC';

%% 5 Tissue BEM

ti = tic;

clear bnd
id = [1 2 3 4 5];
for ii = 1:numel(id)

    tmp = [];
    tmp.tri         = bnds(id(ii)).tri;
    tmp.pos         = bnds(id(ii)).pos;
    tmp.unit        = bnds(id(ii)).unit;
    tmp             = ft_convert_units(tmp,'m');

    bnd(ii) = tmp;

end

ci = [0.33 0.33/cratio .62 .05 .23];
co = [0.33/cratio .23 .23 .23 0];

cfg                     = [];
cfg.method              = 'bem_hbf';
cfg.conductivity        = [ci;co];
vol_bem5c               = ft_prepare_headmodel(cfg,bnd);

cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol_bem5c;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat({fwd_tmp.leadfield{:}});
exec_time{end+1} = toc(ti);
models{end+1} = '5c LC';

%% 5 Tissue BEM with Isolated Source approach

ti = tic;

cfg                     = [];
cfg.method              = 'bem_hbf';
cfg.conductivity        = [ci;co];
cfg.isolatedsource       = 1;
vol_bem5ci               = ft_prepare_headmodel(cfg,bnd);

cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol_bem5ci;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat({fwd_tmp.leadfield{:}});
exec_time{end+1} = toc(ti);
models{end+1} = '5ci LC';

%% Load in Linear Galerkin BEM results

load(fullfile(files.root,'geometries','LFM_LG_upsampled.mat'));

if triaxial
    x = 2;
else
    x = 1;
end

fprintf('%s\n',labels_L_LG{x,1,1})
L{end+1} = L_LG{x,1,1};
models{end+1} = '5c LG';
exec_time{end+1} = NaN;

fprintf('%s\n',labels_L_LG{x,1,2})
L{end+1} = L_LG{x,1,2};
models{end+1} = '5ci LG';
exec_time{end+1} = NaN;


%% Generate FEM Results

ti = tic;

clear bnd
id = [1 2 3 4 5];
for ii = 1:numel(id)

    tmp = [];
    tmp.tri         = bnds(id(ii)).tri;
    tmp.pos         = bnds(id(ii)).pos;
    tmp.unit        = bnds(id(ii)).unit;
    tmp             = ft_convert_units(tmp,'m');

    bnd(ii) = tmp;

end

% to allow for more reliable meshing in the spinal cord, we generate a new
% mesh 4 mm inwards and add that to the meshing.
bnd(1) = go_add_concentric_mesh(bnd(1),0.004);
vol_fem = go_paperbnds2mesh_denserspine(bnd,2);

files.femdir = fullfile(files.root,'fem_tempdir');
% create fem proc director or purge contents (to stop it using the a
% pre-exisitng transfer matrix
if ~exist(files.femdir)
    mkdir(files.femdir);
else
    delete(fullfile(files.femdir,'*'));
end

cratio = 40;

S = [];
S.dir = files.femdir;
S.mesh = vol_fem;
S.grad = grad;
S.src = src;
S.cond = [0.33 0.33/cratio .62 .05 .23];
% For FIL-based workstations you need to put the executable in the wtcnapps
% folder, then point to the directory where it is stored. Uncomment the
% code below and set the the appropriate path!
%S.bindir = 'C:\wtcnapps\duneuro';

% If this works, get a cup of tea, call your mum, read a paper whilst this
% computes.
fwd_tmp = fem_calc_fwds(S);

L{end+1} = fwd_tmp;
exec_time{end+1} = toc(ti);
models{end+1} = 'FEM';


%% Save results

fname_results = fullfile(files.results,'LFM_all_central_sources.mat');
save(fname_results,'L','models','exec_time','src','grad');