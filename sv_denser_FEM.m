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

%% Prepare bits to save

L = [];
models = [];

%% Now try and mesh this

vol_fem = go_paperbnds2mesh_denserspine(bnds,1);

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
models{end+1} = 'FEM Orig';


%% Generate the denser meshing

% Add a concentric mesh into the first compartment 4mm inwards
bnds(1) = go_add_concentric_mesh(bnds(1),0.004);
vol_fem = go_paperbnds2mesh_denserspine(bnds,2);

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
models{end+1} = 'FEM Denser';

%% Load the central sources and see how the eigenfields look

tmp = load(fullfile(files.results,'LFM_all_central_sources.mat'));


L{3} = tmp.L{10};


nmodels = numel(L);
nsources = size(L{1},2) / 3;

ratios = zeros(nmodels,nsources);

for ii = 1:nmodels
    for jj = 1:nsources

        sid = 3*(jj - 1) + (1:3);
        l = L{ii}(:,sid);

        [U,S,V] = svd(l,'econ');

        ratios(ii,jj) = S(1) / S(5);
        V1(ii,jj,:) = V(:,1);

    end
end


figure
plot(ratios(1:end,:)','-','linewidth',2)
grid on
ylim([0 ceil(1.1 * max(ratios(:)))])
ylabel('Ratio of S1 to S2')
xlabel('Location ID')
xlim([1, 60])
set(gcf,'color','w')

axis square

cmap = [0 0 1; 1 0 0; 0 0 0];
set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)
legend({'Original FEM','Denser FEM','5c LGi'},'location','eo')

fname = fullfile(files.results,'fem_meshing_results.png');
exportgraphics(gcf,fname,'Resolution',600);