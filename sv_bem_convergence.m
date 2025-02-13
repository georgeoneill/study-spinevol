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
src.pos = sources_disk.p;
src.inside = ones(length(src.pos),1);
src.unit = 'm';

%% Convergence testing

% To see if the 1-shell BEM is adequately sampled, we shall reduce the mesh
% density of the torso to determine how much the lead fields change -
% hopefully we see that the similarity between patterns is very high at the
% mesh density

ref_bnd.vertices = bnds(5).pos;
ref_bnd.faces = bnds(5).tri;

ntri = length(ref_bnd.faces);

targets = floor(linspace(0.05*ntri,ntri,20));

L1c = cell(numel(targets),1);

parfor ii = 1:numel(targets)

    ds_mesh = spm_mesh_reduce(ref_bnd,targets(ii));

    bnd = [];
    bnd.tri         = ds_mesh.faces;
    bnd.pos         = ds_mesh.vertices;
    bnd.unit        = 'm';
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

    L1c{ii} = cell2mat({fwd_tmp.leadfield{:}});

end

[re cc] = compare_results(L1c);

figure;clf
plot(5:5:100,re(:,end),'linewidth',2)
hold on
plot(5:5:100,cc(:,end),'linewidth',2)
ylabel('Metric')
xlabel('Torso Mesh Density / %')
axis square

set(gcf,'color','w')
grid on

set(gcf,'Position',[ 473.8000  180.2000  652.0000  474.4000])
xlim([5 100])
legend('Relative Error','Correlation^2','Location','eo')

cmap = [28 73 136;
    205 0 0]/255;

set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)

fname = fullfile(files.results,'convergence_1cBEM.png');
exportgraphics(gcf,fname,'resolution',600)


%% Generate the initial field patterns, use as reference

cratio = 40;

bnd = [];

bnd(1).pos = bnds(1).pos;
bnd(1).tri = bnds(1).tri;
bnd(1).unit = 'm';
bnd(2) = bnds(2);
bnd(3) = bnds(5);

ci = [0.33 0.33/cratio 0.23];
co = [0.33/cratio 0.23 0];

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

Lref = cell2mat({fwd_tmp.leadfield{:}});


%% Now adjust the WM and leave the bone intact

cratio = 40;

% Get white matter meshes and target densities
ref_wm.vertices = bnds(1).pos;
ref_wm.faces = bnds(1).tri;
ntri_wm = length(ref_wm.faces);
targets_wm = floor(linspace(0.05*ntri_wm,ntri_wm,20));

% Get bone meshes and target densities
ref_bone.vertices = bnds(2).pos;
ref_bone.faces = bnds(2).tri;
ntri_bone = length(ref_bone.faces);
targets_bone = floor(linspace(0.05*ntri_bone,ntri_bone,20));

re = cell(20,20);
cc = cell(20,20);

files.checkpoints = fullfile(files.results,'convergence_results_5c.mat');
if ~exist(files.checkpoints,'file')
    save(files.checkpoints,'cc','re')
else
    load(files.checkpoints)
end

for ii = 1:numel(targets_wm)
    for jj = 1:numel(targets_bone)
        if isempty(re{ii,jj})

        ds_wm = spm_mesh_reduce(ref_wm,targets_wm(ii));

        bnd = [];

        bnd(1).tri           = ds_wm.faces;
        bnd(1).pos           = ds_wm.vertices;
        bnd(1).unit           =  'm';

        ds_bone = spm_mesh_reduce(ref_bone,targets_bone(jj));

        bnd(2).tri = ds_bone.faces;
        bnd(2).pos = ds_bone.vertices;
        bnd(2).unit           =  'm'; 

        bnd(3) = bnds(5);

        ci = [0.33 0.33/cratio 0.23];
        co = [0.33/cratio 0.23 0];

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

        Lds = cell2mat({fwd_tmp.leadfield{:}});

        tmp = [];
        tmp{1} = Lds;
        tmp{2} = Lref;

        [re_tmp, cc_tmp] = compare_results(tmp);

        re{ii,jj} = re_tmp(2,1);
        cc{ii,jj} = cc_tmp(2,1);

        save(files.checkpoints,'cc','re')

        end
    end
end

%% Generate the plots for above
close all

if iscell(re)
re = cell2mat(re);
end

figure
imagesc(re)
axis equal
axis off
% colorbar
colormap(brewermap(100,'RdPu'))
clim([0 1]);
set(gcf,'color','w')
set(gcf,'position',[   481   411   879   527]);
if iscell(cc)
cc = cell2mat(cc);
end
fname = fullfile(files.results,'convergence_5c_re_mat.png');
exportgraphics(gca,fname,'resolution',600);

figure
imagesc(cc)
axis equal
axis off
% colorbar
colormap(brewermap(100,'Reds'))
clim([0 1]);
set(gcf,'color','w')
set(gcf,'position',[   481   411   879   527]);
fname = fullfile(files.results,'convergence_5c_cc_mat.png');
exportgraphics(gca,fname,'resolution',600);

% Plot metrics for full density wm mesh, adjusted bone mesh
figure
plot(5:5:100,re(end,:),'linewidth',2);
hold on
plot(5:5:100,cc(end,:),'linewidth',2);
ylabel('Metric')
xlabel('Bone Mesh Density / %')
axis square
set(gcf,'color','w')
grid on
set(gcf,'Position',[ 473.8000  180.2000  652.0000  474.4000])
xlim([5 100])
legend('Relative Error','Correlation^2','Location','eo')
cmap = [28 73 136;
    205 0 0]/255;
set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)
fname = fullfile(files.results,'convergence_5c_adjust_bone.png');
exportgraphics(gca,fname,'resolution',600);

% Plot metrics for full density bone mesh, adjusted wm mesh
figure
plot(5:5:100,re(:,end),'linewidth',2);
hold on
plot(5:5:100,cc(:,end),'linewidth',2);
ylabel('Metric')
xlabel('Spinal Cord Mesh Density / %')
axis square
set(gcf,'color','w')
grid on
set(gcf,'Position',[ 473.8000  180.2000  652.0000  474.4000])
xlim([5 100])
legend('Relative Error','Correlation^2','Location','eo')
cmap = [28 73 136;
    205 0 0]/255;
set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)
fname = fullfile(files.results,'convergence_5c_adjust_sc.png');
exportgraphics(gca,fname,'resolution',600);

figure
imagesc(rand(size(re)))
axis equal
axis off
% colorbar
colormap(brewermap(100,'RdPu'))
clim([0 1]);
set(gcf,'color','w')
set(gcf,'position',[   481   411   879   527]);
if iscell(cc)
cc = cell2mat(cc);
end
fname = fullfile(files.results,'random_alignment_mat.png');
exportgraphics(gca,fname,'resolution',600);