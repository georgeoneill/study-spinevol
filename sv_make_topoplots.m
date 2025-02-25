clearvars
close all
clc

proj_init

source_id = 30; % Approx. T9 source

triaxial = 1; % Turn on or off triaxial sensors

files = [];
files.root = proj_path;

files.results = fullfile(files.root,'results','3axis');
files.data = fullfile(files.results,'LFM_all_central_sources.mat');

load(files.data);

models = {'Inf','SS','CS','1c LC','3c LC','4c LC','5c LC','5c LCi',...
    '5c LG','5c LGi','FEM'};

%% preprare mesh, sensors etc for plotters

load(fullfile(files.root,'geometries','all_geometries_v2.mat'));
% make a mesh for the sensor plane

sens_xy = grad.coilpos(1:250,[1 3]);
tri = delaunay(sens_xy);

% go thorough each face, calculate edge length, if any are over 50 mm,
% eliminate those faces
bads =[];
for ii = 1:size(tri,1)
    M.vertices = grad.coilpos;
    M.faces = tri(ii,:);
    D = spm_mesh_distmtx(M);

    d = full(unique(D));
    if any(d > 0.048)
        bads = [bads ii];
    end
end
tri(bads,:) = [];

grad_plot.p = grad.coilpos(1:250,:);
grad_plot.e = tri;

%% Plot one model for now

for model_id = 1:numel(models) % 1:numel(models) % for all plots!

    va = [-180 0];

    slice_id = 3*(source_id-1)+(1:3);
    l = L{model_id}(:,slice_id);
    l3d = reshape(l,250,3,3);
    plotscale(1) = max(max(abs(l3d(:,:,1))));
    plotscale(2) = max(max(abs(l3d(:,:,2))));
    plotscale(3) = max(max(abs(l3d(:,:,3))));

    figure;
    hold on
    axis off
    plotcline = 1;

    % AP plots
    vec = l3d(:,1,2);
    spine_topoplot_wrapper(vec,grad_plot,mesh_torso,...
        [3 3 1],va,plotscale(2),plotcline)
    vec = l3d(:,2,2);
    spine_topoplot_wrapper(vec,grad_plot,mesh_torso,...
        [3 3 2],va,plotscale(2),plotcline)
    vec = l3d(:,3,2);
    spine_topoplot_wrapper(vec,grad_plot,mesh_torso,...
        [3 3 3],va,plotscale(2),plotcline)

    % % RL plots
    vec = l3d(:,1,1);
    spine_topoplot_wrapper(vec,grad_plot,mesh_torso,...
        [3 3 4],va,plotscale(1),plotcline)
    vec = l3d(:,2,1);
    spine_topoplot_wrapper(vec,grad_plot,mesh_torso,...
        [3 3 5],va,plotscale(1),plotcline)
    vec = l3d(:,3,1);
    spine_topoplot_wrapper(vec,grad_plot,mesh_torso,...
        [3 3 6],va,plotscale(1),plotcline)

    % % IS plots
    vec = l3d(:,1,3);
    spine_topoplot_wrapper(vec,grad_plot,mesh_torso,...
        [3 3 7],va,plotscale(3),plotcline)
    vec = l3d(:,2,3);
    spine_topoplot_wrapper(vec,grad_plot,mesh_torso,...
        [3 3 8],va,plotscale(3),plotcline)
    vec = l3d(:,3,3);
    spine_topoplot_wrapper(vec,grad_plot,mesh_torso,...
        [3 3 9],va,plotscale(3),plotcline)

    set(gcf,'position',[ 744.0000   75.4000  592.2000  909.600], ...
        'color','w')

    plotscale = plotscale*1e6;
    modname = models{model_id};
    modname = strrep(modname,' ','_');
    fname = deblank(sprintf('fieldmaps_3by3_%s_%.2e_%.2e_%.2e.png',modname,plotscale(2),plotscale(1),plotscale(3)));
    fname = fullfile(files.results,fname);
    exportgraphics(gcf,fname,'resolution',600)

end