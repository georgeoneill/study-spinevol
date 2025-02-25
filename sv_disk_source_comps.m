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

files.data = fullfile(files.results,'LFM_all_disk_sources.mat');
load(files.data);

% Update the models names to match paper
models = {'Inf','SS','CS','1c LC','3c LC','4c LC','5c LC','5c LCi',...
    '5c LG','5c LGi','FEM'};

%% Determine current flow along the orientation of the spine

% get difference in positions along spine and normalise
pos_line = src.pos(1:8:end,:);
ori = diff(pos_line);
ori = ori ./vnorm(ori,2);
ori(end+1,:) = ori(end,:);

%% Get the stats between pairs of dipoles

clear re cc

labels_disk = {'pc','pl','cl','al','ac','ar','cr','pr'};
disk_id = [3 7]; % 3-> central-left; 7-> central-right

for ii = 1:numel(L)

    Lmodel = L{ii};
    Lmodel = reshape(Lmodel,size(Lmodel,1),3,[]);
    Lset1 = Lmodel(:,:,disk_id(1):8:end);
    Lset2 = Lmodel(:,:,disk_id(2):8:end);

    for jj = 1:size(Lset1,3)

        tmp1 = Lset1(:,:,jj);
        tmp2 = Lset2(:,:,jj);

        % Reduce down to dipoles oriented along spine for clarity
        tmp1a = tmp1*ori(jj,:)';
        tmp2a = tmp2*ori(jj,:)';

        cc(jj,ii) = corr(tmp1a(:),tmp2a(:)).^2;
        re(jj,ii) = relative_error(tmp1a(:),tmp2a(:));

    end

end

%% Plot these using the Gramm toolbox

mods = [];

for ii = 1:numel(L)
    [mods{((ii-1)*61)+(1:61)}] = deal(models{ii});
end

for type = {'re','cc'}

    data = eval(type{:});

    figure
    g = gramm('x',mods,'y',data(:),'color',mods);
    g.geom_jitter('width',0.5,'edgewidth',0.1,'edgecolor','none'); %set edgecolor to 'none' to remove any outlines of points
    g.stat_boxplot('width',5,'alpha',0,'dodge',1,'linewidth',1.5,'drawoutlier',0);
    g.set_order_options('x',0,'color',0);

    switch type{:}
        case 'cc'
            g.set_names('x',[],'y','Correlation^2');
        case 're'
            g.set_names('x',[],'y','Relative Error');
    end

    g.set_color_options("map",turbo(11));

    switch type{:}
        case 'cc'
            g.axe_property('xgrid','on','ygrid','on','xticklabelrotation',30,'ylim',[0 1]);
        case 're'
            g.axe_property('xgrid','on','ygrid','on','xticklabelrotation',30,'ylim',[0 1]);
    end

    g.set_text_options('base_size',14,'font','atkinson hyperlegible');
    set(gcf,'position',[744.0000  430.6000  582.6000  519.4000]);
    g.no_legend();
    g.draw();

        fname = sprintf('diff_fields_allmodels_%s.png',type{:});
    
    fname = fullfile(files.results,fname);
    exportgraphics(gcf,fname,'Resolution',600);

end

%% Generate the difference plots between two models

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

model_ids = [4 7];
va = [-180 0];

src_id = 30; % its always 30!

for ii = 1:numel(model_ids)

    Lmodel = L{model_ids(ii)};
    Lmodel = reshape(Lmodel,size(Lmodel,1),3,[]);

    Lset1 = Lmodel(:,:,disk_id(1):8:end);
    Lset2 = Lmodel(:,:,disk_id(2):8:end);

    % Just take the radial sensors, and the specific source orientation
    L1 = Lset1(251:500,:,src_id)*ori(src_id,:)';
    L2 = Lset2(251:500,:,src_id)*ori(src_id,:)';
    dL = L2 - L1;

    plotscale = max(max([abs(L1) abs(L2) dL(:)]));

    figure
    axis off
    hold on
    spine_topoplot_wrapper(L1,grad_plot,mesh_torso,...
        [1 3 1],va,plotscale,1)
    spine_topoplot_wrapper(L2,grad_plot,mesh_torso,...
        [1 3 2],va,plotscale,1)
    spine_topoplot_wrapper(dL,grad_plot,mesh_torso,...
        [1 3 3],va,plotscale,1)

    set(gcf,'color','w')

    fname = sprintf('diff_fields_%s_%.2e.png',models{model_ids(ii)},...
        plotscale*1e6);
    
    fname = fullfile(files.results,fname);
    exportgraphics(gcf,fname,'Resolution',600);

end

%% Auxillary function
function e = relative_error(La,Lb)
e = norm(Lb-La)./(norm(La)+norm(Lb));
end