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


%% 5 Tissue BEM

L = [];
exec_time = [];
models = [];

cratio = [1 2 3 4 5:5:40];


for jj = 1:numel(cratio)

ti = tic;


ci = [0.33 0.33/cratio(jj) .62 .05 .23];
co = [0.33/cratio(jj) .23 .23 .23 0];

cfg                     = [];
cfg.method              = 'bem_hbf';
cfg.conductivity        = [ci;co];
vol_bem5c               = ft_prepare_headmodel(cfg,bnds);

cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol_bem5c;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{jj} = cell2mat({fwd_tmp.leadfield{:}});
models{jj} = sprintf('cratio  = %02d',cratio(jj));

end

%% Compare fields

[re cc] = compare_results(L);

figure
plot(cratio,re(end,:),'linewidth',2);
hold on
plot(cratio,cc(end,:),'linewidth',2);
ylabel('Metric')
xlabel('Conductivity Ratio')
axis square
set(gcf,'color','w')
grid on
set(gcf,'Position',[ 315   234   727   475])
xlim([1 40])
ylim([0 1])
legend('Relative Error','Correlation^2','Location','eo')
cmap = [28 73 136;
    205 0 0]/255;
set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)
fname = fullfile(files.results,'cond_adj_metrics.png');
exportgraphics(gca,fname,'resolution',600);

%% Look at ratio of eigenfields along spine

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

for ii = 1:numel(cratio)
    ccell{ii} = num2str(cratio(ii));
end

cmap = viridis(12);
set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)
legend(ccell,'location','eo')
set(gcf,'Position',[ 473.8000  180.2000  652.0000  474.4000])
fname = fullfile(files.results,'cond_adj_eigenfields.png');
exportgraphics(gca,fname,'resolution',600);

%% Plot something here

figure
scatter(cratio,ratios(:,30),'o','filled');
xlabel('Condtivity Ratio')
ylabel('Ratio of S1 to S2')
axis square
grid on
set(gcf,'Position',[ 473.8000  180.2000  652.0000  474.4000],'color','w')
set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)

