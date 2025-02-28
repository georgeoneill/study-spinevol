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

%% For each model, perform SVD on field triplets to see ratio

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
plot(ratios(1:6,:)','-','linewidth',2)
hold on
plot(ratios(7:end,:)','-','linewidth',2)
grid on
ylim([0 ceil(1.1 * max(ratios(:)))])
ylabel('Ratio of S1 to S2')
xlabel('Location ID')
xlim([1, 60])
set(gcf,'color','w')

axis square

cmap = turbo(11);
set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)
legend(models,'location','eo')
fname = fullfile(files.results,'field_eigenvalue_ratios.png');
exportgraphics(gcf,fname,'Resolution',600);

%% Determine orientation error
% for now see how much the orientation error is detemined by the basic
% cardinal orientation of the world or perhaps the spine direction, so
% first we need to work out the general flow of the spine

result = 'si'; %options 'spine', 'si', 'lr'

dpos = diff(src.pos);
spine_ori = dpos./ vnorm(dpos,2);
nlocs = length(src.pos);
is_ori = repmat([0 0 1],nlocs-1,1);
lr_ori = repmat([1 0 0],nlocs-1,1);


for ii = 1:numel(L)
    absv = squeeze(V1(ii,1:nlocs-1,:));
    is_err(ii,:) = acosd(dot(absv,is_ori,2));
    lr_err(ii,:) = acosd(dot(absv,lr_ori,2));
    spine_err(ii,:) = acosd(dot(absv,spine_ori,2));
end

is_err = min(is_err,abs(180-is_err));
lr_err = min(lr_err,abs(180-lr_err));

spine_err = min(spine_err,abs(180-spine_err));


btwn_err = acosd(dot(absv,is_ori,2));
btwn_err = min(btwn_err,abs(180-btwn_err));


figure(100);clf

switch result
    case 'lr'
        plot(lr_err','linewidth',2)
        % title('Angle between V1 and LR')
        ylabel('Angle between V1 and RL axis / degrees')

    case 'si'
        plot(is_err','linewidth',2)
        % title('Angle between V1 and SI')
        ylabel('Angle between V1 and IS axis / degrees')

    case 'spine'
        plot(spine_err','linewidth',2)
        title('Angle between V1 and Spine')

end

% legend(models,'Location','eo')
set(gcf,'color','w')
xlabel('Location ID')
xlim([1, 60])
ylim([0 90])
grid on
hold on
cmap = turbo(11);
set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)
switch result
    case 'si'
        plot(btwn_err,'m--','linewidth',3)
end
axis square
legend(models,'location','eo')
fname = fullfile(files.results,'field_eigenvector_orientations.png');
exportgraphics(gcf,fname,'Resolution',600);

%% Plot an eigenfield or two?

% prep meshes

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

source_no = 30;

for model_no = [4 7]

    for ii = 1:model_no
        for jj = source_no

            sid = 3*(jj - 1) + (1:3);
            l = L{ii}(:,sid);

            [U,S,V] = svd(l,'econ');

            % ratios(ii,jj) = S(1) / S(5);
            % V1(ii,jj,:) = V(:,1);

        end
    end


    Lrot = U * S * eye(3);
    ca = max(max(abs(Lrot(251:500,:))));


    figure
    axis off
    hold on

    ps = 3;
    va = [180 0];

    % cmap = flip(brewermap(35,'RdBu'));

    spine_topoplot_wrapper(Lrot(251:500,1),grad_plot,mesh_torso,...
        [1 2 1],va,ca,1)

    spine_topoplot_wrapper(Lrot(251:500,2),grad_plot,mesh_torso,...
        [1 2 2],va,ca,1)

    % Just want to be able to adequately highlight the contours for publication
    if model_no == 7
        ca2 = ca;
        ca = max(abs(Lrot(251:500,2)));
    else
        ca2 = ca;
    end
    % clines = linspace(-ca,ca,12);
    % ft_plot_topo3d_v2(grad_plot.p,grad_plot.e,Lrot(251:500,2),...
    %     'isolines',clines,'contourstyle','black',...
    %     'neighbourdist',inf,'topostyle',false)


    set(gcf,'color','w')

    fname = fullfile(files.results,sprintf('eigenfields_%s_%.2e_ft.png',models{model_no},ca2*1e6));
    exportgraphics(gcf,fname,'Resolution',600)

end