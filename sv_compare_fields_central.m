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

files.data = fullfile(files.results,'LFM_all_central_sources.mat');

load(files.data);

% Update the models names to match paper
models = {'Inf','SS','CS','1c LC','3c LC','4c LC','5c LC','5c LCi',...
    '5c LG','5c LGi','FEM'};

%% Perform the intital comparison - make plots
% Plot RE and CC matricies for all combined orientations

[re, cc] = compare_results(L);
export_fieldcomp_figures(models,re,cc,files.results,'all');

%% Repeat but for differing orientations

% Right-Left

ori = repmat([1 0 0],size(L{1},2)/3,1);

Lrl = [];
for ii = 1:numel(L)
    tmp = [];
    lmodel = L{ii};
    lmodel = reshape(lmodel,size(lmodel,1),3,[]);
    for jj = 1:length(ori)
        tmp(jj,:) = lmodel(:,:,jj) * ori(jj,:)';
    end
    Lrl{ii} = tmp;
end

[re, cc] = compare_results(Lrl);
export_fieldcomp_figures(models,re,cc,files.results,'rl');


% Inferior-Superior

ori = repmat([0 0 1],size(L{1},2)/3,1);

Lis = [];
for ii = 1:numel(L)
    tmp = [];
    lmodel = L{ii};
    lmodel = reshape(lmodel,size(lmodel,1),3,[]);
    for jj = 1:length(ori)
        tmp(jj,:) = lmodel(:,:,jj) * ori(jj,:)';
    end
    Lis{ii} = tmp;
end

[re, cc] = compare_results(Lis);
export_fieldcomp_figures(models,re,cc,files.results,'is');

% Aterior-Posterior

ori = repmat([0 1 0],size(L{1},2)/3,1);

Lap = [];
for ii = 1:numel(L)
    tmp = [];
    lmodel = L{ii};
    lmodel = reshape(lmodel,size(lmodel,1),3,[]);
    for jj = 1:length(ori)
        tmp(jj,:) = lmodel(:,:,jj) * ori(jj,:)';
    end
    Lap{ii} = tmp;
end

[re, cc] = compare_results(Lap);
export_fieldcomp_figures(models,re,cc,files.results,'ap');

% Spine curvature
pos_line = src.pos(1:end,:);
ori = diff(pos_line);
ori = ori ./vnorm(ori,2);
ori(end+1,:) = ori(end,:);

Lcurv = [];
for ii = 1:numel(L)
    tmp = [];
    lmodel = L{ii};
    lmodel = reshape(lmodel,size(lmodel,1),3,[]);
    for jj = 1:length(ori)
        tmp(jj,:) = lmodel(:,:,jj) * ori(jj,:)';
    end
    Lcurv{ii} = tmp;
end

[re, cc] = compare_results(Lcurv);
export_fieldcomp_figures(models,re,cc,files.results,'spinecurv');
