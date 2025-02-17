function export_fieldcomp_figures(models,re,cc,resdir,tag)

files.results = resdir;

figure;clf
imagesc(re)
axis equal
axis off
% colorbar
colormap(brewermap(100,'RdPu'))
clim([0 1]);

for ii = 1:numel(re)
    [x,y] = ind2sub(size(re),ii);
    if x ~= y

        if abs(re(ii)) < 0.52
            col = 'k';
        else
            col = 'w';
        end
        text(y,x,sprintf('%.2f',re(ii)),...
            'horizontalalignment','center','color',col,...
            'fontsize',12,'FontName',proj_font)
    elseif re(ii) > 0
        text(y,x,sprintf('%1.e',re(ii)),...
            'horizontalalignment','center','color','w',...
            'fontsize',12,'FontName',proj_font)
    end


end
set(gcf,'color','w')
set(gcf,'position',[   481   411   879   527]);
fname = fullfile(files.results,['re_mat_' tag '.png']);
exportgraphics(gcf,fname,'Resolution',600);

% Generate Correlation Matrix

figure;clf
imagesc(cc)
axis equal
axis off
% colorbar
colormap(brewermap(100,'Reds'))
caxis([0 1])

for ii = 1:numel(cc)

    [x,y] = ind2sub(size(cc),ii);
    if x ~= y
        if abs(cc(ii)) > 0.58;
            col = 'w';
        else
            col = 'k';
        end
        text(y,x,sprintf('%.2f',cc(ii)),...
            'horizontalalignment','center','color',col,...
            'fontsize',12,'FontName',proj_font)

    end


end

set(gcf,'color','w')
set(gcf,'position',[481   411   879   527]);
fname = fullfile(files.results,['cc_mat_' tag '.png']);
exportgraphics(gcf,fname,'Resolution',600);

% Plot Correlation, Errors relative to FEM

figure;clf
plot(1:11,re(1:11,end),'linewidth',2)
hold on
plot(1:11,cc(1:11,end),'linewidth',2)
ylabel('Metric')
xticklabels(models)
axis square

set(gcf,'color','w')
grid on

set(gcf,'Position',[ 473.8000  180.2000  652.0000  474.4000])
xlim([1 11])
ylim([0 1])
xline(6.5,'--','','linewidth',2)
% legend('Relative Error','Correlation^2','Location','eo')

cmap = [28 73 136;
    205 0 0]/255;

set(gca,'colororder',cmap,'fontsize',14,'FontName',proj_font)
fname = fullfile(files.results,...
    ['comp2fem_' tag '.png']);
exportgraphics(gca,fname,'Resolution',600)