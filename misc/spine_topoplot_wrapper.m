function spine_topoplot_wrapper(map,sens,torso,ax,va,plotscale,plotcline)

if isempty(which('CreateAxes'))
    error('You do not have the plotting tools due to them not being licensed for sharing.')
end

cmap = flip(brewermap(25,'RdBu'));

CreateAxes(ax(1),ax(2),ax(3)); hold on;
PlotMesh(torso, 'facecolor', [.8 .8 .8],'facealpha',.2,'view',va);
PlotDataOnMesh(sens, map,'caxis',[-plotscale plotscale],...
    'view',va,'colorbar',0,'colormap',cmap,'pointset',...
   sens.p+[0 0.01 0],'pointsize',4);

% lim = max(abs(vec));
clines = linspace(-plotscale,plotscale,12);

if plotcline
    try
    ft_plot_topo3d_v2(sens.p,sens.e,map,...
    'isolines',clines,'contourstyle','black',...
    'neighbourdist',inf,'topostyle',false)
    end
end