function makemodelplots(savedir,var,x,y,t,xhovind,yhovind,varname,timeunits,distunits)
close all

%% vertical hovmoller
hovdata = squeeze(var(xhovind,:,:));

figure(); clf;
ax = gca;

pcolor(ax,t,y,hovdata); 
colormap(ax,jet); 
shading(ax, 'interp')

xlabel(ax,'Time (days)')
ylabel(ax,'Y (km)')
box(ax,'on')
grid(ax,'on')
set(ax,'layer','top')
yc = colorbar(ax);
ylabel(yc,varname)
saveas(gcf,[savedir,'_verthov'],'png');

% close all

%% horizontal hovmoller
hovdata = squeeze(var(:,yhovind,:));

figure(); clf;
ax = gca;

pcolor(ax,t,x,hovdata); 
colormap(ax,jet); 
shading(ax, 'interp')

xlabel(ax,'Time (days)')
ylabel(ax,'X (km)')
box(ax,'on')
grid(ax,'on')
set(ax,'layer','top')
yc = colorbar(ax);
ylabel(yc,varname)
saveas(gcf,[savedir,'_horzhov'],'png');

