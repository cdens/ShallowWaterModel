function makemodelplots(savedir,var,x,y,t,xhovind,yhovind,varname,strtype,timeunits,makegif)
close all

%% vertical hovmoller
hovdata = squeeze(var(xhovind,:,:));

figure(); clf;
ax = gca;

pcolor(ax,t,y,hovdata); 
colormap(ax,jet); 
shading(ax, 'interp')

xlabel(ax,['Time (',timeunits,')'])
ylabel(ax,'Y (m)')
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

xlabel(ax,['Time (',timeunits,')'])
ylabel(ax,'X (m)')
box(ax,'on')
grid(ax,'on')
set(ax,'layer','top')
yc = colorbar(ax);
ylabel(yc,varname)
saveas(gcf,[savedir,'_horzhov'],'png');

% close all

%% animation + final snapshot

figure(); 
cbounds = [min(var(:)),max(var(:))];

if makegif
    framestodo = 1:5:length(t);
else
    framestodo = length(t);
end

c = 0;
for i = framestodo
    c = c + 1;
    
    clf
    ax = gca;

    pcolor(ax,x,y,squeeze(var(:,:,i))'); 
    colormap(ax,jet); 
    caxis(ax,cbounds)
    shading(ax, 'interp')

    xlabel(ax,'X (m)')
    ylabel(ax,'Y (m)')
    box(ax,'on')
    grid(ax,'on')
    set(ax,'layer','top')
    yc = colorbar(ax);
    ylabel(yc,varname)
    
    title(ax,['Time: ',sprintf(strtype,t(i)),' ',timeunits])
    
    if makegif
        M(c) = getframe(gcf);
    end
end

%saving last frame
saveas(gcf,[savedir,'_finalframe'],'png');

%saving gif
if makegif
    movie2gif(M,[savedir,'_animation.gif'],'DelayTime',0.03,'LoopCount',Inf)
end

% close all
