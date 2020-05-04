function makemodelplots(savedir,var,x,y,t,xhovind,yhovind,varname,strtype,timeunits,distunits,makegif,framestodo)
close all

%% vertical hovmoller
hovdata = squeeze(var(xhovind,:,:));

figure(); clf;
ax = gca;

pcolor(ax,t,y,hovdata); 
colormap(ax,jet); 
shading(ax, 'interp')

xlabel(ax,['Time (',timeunits,')'])
ylabel(ax,['Y (',distunits,')'])
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
ylabel(ax,['X (',distunits,')'])
box(ax,'on')
grid(ax,'on')
set(ax,'layer','top')
yc = colorbar(ax);
ylabel(yc,varname)
saveas(gcf,[savedir,'_horzhov'],'png');

% close all

%% animation + final snapshot

figure(); 
set(gcf,'position',[1,1,1.8,1].*get(gcf,'position'))

startvar = squeeze(var(:,:,1));
cbounds = [min(startvar(:)),max(startvar(:))];

if ~makegif
    framestodo = length(t);
end

c = 0;
for i = framestodo
    c = c + 1;
    
    clf
    ax = gca;
    
    axis equal

    pcolor(ax,x,y,squeeze(var(:,:,i))'); 
    colormap(ax,jet); 
    caxis(ax,cbounds)
    shading(ax, 'interp')

    xlabel(ax,['X (',distunits,')'])
    ylabel(ax,['Y (',distunits,')'])
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
