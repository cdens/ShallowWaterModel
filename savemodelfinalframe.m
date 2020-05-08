function savemodelfinalframe(savedir,zeta,U,V,x,y,t,varlabel,caxisval)
close all


%% animation + final snapshot

figure();
set(gcf,'position',[1,1,1.8,1].*get(gcf,'position'))

if ~isnan(caxisval)
    cbounds = caxisval;
else
    cbounds = [min(zeta(:)),max(zeta(:))];
end

xbounds = [min(x),max(x)];
ybounds = [min(y),max(y)];

i = length(t);

clf
ax = gca;

axis equal

dvx = 10;
dvy = 5;

Us = squeeze(U(:,:,i));
Vs = squeeze(V(:,:,i));

vel = sqrt(Us.^2 + Vs.^2);
maxvel = max(vel(:));

hold on
pcolor(ax,x,y,squeeze(zeta(:,:,i))');
quiver(ax,x(1:dvx:end),y(1:dvy:end),Us(1:dvx:end,1:dvy:end)',Vs(1:dvx:end,1:dvy:end)','color','w')

colormap(ax,jet);
caxis(ax,cbounds)
shading(ax, 'interp')

xlim(xbounds)
ylim(ybounds)
xlabel(ax,'X (km)')
ylabel(ax,'Y (km)')
box(ax,'on')
grid(ax,'on')
set(ax,'layer','top')
yc = colorbar(ax);
ylabel(yc,varlabel)
caxis(ax,cbounds)

title(ax,['Time: ',sprintf('%04.2f',t(i)),' days, U_{MAX}=',sprintf('%04.2f',maxvel),'m/s'])

saveas(gcf,[savedir,'_finalframe'],'png')

