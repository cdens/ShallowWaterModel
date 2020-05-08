% FinalProjectDriver.m
% Casey Densmore
% Driving code to run final project (shallow water equation model) for
% 12.850

clearvars
close all
clc


%% Initializing variables

dx = 20000; %dx = 20 km
dy = 20000; %dy = 20 km
% dt = 21600; %dt = 0.25 days (6 hours)
dt = 2000; %2000 seconds- satisfies wave CFL condition for long Rossby waves
% dt = 60; % 1 minute

nx = 200; %4000 km domain w/ dx = 20 km
ny = 100; %2000 km domain w/ dy = 20 km
nt = 1000; %200 days time domain @ 6 hrly, 23 days for 2000 sec

%CFL parameter
cfl = 7.27E-5*cosd(30)*nx^2*dx*dt/(2*pi^2*6536000);

minresidual = 1E-8;
maxiterations = 1000;

Ftau = 0; %default no forcing anywhere
kx = 0; %no diffusion
ky = 0; %no diffusion

caxisval = NaN; %default- auto determine color axis for gifs

ctrlat = 30; %run model at 30 degrees N
planetype = 'real'; %beta changes with latitude
xbound = 'periodic'; %periodic domain in x
gridtype = 'D'; %either C or D grid (determines which faces velocities are on)

friction = 1/(3*24*3600); % 3 days

% run base test, build x/y arrays
zeta_initial_empty = zeros(nx,ny); %no flow
[~,~,~,~,x,y,~] = runShallowWaterModel(zeta_initial_empty,nx,dx,ny,dy,3,dt,...
    minresidual,maxiterations,Ftau,kx,ky,ctrlat,planetype,xbound,gridtype,friction);



%% building initial conditions for Rossby wave

%generating psi
[yy,xx] = meshgrid(y,x);
xctr = 3250000; %starting x in m
yctr = dy*ny/2; %starting y = halfway up domain
r = 50000; %radius of wave in m
rr = sqrt((xx-xctr).^2 + (yy-yctr).^2);

%creating psi initial
zetamax = -5E-6; 
zeta_initial_rossby = zeros(nx,ny); %no flow
zeta_initial_rossby(rr <= r) = cos((pi/(2*r)).*rr(rr <= r)); %cosine tapered to 0 at distance r
zeta_initial_rossby = zetamax.*zeta_initial_rossby./max(zeta_initial_rossby(:));

%% f-plane (with or without friction)

[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_rossby,nx,dx,ny,dy,nt,dt,... 
    minresidual,maxiterations,Ftau,kx,ky,ctrlat,'f','periodic','a',friction);

makemodelplots('Output/fplane_withfriction_psi',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\Psi (m^2/s)')

[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_rossby,nx,dx,ny,dy,nt,dt,... 
    minresidual,maxiterations,Ftau,kx,ky,ctrlat,'f','periodic','a',0);

makemodelplots('Output/fplane_frictionless_psi_test',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\Psi (m^2/s)') 


%% beta-plane

%beta plane- periodic and no friction
[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_rossby,nx,dx,ny,dy,nt,dt,... 
    minresidual,maxiterations,Ftau,kx,ky,ctrlat,'beta','periodic','a',0);

makemodelplots('Output/betaplane_periodic_psi',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\psi (m^2/s)') 
makemodelgif('Output/betaplane_periodic_zeta',zeta,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\zeta (s^{-1})',0.05,caxisval)
makemodelgif('Output/betaplane_periodic_psi',psi,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\psi (m s^{-1})',0.05,caxisval)


%beta plane at 60N- periodic and no friction
[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_rossby,nx,dx,ny,dy,nt,dt,... 
    minresidual,maxiterations,Ftau,kx,ky,60,'beta','periodic','a',0);

makemodelplots('Output/betaplane60_periodic_psi',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\psi (m^2/s)') 
% makemodelgif('Output/betaplane60_periodic_zeta',zeta,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\zeta (s^{-1})',0.05,caxisval)
makemodelgif('Output/betaplane60_periodic_psi',psi,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\psi (m s^{-1})',0.05,[0,3500])


%beta plane at 0N- periodic and no friction
[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_rossby,nx,dx,ny,dy,nt,dt,... 
    minresidual,maxiterations,Ftau,kx,ky,0,'beta','periodic','a',0);

makemodelplots('Output/betaplane0_periodic_psi',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\psi (m^2/s)') 
% makemodelgif('Output/betaplane0_periodic_zeta',zeta,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\zeta (s^{-1})',0.05,caxisval)
% makemodelgif('Output/betaplane0_periodic_psi',psi,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\psi (m s^{-1})',0.05,caxisval)


%beta plane- wall and friction
[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_rossby,nx,dx,ny,dy,nt,dt,... 
    minresidual,maxiterations,Ftau,kx,ky,ctrlat,'beta','wall','a',friction);

makemodelplots('Output/betaplane_wall_psi',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\psi (m^2/s)') 
makemodelgif('Output/betaplane_wall_zeta',zeta,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\zeta (s^{-1})',0.05,caxisval)
makemodelgif('Output/betaplane_wall_psi',psi,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\psi (m s^{-1})',0.05,caxisval)


%% adding forcing- basin spinup

Ftau_spinup = -2E-13;
kspinup = 10^1;
spinupfriction = friction;

nxspin = nx; %same number of points in x
nyspin = ny; %same number of points in y
dxspin = 2*dx; %make the basin twice as wide
dyspin = 2*dy; %make the basin twice as wide
ntspin = 1500; %run it out a little longer

[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_empty,nxspin,dxspin,nyspin,dyspin,ntspin,dt,... 
    minresidual,maxiterations,Ftau_spinup,kspinup,kspinup,ctrlat,'beta','wall','a',friction);

makemodelplots('Output/spinup_psi',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\psi (m^2/s)') 
makemodelplots('Output/spinup_zeta',zeta,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\zeta (s^{-1})') 
makemodelgif('Output/spinup_zeta',zeta,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\zeta (s^{-1})',0.05,caxisval)
makemodelgif('Output/spinup_psi',psi,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\psi (m s^{-1})',0.05,caxisval)
savemodelfinalframe('Output/spinup_psi',psi,U,V,x/1000,y/1000,t/(24*3600),'\psi (m s^{-1})',caxisval)
% makemodelgif('Output/spinup_zeta_noUlabel',zeta,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\zeta (s^{-1})',0.05,caxisval)
% makemodelgif('Output/spinup_psi_noUlabel',psi,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\psi (m s^{-1})',0.05,caxisval)



%% time-variant 

nxvary = nx; %same number of points in x
nyvary = ny; %same number of points in y
dxvary = 2*dx; %make the basin twice as wide
dyvary = 2*dy; %make the basin twice as wide
ntvary = 1000; %run it out a little longer

%setting up wind stress forcing- applying a sinusoidal forcing to a square
%region on the west side of the domain
Ftaumax = 2E-13; %amplitude forcing
Ftau_timevary = zeros(nx,ny,nt);
forcesquare = zeros(nx,ny); %x/y region for forcing
forcesquare(190:198,46:54) = 1; %vertically centered near eastern bound
tvarying = (1:nt)*dt/(24*3600); %time in days
omega_force = 2*pi/9; %period of forcing = 9 days
for i = 1:nt
    Ftau_timevary(:,:,i) = forcesquare.*Ftaumax.*sin(omega_force*i*(dt/(24*3600)));
end
    
%first case- 30 N (with matching forcing period)
[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_empty,nxvary,dxvary,nyvary,dyvary,ntvary,dt,... 
    minresidual,maxiterations,Ftau_timevary,kspinup,kspinup,ctrlat,'real','wall','a',friction);
makemodelplots('Output/timevarying_psi',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\psi (m^2/s)')
makemodelgif('Output/timevarying_zeta',zeta,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\zeta (s^{-1})',0.05,caxisval)
makemodelgif('Output/timevarying_psi',psi,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\psi (m s^{-1})',0.05,[-2,2].*1E2)


%wave frequency lower than forcing frequency- 80 N
[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_empty,nxvary,dxvary,nyvary,dyvary,ntvary,dt,... 
    minresidual,maxiterations,Ftau_timevary,kspinup,kspinup,80,'real','wall','a',friction);
makemodelplots('Output/timevarying80N_psi',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\psi (m^2/s)')
makemodelgif('Output/timevarying80N_zeta',zeta,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\zeta (s^{-1})',0.05,caxisval)
makemodelgif('Output/timevarying80N_psi',psi,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\psi (m s^{-1})',0.05,[-2,2].*1E2)


%wave frequency higher than forcing frequency- 0 N
[U,V,psi,zeta,x,y,t] = runShallowWaterModel(zeta_initial_empty,nxvary,dxvary,nyvary,dyvary,ntvary,dt,... 
    minresidual,maxiterations,Ftau_timevary,kspinup,kspinup,0,'real','wall','a',friction);
makemodelplots('Output/timevarying0N_psi',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\psi (m^2/s)')
makemodelgif('Output/timevarying0N_zeta',zeta,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\zeta (s^{-1})',0.05,caxisval)
makemodelgif('Output/timevarying0N_psi',psi,U,V,x/1000,y/1000,t/(24*3600),1:5:length(t),'\psi (m s^{-1})',0.05,[-2,2].*1E2)
