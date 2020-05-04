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
% dt = 21600/2; %dt = 0.125 days (3 hours)
dt = 60; % 1 minute

nx = 200; %4000 km domain w/ dx = 20 km
ny = 100; %2000 km domain w/ dy = 20 km
nt = 1000; %200 days time domain @ 6 hrly, 100 days @ 3 hrly

theta = 0.5; %Crank-Nicholson semi-implicit spatial discretization
minresidual = 1E-8;
maxiterations = 1000;

Ftau = 0; %default no forcing anywhere
kx = 0; %no diffusion
ky = 0; %no diffusion

ctrlat = 30; %run model at 30 degrees N
planetype = 'real'; %beta changes with latitude
xbound = 'periodic'; %periodic domain in x
gridtype = 'D'; %either C or D grid (determines which faces velocities are on)

isnonlinear = 0;


%% Negative control run- no flow, no fluxes
zeta_initial = zeros(nx,ny); %no flow
% [~,~,psi,~,~,~,~,~] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,10,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,planetype,xbound,'c',isnonlinear);
% disp(['Test Run Grid C: AVG psi_f = ',num2str(nanmean(squeeze(psi(:,:,end)),'all')),...
%     ',  STD psi_f = ',num2str(nanstd(squeeze(psi(:,:,end)),0,'all'))])
[~,~,~,~,x,y,~,~] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,10,dt,...
    theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,planetype,xbound,'d',isnonlinear);
% disp(['Test Run Grid D: AVG psi_f = ',num2str(nanmean(squeeze(psi(:,:,end)),'all')),...
%     ',  STD psi_f = ',num2str(nanstd(squeeze(psi(:,:,end)),0,'all'))])


%% Basic model run: midlatitude Rossby wave (testing configurations)

%generating psi
[yy,xx] = meshgrid(y,x);
xctr = 3250000; %starting x = 1750 km
yctr = dy*ny/2; %starting y = halfway up domain
r = 500000; %radius of wave = 100 km
rr = sqrt((xx-xctr).^2 + (yy-yctr).^2);

%creating psi initial
zetamax = -5E-6; 
zeta_initial = zeros(nx,ny); %no flow
zeta_initial(rr <= r) = cos((pi/(2*r)).*rr(rr <= r)); %cosine tapered to 0 at distance r
zeta_initial = zetamax.*zeta_initial./max(zeta_initial(:));

%f-plane -ALL WORK
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'f','periodic','c',isnonlinear);
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'f','periodic','d',isnonlinear);
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'f','wall','c',isnonlinear);
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'f','wall','d',isnonlinear);

RossbySpeed = (2.*7.27E-5.*cosd(ctrlat)./6356000)*4*pi^2/(2*r*1E3)^2; %Rossby wave speed
TimeToPropagate = (nx*dx/RossbySpeed)/(24*3600); %time to propagate across dom ain in days

%beta plane
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'beta','periodic','a',isnonlinear);
[U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
    theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'beta','periodic','c',isnonlinear);
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     0,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'beta','periodic','c',isnonlinear); %forward Euler method
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     1,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'beta','periodic','c',isnonlinear); %backward Euler method
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'beta','periodic','d',isnonlinear);
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'beta','wall','c',isnonlinear);
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'beta','wall','d',isnonlinear);

%real
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'real','periodic','c',isnonlinear);
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'real','periodic','d',isnonlinear);
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'real','wall','c',isnonlinear);
% [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta_initial,nx,dx,ny,dy,500,dt,...
%     theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'real','wall','d',isnonlinear);

makemodelplots('temp',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\Psi (m^2/s)','%04.2f','days','km',true,1:501) %streamfunction
makemodelplots('vorticity',zeta,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\zeta (s^{-1})','%04.2f','days','km',true,1:501) %voriticity

%% Racing Rossby waves



%% Midlatitude Kelvin waves



%% Basin spinup



%% Equatorial Kelvin + Rossby (delayed oscillator)



%% Equatorial- noisy forcing


