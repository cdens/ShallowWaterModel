% FinalProjectDriver.m
% Casey Densmore
% Driving code to run final project (shallow water equation model) for
% 12.850

clearvars
close all
clc


%% Initializing variables

dx = 5000; %dx = 5 km
dy = 5000; %dy = 5 km
dt = 3600; %dt = 1 hour

nx = 400; %2000 km domain w/ dx = 5 km
ny = 200; %1000 km domain w/ dy = 5 km
nt = 1000; %1,000 hours time domain

theta = 0.5; %Crank-Nicholson semi-implicit spatial discretization
minresidual = 1E-8;
maxiterations = 1000;

Ftau = 0; %default no forcing anywhere
kx = 0; %no diffusion
ky = 0; %no diffusion

ctrlat = 30; %run model at 30 degrees N
planetype = 'real'; %beta changes with latitude
xbound = 'periodic'; %periodic domain in x


%% Negative control run- no flow, no fluxes
psi_initial = zeros(nx,ny); %no flow
[~,~,psi,~,x,y,~,~] = runShallowWaterModel(psi_initial,nx,dx,ny,dy,10,dt,...
    theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,planetype,xbound);
disp(['Test Run: AVG psi_f = ',num2str(nanmean(squeeze(psi(:,:,end)),'all')),...
    ',  STD psi_f = ',num2str(nanstd(squeeze(psi(:,:,end)),0,'all'))])


%% Basic model run: midlatitude Rossby wave

%generating psi
[yy,xx] = meshgrid(y,x);
xctr = 1750000; %starting x = 1750 km
yctr = dy*ny/2; %starting y = halfway up domain
r = 100000; %radius of wave = 100 km
rr = sqrt((xx-xctr).^2 + (yy-yctr).^2);

%creating psi- maximum value corresponds to an increase in velocity
%from 0 to 1 m/s over radius
psimax = 1/r; 
psi_initial = zeros(nx,ny); %no flow
psi_initial(rr <= r) = psimax.*cos((pi/(2*r)).*rr(rr <= r)); %cosine tapered to 0 at distance r

%running and making output plots
[U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(psi_initial,nx,dx,ny,dy,50,dt,...
    theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,planetype,xbound);
makemodelplots('../Figs/InitialRossby',psi,x,y,t,nx/2,ny/2,'\Psi (m^2/s)','%04.2f','days',true)


%% Racing Rossby waves



%% Midlatitude Kelvin waves



%% Basin spinup



%% Equatorial Kelvin + Rossby (delayed oscillator)



%% Equatorial- noisy forcing


