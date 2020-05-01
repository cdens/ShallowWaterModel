% FinalProjectDriver.m
% Casey Densmore
% Driving code to run final project (shallow water equation model) for
% 12.850

clearvars
close all
clc


%% Initializing variables

dx = 10000; %dx = 10 km
dy = 10000; %dy = 10 km
dt = 1800; %dt = 30 mins

nx = 200; %2000 km domain w/ dx = 10 km
ny = 100; %1000 km domain w/ dy = 10 km
nt = 1000; %500 hours time domain

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


%% Negative control run- no flow, no fluxes
psi_initial = zeros(nx,ny); %no flow
[~,~,psi,~,~,~,~,~] = runShallowWaterModel(psi_initial,nx,dx,ny,dy,10,dt,...
    theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,planetype,xbound,'c');
disp(['Test Run Grid C: AVG psi_f = ',num2str(nanmean(squeeze(psi(:,:,end)),'all')),...
    ',  STD psi_f = ',num2str(nanstd(squeeze(psi(:,:,end)),0,'all'))])
[~,~,psi,~,x,y,~,~] = runShallowWaterModel(psi_initial,nx,dx,ny,dy,10,dt,...
    theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,planetype,xbound,'d');
disp(['Test Run Grid D: AVG psi_f = ',num2str(nanmean(squeeze(psi(:,:,end)),'all')),...
    ',  STD psi_f = ',num2str(nanstd(squeeze(psi(:,:,end)),0,'all'))])


%% Basic model run: midlatitude Rossby wave

%generating psi
[yy,xx] = meshgrid(y,x);
xctr = 1750000; %starting x = 1750 km
yctr = dy*ny/2; %starting y = halfway up domain
r = 100000; %radius of wave = 100 km
rr = sqrt((xx-xctr).^2 + (yy-yctr).^2);

%creating psi initial
psimax = 10000; 
psi_initial = zeros(nx,ny); %no flow
psi_initial(rr <= r) = cos((pi/(2*r)).*rr(rr <= r)); %cosine tapered to 0 at distance r
psi_initial = psimax.*psi_initial./max(psi_initial(:));
[U,V] = getvelocityfrompsi(psi_initial,dx,dy,xbound,'c');
% % % Uctr = 0.5.*(U(:,1:end-1) + U(:,2:end)); %D grid
% % % Vctr = 0.5.*(V(1:end-1,:) + V(2:end,:));
% % % % Vctr = 0.5.*(V(:,1:end-1) + V(:,2:end)); %C grid
% % % % Uctr = 0.5.*(U(1:end-1,:) + U(2:end,:));
% % % disp(max(U(:)));
% % % quiver(x,y,Uctr',Vctr')

%running and making output plots
[U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(psi_initial,nx,dx,ny,dy,50,dt,...
    theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,planetype,'wall','c');

[U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(psi_initial,nx,dx,ny,dy,50,dt,...
    theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,'f','wall','c');

makemodelplots('temp',psi,x/1000,y/1000,t/(24*3600),nx/2,ny/2,'\Psi (m^2/s)','%04.2f','days','m',true)
%% Racing Rossby waves



%% Midlatitude Kelvin waves



%% Basin spinup



%% Equatorial Kelvin + Rossby (delayed oscillator)



%% Equatorial- noisy forcing


