% runShallowWaterModel.m
% C.R. Densmore 28APR2020
%
% Driver function to run 2D hyperbolic shallow water model

function [U,V,psi,zeta,x,y,t,residual] = runShallowWaterModel(zeta,nx,dx,ny,dy,nt,dt,...
    theta,minresidual,maxiterations,Ftau,kx,ky,ctrlat,planetype,xbound,gridtype)


%% configuring wind stress curl flux 
%Ftau = 1/(rho_o*H) * dTau_y/dx - dTau_x/dy

if ismatrix(Ftau)
    if size(Ftau,1) == 1 && size(Ftau,2) == 1
        Ftau(1:nx,1:ny) = Ftau;
    elseif size(Ftau,1) ~= nx || size(Ftau,2) ~= ny
        error('Invalid initial Ftau matrix')
    end
    
    Ftau2D = Ftau; clear Ftau;
    for i = 1:nt
        Ftau(:,:,i) = Ftau2D;
    end
end

%% setting up model grid 

%preallocating 
switch(lower(gridtype))
    case 'c'
        Uout = NaN.*ones(nx+1,ny,nt+1);
        Vout = NaN.*ones(nx,ny+1,nt+1);
        
    case 'd'
        Uout = NaN.*ones(nx,ny+1,nt+1);
        Vout = NaN.*ones(nx+1,ny,nt+1);
        
    otherwise
        error('Invalid grid type selected')
        
end
psiout = NaN.*ones(nx,ny,nt+1);
zetaout = NaN.*ones(nx,ny,nt+1);
residual = NaN.*ones(1,nt+1);

%getting x/y/t
x = (0:nx-1)*dx;
y = (0:ny-1)*dy;
t = (0:nt)*dt;

%calculating beta 
ctry = y(round(ny/2));
switch(lower(gridtype))
    case 'c'
        yface = [y(1)-dy,y+dy];
    case 'd'
        yface = y;
end
% yface = [y(1)-dy,y+dy]; %y values on faces (matching V locations)
switch(lower(planetype))
    case 'real' %beta changing with latitude
        lat = ctrlat + km2deg((yface - ctry)/1000);
        beta = 2.*7.27E-5.*cosd(lat)./6356000;
    case 'beta'
        beta(1:length(yface)) = 2.*7.27E-5.*cosd(ctrlat)./6356000; 
    case 'f'
        beta(1:length(yface)) = 0;
    otherwise
        warning('Invalid plane type specified- reverting to beta plane')
        beta(1:length(yface)) = 2.*7.27E-5.*cosd(ctrlat).*dphidy; 
end


%% configuring initial fields for psi, U, V, zeta

% generating C matrix for psi 
C_psi = zeros(nx,ny,6);
cpall = [-2/dx^2 - 2/dy^2,  1/dx^2,  1/dy^2,  1/dx^2,  1/dy^2,  0];
for i = 1:6
    C_psi(:,:,i) = cpall(i);
end

%getting psi from zeta (initialize psi with zeros)
[psi,~] = iterate_cgm(C_psi,zeros(nx,ny),zeta,nx,ny,xbound,minresidual,maxiterations);

% getting initial U/V from psi
[U,V] = getvelocityfrompsi(psi,dx,dy,xbound,gridtype);

% recording initial states
Uout(:,:,1) = U;
Vout(:,:,1) = V;
psiout(:,:,1) = psi;
zetaout(:,:,1) = zeta;

%% time iteration loop
for i = 1:nt
    
    %calculate C matrix (eg Ckk)
    switch(lower(gridtype))
        case 'c'
            C = getiterationmatrix_Cgrid(U,V,Ftau,beta,kx,ky,nx,dx,ny,dy);
        case 'd'
            C = getiterationmatrix_Dgrid(U,V,Ftau,beta,kx,ky,nx,dx,ny,dy,xbound);
    end
    
    %calculate LHS and RHS iteration matrices
    C_L = -1.*C.*dt.*theta;
    C_R = C.*dt.*(1 - theta);
    C_L(:,:,1) = 1 + C_L(:,:,1);
    C_R(:,:,1) = 1 + C_R(:,:,1);

    %explicitly solve RHS of problem
    zeta_R = iterate2Dexplicit(zeta,nx,ny,C_R,xbound);

    %iterate voriticity equation forward in time with CGKM (hyperbolic problem)
    [zeta,curresidual] = iterate_cgm(C_L,zeta,zeta_R,nx,ny,xbound,minresidual,maxiterations);
    
    %calculate streamfunction from voriticity with CGKM (elliptic problem)
    [psi,~] = iterate_cgm(C_psi,psi,zeta,nx,ny,xbound,minresidual,maxiterations);
    
    %calculate velocities from streamfunction
    [U,V] = getvelocityfrompsi(psi,dx,dy,xbound,gridtype);
            
    %add current iteration to output variables
    Uout(:,:,i+1) = U;
    Vout(:,:,i+1) = V;
    psiout(:,:,i+1) = psi;
    zetaout(:,:,i+1) = zeta;
    residual(i) = curresidual;
    
    %print timestep number to command line
    if rem(i,10) == 0
        disp(['Timestep = ',num2str(i)])
    end
    
end


%% renaming variables for return
switch lower(gridtype)
    case 'c'
        U = 0.5.*(Uout(1:end-1,:,:) + Uout(2:end,:,:));
        V = 0.5.*(Vout(:,1:end-1,:) + Vout(:,2:end,:));
        
    case 'd'
        U = 0.5.*(Uout(:,1:end-1,:) + Uout(:,2:end,:));
        V = 0.5.*(Vout(1:end-1,:,:) + Vout(2:end,:,:));
        
end
psi = psiout;
zeta = zetaout;

end %end of model