%generate C matrix to iterate vaariables towards steady state solution

function C = getiterationmatrix_Dgrid(U,V,Ftau,beta,kx,ky,nx,dx,ny,dy,xbound)
%   o C is size(nx,ny,6), with iteration values for each term in iteration
%       equation: 0 = c1Xij + c2Xi+1j + c3Xij+1 + c4Xi-1j + c5Xij-1 + F,
%       where c1-5 are stored in C(i,j,1:5) and F is stored as C(i,j,6)


%% preallocating
C = zeros(nx,ny,6);


%% identifying U+/U-, V+/V- for 1st order upstream spatial discretization

Up = (U + abs(U))./2;
Um = (U - abs(U))./2;

Vp = (V + abs(V))./2;
Vm = (V - abs(V))./2;

%getting grid-centered U+, U-, V+, V-
Upc = 0.5.*(Up(:,1:end-1) + Up(:,2:end));
Umc = 0.5.*(Um(:,1:end-1) + Um(:,2:end));
Vpc = 0.5.*(Vp(1:end-1,:) + Vp(2:end,:));
Vmc = 0.5.*(Vm(1:end-1,:) + Vm(2:end,:));

%recentering onto C-grid faces
Um = zeros(nx+1,ny);
Up = Um;
Vm = zeros(nx,ny+1);
Vp = Vm;

Um(2:end-1,:) = 0.5.*(Umc(1:end-1,:) + Umc(2:end,:));
Up(2:end-1,:) = 0.5.*(Upc(1:end-1,:) + Upc(2:end,:));
Vm(:,2:end-1) = 0.5.*(Vmc(:,1:end-1) + Vmc(:,2:end));
Vp(:,2:end-1) = 0.5.*(Vpc(:,1:end-1) + Vpc(:,2:end));

%periodic u boundary
if strcmpi(xbound,'periodic')
    Um(1,:) = Um(end-1,:);
    Up(1,:) = Up(end-1,:);
    Um(end,:) = Um(1,:);
    Up(end,:) = Up(1,:);
end

%% looping for all points 
for i = 1:nx
    iface = i + 1;
    
    for j = 1:ny
        jface = j + 1;
        
        C(i,j,1) = (-Up(iface-1,j) + Um(iface,j))/dx + (-Vp(i,jface-1) + Vm(i,jface))/dy - 2*kx/dx^2 - 2*ky/dy^2;
        C(i,j,2) = -Um(iface,j)/dx + kx/dx^2;
        C(i,j,3) = -Vm(i,jface)/dy + ky/dy^2;
        C(i,j,4) = Up(iface-1,j)/dx + kx/dx^2;
        C(i,j,5) = Vp(i,jface-1)/dy + ky/dy^2;
        C(i,j,6) = Ftau(i,j) - beta(j)*Vp(i,jface-1) - beta(j)*Vm(i,jface);
    end
end


end %end of function