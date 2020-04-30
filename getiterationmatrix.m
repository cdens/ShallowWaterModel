%generate C matrix to iterate vaariables towards steady state solution

function C = getiterationmatrix(U,V,Ftau,beta,kx,ky,nx,dx,ny,dy)
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
        C(i,j,6) = Ftau(i,j) - beta(jface-1)*Vp(i,jface-1) - beta(jface)*Vm(i,jface);
    end
end


end %end of function