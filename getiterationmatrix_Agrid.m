%generate C matrix to iterate vaariables towards steady state solution

function C = getiterationmatrix_Agrid(~,V,Ftau,beta,kx,ky,nx,dx,ny,dy,~)
%   o C is size(nx,ny,6), with iteration values for each term in iteration
%       equation: 0 = c1Xij + c2Xi+1j + c3Xij+1 + c4Xi-1j + c5Xij-1 + F,
%       where c1-5 are stored in C(i,j,1:5) and F is stored as C(i,j,6)


%% preallocating
C = zeros(nx,ny,6);


%% looping for all points 
for i = 1:nx    
    for j = 1:ny
        C(i,j,1) = - 2*kx/dx^2 - 2*ky/dy^2;
        C(i,j,2) = + kx/dx^2;
        C(i,j,3) = + ky/dy^2;
        C(i,j,4) = + kx/dx^2;
        C(i,j,5) = + ky/dy^2;
        C(i,j,6) = Ftau(i,j) - beta(j)*V(i,j); %interpolated to grid centers
    end
end


end %end of function