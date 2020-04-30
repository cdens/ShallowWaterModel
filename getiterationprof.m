function P_iterate = getiterationprof(P,nx,ny)

%   o P is size(nx,ny), and P_iterate adds boundary conditions on y axis (x
%       axis is periodic) so is size(nx,ny+2)

%preallocation, top + bottom rows = zeros
P_iterate = zeros(nx,ny+2);

%interior
P_iterate(:,2:end-1) = P;
